#!/usr/bin/env python3

import argparse
import math
import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    from openfermion import (
        FermionOperator,
        QubitOperator,
        InteractionOperator,
        get_fermion_operator,
        normal_ordered,
        jordan_wigner,
        bravyi_kitaev,
    )
except ImportError:
    from openfermion.ops import FermionOperator, QubitOperator, InteractionOperator
    from openfermion.transforms import (
        get_fermion_operator,
        normal_ordered,
        jordan_wigner,
        bravyi_kitaev,
    )


@dataclass(frozen=True)
class FermionicJob:
    label: str
    terms: tuple
    coeffs: tuple
    support: frozenset[int]


@dataclass(frozen=True)
class PauliJob:
    label: str
    key: tuple
    coeff: complex


def load_ham2(path: Path):
    if path.suffix == ".npz":
        data = np.load(path, allow_pickle=False)
        constant = complex(np.asarray(data["constant"]).item())
        one_body = data["one_body"]
        two_body = data["two_body"]
        return InteractionOperator(constant, one_body, two_body)

    with open(path, "rb") as f:
        return pickle.load(f)


def fermion_key_has_duplicate_ladder(key):
    return len(set(key)) < len(key)


def clean_fermion_operator(op, tol):
    op = normal_ordered(op)
    op.compress(abs_tol=tol)

    cleaned = FermionOperator()

    for key, coeff in op.terms.items():
        if abs(coeff) <= tol:
            continue

        # Terms such as a0^ a0^ a0 a0 are algebraically zero.
        if key != () and fermion_key_has_duplicate_ladder(key):
            continue

        cleaned += FermionOperator(key, coeff)

    cleaned = normal_ordered(cleaned)
    cleaned.compress(abs_tol=tol)
    return cleaned


def dagger_key(key):
    return tuple((mode, 1 - action) for mode, action in reversed(key))


def canonical_block_key(key):
    return min(key, dagger_key(key))


def support_of_fermion_key(key):
    return frozenset(mode for mode, _action in key)


def fermion_operator_to_jobs(fermion_op, tol):
    fermion_op = clean_fermion_operator(fermion_op, tol)

    grouped = {}

    for key, coeff in fermion_op.terms.items():
        if key == ():
            continue
        if abs(coeff) <= tol:
            continue

        block_key = canonical_block_key(key)
        grouped.setdefault(block_key, {})
        grouped[block_key][key] = grouped[block_key].get(key, 0.0) + coeff

    jobs = []

    for block_key in sorted(grouped.keys()):
        op = FermionOperator()

        for key, coeff in grouped[block_key].items():
            op += FermionOperator(key, coeff)

        op = clean_fermion_operator(op, tol)

        if len(op.terms) == 0:
            continue

        keys = tuple(sorted(op.terms.keys()))
        coeffs = tuple(op.terms[k] for k in keys)

        support = frozenset()
        for k in keys:
            support = support.union(support_of_fermion_key(k))

        jobs.append(
            FermionicJob(
                label=f"F{len(jobs):04d}",
                terms=keys,
                coeffs=coeffs,
                support=support,
            )
        )

    return jobs


def job_to_fermion_operator(job, tol):
    op = FermionOperator()
    for key, coeff in zip(job.terms, job.coeffs):
        op += FermionOperator(key, coeff)
    return clean_fermion_operator(op, tol)


def exact_fermion_commutes(job_i, job_j, tol):
    op_i = job_to_fermion_operator(job_i, tol)
    op_j = job_to_fermion_operator(job_j, tol)

    comm = clean_fermion_operator(op_i * op_j - op_j * op_i, tol)
    return len(comm.terms) == 0


def clean_qubit_operator(op, tol):
    op.compress(abs_tol=tol)

    cleaned = QubitOperator()

    for key, coeff in op.terms.items():
        if abs(coeff) <= tol:
            continue
        cleaned += QubitOperator(key, coeff)

    cleaned.compress(abs_tol=tol)
    return cleaned


def qubit_operator_to_pauli_jobs(qubit_op, tol):
    qubit_op = clean_qubit_operator(qubit_op, tol)

    jobs = []

    for key in sorted(qubit_op.terms.keys()):
        if key == ():
            continue

        jobs.append(
            PauliJob(
                label=f"P{len(jobs):04d}",
                key=key,
                coeff=qubit_op.terms[key],
            )
        )

    return jobs


def pauli_strings_commute(key_a, key_b):
    dict_a = dict(key_a)
    dict_b = dict(key_b)

    overlap = set(dict_a).intersection(dict_b)

    anti_count = 0

    for q in overlap:
        if dict_a[q] != dict_b[q]:
            anti_count += 1

    return anti_count % 2 == 0


def build_commutation_matrix(jobs, commute_func):
    n = len(jobs)
    commutes = [[False for _ in range(n)] for _ in range(n)]

    for i in range(n):
        commutes[i][i] = True

    for i in range(n):
        for j in range(i + 1, n):
            val = commute_func(jobs[i], jobs[j])
            commutes[i][j] = val
            commutes[j][i] = val

    return commutes


def noncommuting_edges(commutes):
    n = len(commutes)
    return [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if not commutes[i][j]
    ]


def graph_components(n, edges):
    adj = {i: set() for i in range(n)}

    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)

    seen = set()
    components = []

    for start in range(n):
        if start in seen:
            continue

        stack = [start]
        seen.add(start)
        comp = []

        while stack:
            v = stack.pop()
            comp.append(v)

            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)

        components.append(sorted(comp))

    return components


def count_acyclic_orientations_bruteforce(vertices, edges):
    """
    Count acyclic orientations by trying all 2^e edge orientations.

    This is fine for small edge counts, e.g. Pauli H2.
    """
    vertex_map = {v: k for k, v in enumerate(vertices)}
    local_edges = [(vertex_map[i], vertex_map[j]) for i, j in edges]

    n = len(vertices)
    e = len(local_edges)

    count = 0

    for mask in range(1 << e):
        adj = [[] for _ in range(n)]
        indeg = [0 for _ in range(n)]

        for bit, (u, v) in enumerate(local_edges):
            if (mask >> bit) & 1:
                a, b = u, v
            else:
                a, b = v, u

            adj[a].append(b)
            indeg[b] += 1

        queue = [i for i in range(n) if indeg[i] == 0]
        seen_count = 0

        while queue:
            x = queue.pop()
            seen_count += 1

            for y in adj[x]:
                indeg[y] -= 1
                if indeg[y] == 0:
                    queue.append(y)

        if seen_count == n:
            count += 1

    return count


def summarize_graph(name, jobs, commutes, orientation_limit):
    n = len(jobs)
    total_pairs = n * (n - 1) // 2

    edges = noncommuting_edges(commutes)
    e = len(edges)
    commuting = total_pairs - e

    print()
    print("=" * 80)
    print(name)
    print("=" * 80)
    print(f"Jobs/terms being ordered: {n}")
    print(f"Raw orders: {n}! = {math.factorial(n)}")
    print(f"Pairwise checks: C({n},2) = {total_pairs}")
    print(f"Commuting pairs: {commuting}")
    print(f"Noncommuting pairs: {e}")

    if total_pairs:
        print(f"Commutation density: {commuting / total_pairs:.6f}")
        print(f"Noncommutation density: {e / total_pairs:.6f}")

    print()
    print(f"Trace-class upper bound: 2^{e} = {2 ** e}")

    comps = graph_components(n, edges)

    print()
    print("Noncommutation graph components:")
    exact_total = 1
    exact_known = True

    for k, comp in enumerate(comps):
        comp_set = set(comp)
        comp_edges = [
            (i, j)
            for i, j in edges
            if i in comp_set and j in comp_set
        ]

        labels = [jobs[i].label for i in comp]

        print(f"  component {k}: vertices={len(comp)}, edges={len(comp_edges)}, labels={labels}")

        if len(comp_edges) <= orientation_limit:
            c = count_acyclic_orientations_bruteforce(comp, comp_edges)
            exact_total *= c
            print(f"    exact acyclic orientations: {c}")
        else:
            exact_known = False
            print(f"    exact acyclic orientations: skipped because edges > {orientation_limit}")

    print()

    if exact_known:
        print(f"Exact trace classes: {exact_total}")
        print(f"Trotter-error calculations needed if exhaustive: {exact_total}")
        print(f"Reduction factor from raw orders: {math.factorial(n) / exact_total:.6e}")
    else:
        print("Exact trace classes: not computed")
        print(f"Trotter-error calculations needed if exhaustive: <= {2 ** e}")
        print(f"Guaranteed reduction factor from raw orders: {math.factorial(n) / (2 ** e):.6e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ham2_file")
    parser.add_argument("--mapping", choices=["jw", "bk"], default="jw")
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument(
        "--orientation-limit",
        type=int,
        default=24,
        help="Brute-force exact acyclic orientation count only for components with <= this many edges.",
    )
    args = parser.parse_args()

    ham2 = load_ham2(Path(args.ham2_file))

    fermion_op = get_fermion_operator(ham2)
    fermion_op = clean_fermion_operator(fermion_op, args.tol)

    print()
    print(f"Loaded ham2 file: {args.ham2_file}")
    print(f"Cleaned FermionOperator terms including constant: {len(fermion_op.terms)}")
    print(f"Cleaned nonconstant FermionOperator monomials: {sum(1 for k in fermion_op.terms if k != ())}")

    fermion_jobs = fermion_operator_to_jobs(fermion_op, args.tol)

    fermion_commutes = build_commutation_matrix(
        fermion_jobs,
        lambda a, b: exact_fermion_commutes(a, b, args.tol),
    )

    if args.mapping == "jw":
        qubit_op = jordan_wigner(fermion_op)
        mapping_name = "JW Pauli-string ordering"
    else:
        qubit_op = bravyi_kitaev(fermion_op)
        mapping_name = "BK Pauli-string ordering"

    qubit_op = clean_qubit_operator(qubit_op, args.tol)
    pauli_jobs = qubit_operator_to_pauli_jobs(qubit_op, args.tol)

    pauli_commutes = build_commutation_matrix(
        pauli_jobs,
        lambda a, b: pauli_strings_commute(a.key, b.key),
    )

    summarize_graph(
        "Fermionic exact-commutation ordering",
        fermion_jobs,
        fermion_commutes,
        args.orientation_limit,
    )

    summarize_graph(
        mapping_name,
        pauli_jobs,
        pauli_commutes,
        args.orientation_limit,
    )


if __name__ == "__main__":
    main()
