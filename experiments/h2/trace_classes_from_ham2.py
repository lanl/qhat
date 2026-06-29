#!/usr/bin/env python3

import argparse
import math
import pickle
import random
from collections import Counter
from dataclasses import dataclass
from heapq import heappop, heappush
from itertools import permutations
from pathlib import Path

import numpy as np

try:
    from openfermion import (
        FermionOperator,
        InteractionOperator,
        get_fermion_operator,
        normal_ordered,
    )
except ImportError:
    from openfermion.ops import FermionOperator, InteractionOperator
    from openfermion.transforms import get_fermion_operator, normal_ordered


@dataclass(frozen=True)
class FermionicJob:
    label: str
    terms: tuple
    coeffs: tuple
    support: frozenset[int]


def load_ham2(path: Path):
    """
    Load QHAT ham2 active-space data.

    Supports:
      1. ham2 pickle:       *_as-XXX-YYY.pickle
      2. ham2 tensor npz:   *_as-XXX-YYY.tensors.npz

    The npz file is safer if pickle loading fails due Python/OpenFermion version issues.
    """
    if path.suffix == ".npz":
        data = np.load(path, allow_pickle=False)

        print(f"Loaded NPZ file: {path}")
        print(f"NPZ keys: {list(data.keys())}")

        constant = complex(np.asarray(data["constant"]).item())
        one_body = data["one_body"]
        two_body = data["two_body"]

        print(f"constant: {constant}")
        print(f"one_body shape: {one_body.shape}")
        print(f"two_body shape: {two_body.shape}")

        return InteractionOperator(constant, one_body, two_body)

    with open(path, "rb") as f:
        obj = pickle.load(f)

    print(f"Loaded pickle file: {path}")
    print(f"Loaded object type: {type(obj)}")

    if hasattr(obj, "constant"):
        print(f"constant: {obj.constant}")
    if hasattr(obj, "one_body_tensor"):
        print(f"one_body_tensor shape: {obj.one_body_tensor.shape}")
    if hasattr(obj, "two_body_tensor"):
        print(f"two_body_tensor shape: {obj.two_body_tensor.shape}")
    if hasattr(obj, "n_qubits"):
        print(f"n_qubits / spin orbitals: {obj.n_qubits}")

    return obj


def dagger_key(key):
    """
    Hermitian conjugate of an OpenFermion FermionOperator key.

    Example:
        ((0, 1), (1, 0)) -> ((1, 1), (0, 0))

    action = 1 means creation.
    action = 0 means annihilation.
    """
    return tuple((mode, 1 - action) for mode, action in reversed(key))


def canonical_block_key(key):
    """
    Use the lexicographically smaller of a term and its Hermitian conjugate
    to group them into one Hermitian block/job.
    """
    return min(key, dagger_key(key))


def support_of_key(key):
    return frozenset(mode for mode, _action in key)


def format_key(key):
    if key == ():
        return "I"

    pieces = []
    for mode, action in key:
        if action == 1:
            pieces.append(f"a{mode}^")
        else:
            pieces.append(f"a{mode}")
    return " ".join(pieces)


def interaction_to_fermionic_jobs(ham2, tol=1e-12):
    """
    Convert an OpenFermion InteractionOperator into Hermitian fermionic jobs.

    We first convert ham2 to a FermionOperator. Then we group each monomial
    with its Hermitian conjugate.

    Constant terms are skipped because they commute with everything and do not
    affect ordering.
    """
    fermion_op = get_fermion_operator(ham2)

    raw_terms = {
        key: coeff
        for key, coeff in fermion_op.terms.items()
        if key != () and abs(coeff) > tol
    }

    print()
    print(f"Raw nonconstant FermionOperator monomials: {len(raw_terms)}")

    grouped = {}

    for key, coeff in raw_terms.items():
        block_key = canonical_block_key(key)

        if block_key not in grouped:
            grouped[block_key] = {}

        grouped[block_key][key] = grouped[block_key].get(key, 0.0) + coeff

    jobs = []

    for idx, block_key in enumerate(sorted(grouped.keys())):
        term_dict = grouped[block_key]

        keys = tuple(sorted(term_dict.keys()))
        coeffs = tuple(term_dict[k] for k in keys)

        support = frozenset()
        for k in keys:
            support = support.union(support_of_key(k))

        jobs.append(
            FermionicJob(
                label=f"F{idx:04d}",
                terms=keys,
                coeffs=coeffs,
                support=support,
            )
        )

    return jobs


def job_to_fermion_operator(job: FermionicJob):
    op = FermionOperator()
    for key, coeff in zip(job.terms, job.coeffs):
        op += FermionOperator(key, coeff)
    return op


def support_commutes(job_i: FermionicJob, job_j: FermionicJob):
    """
    Conservative sufficient condition:
    disjoint-support even fermionic blocks commute.

    This can miss some commuting pairs, but any pair marked commuting here is safe.
    """
    return job_i.support.isdisjoint(job_j.support)


def exact_commutes(job_i: FermionicJob, job_j: FermionicJob, tol=1e-10):
    """
    Exact symbolic fermionic commutator check using OpenFermion algebra.

    This is stronger than the support-disjoint rule.
    """
    op_i = job_to_fermion_operator(job_i)
    op_j = job_to_fermion_operator(job_j)

    comm = normal_ordered(op_i * op_j - op_j * op_i)

    for coeff in comm.terms.values():
        if abs(coeff) > tol:
            return False

    return True


def build_commutation_matrix(jobs, mode="support", tol=1e-10):
    n = len(jobs)
    commutes = [[False for _ in range(n)] for _ in range(n)]

    for i in range(n):
        commutes[i][i] = True

    for i in range(n):
        for j in range(i + 1, n):
            if mode == "support":
                val = support_commutes(jobs[i], jobs[j])
            elif mode == "exact":
                val = exact_commutes(jobs[i], jobs[j], tol=tol)
            else:
                raise ValueError(f"Unknown commutation mode: {mode}")

            commutes[i][j] = val
            commutes[j][i] = val

    return commutes


def canonical_order(order, commutes):
    """
    Build a DAG preserving the relative order of every noncommuting pair.
    Then return the lexicographically smallest topological sort.
    """
    out_edges = {x: set() for x in order}
    indeg = {x: 0 for x in order}

    for a, i in enumerate(order):
        for j in order[a + 1:]:
            if not commutes[i][j]:
                if j not in out_edges[i]:
                    out_edges[i].add(j)
                    indeg[j] += 1

    heap = []
    for x in sorted(order):
        if indeg[x] == 0:
            heappush(heap, x)

    canon = []

    while heap:
        x = heappop(heap)
        canon.append(x)

        for y in sorted(out_edges[x]):
            indeg[y] -= 1
            if indeg[y] == 0:
                heappush(heap, y)

    if len(canon) != len(order):
        raise RuntimeError("Cycle detected. This should not happen for an order-induced DAG.")

    return tuple(canon)


def exhaustive_classes(n, commutes):
    classes = Counter()

    for pi in permutations(range(n)):
        canon = canonical_order(pi, commutes)
        classes[canon] += 1

    return classes


def sampled_classes(n, commutes, samples, seed):
    rng = random.Random(seed)
    classes = Counter()
    base = list(range(n))

    for _ in range(samples):
        order = base[:]
        rng.shuffle(order)
        canon = canonical_order(tuple(order), commutes)
        classes[canon] += 1

    return classes


def print_jobs(jobs, max_jobs=20):
    print()
    print(f"Number of Hermitian fermionic jobs/blocks: {len(jobs)}")
    print()

    for job in jobs[:max_jobs]:
        print(f"{job.label}: support={sorted(job.support)}")
        for coeff, key in zip(job.coeffs, job.terms):
            print(f"    {coeff:+.12e}   {format_key(key)}")
        print()

    if len(jobs) > max_jobs:
        print(f"... skipped {len(jobs) - max_jobs} additional jobs")
        print()


def print_commutation_summary(jobs, commutes, mode):
    n = len(jobs)
    total_pairs = n * (n - 1) // 2

    commuting_pairs = sum(
        1
        for i in range(n)
        for j in range(i + 1, n)
        if commutes[i][j]
    )

    noncommuting_pairs = total_pairs - commuting_pairs

    print()
    print("=== Commutation summary ===")
    print(f"Commutation mode: {mode}")
    print(f"Number of jobs: {n}")
    print(f"Total pairs: {total_pairs}")
    print(f"Commuting pairs: {commuting_pairs}")
    print(f"Noncommuting pairs: {noncommuting_pairs}")

    if total_pairs > 0:
        print(f"Commutation density: {commuting_pairs / total_pairs:.6f}")
        print(f"Noncommutation density: {noncommuting_pairs / total_pairs:.6f}")


def print_class_summary(jobs, classes, exhaustive, samples):
    n = len(jobs)

    print()
    print("=== Trace-class summary ===")

    if exhaustive:
        raw = math.factorial(n)
        print(f"Raw permutations: {raw}")
        print(f"Exact trace classes: {len(classes)}")
        print(f"Exact compression ratio: {raw / len(classes):.6f}")
    else:
        print(f"Sampled raw permutations: {samples}")
        print(f"Observed trace classes: {len(classes)}")
        print(f"Sample compression ratio: {samples / len(classes):.6f}")

    print()
    print("Largest trace classes:")
    for canon, count in classes.most_common(10):
        label_order = tuple(jobs[i].label for i in canon)
        print(f"  size={count:8d}  canon={label_order}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ham2_file",
        help="QHAT ham2 pickle or ham2 tensors npz",
    )
    parser.add_argument(
        "--commutation",
        choices=["support", "exact"],
        default="support",
        help="support = disjoint-support sufficient rule; exact = symbolic fermion commutator",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=10000,
        help="Number of random raw orders to sample when exhaustive enumeration is too large.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--exhaustive-limit",
        type=int,
        default=8,
        help="Enumerate all permutations only if number of jobs <= this limit.",
    )
    parser.add_argument(
        "--max-print-jobs",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-10,
    )

    args = parser.parse_args()

    ham2_path = Path(args.ham2_file)

    ham2 = load_ham2(ham2_path)
    jobs = interaction_to_fermionic_jobs(ham2)

    print_jobs(jobs, max_jobs=args.max_print_jobs)

    commutes = build_commutation_matrix(
        jobs,
        mode=args.commutation,
        tol=args.tol,
    )

    print_commutation_summary(jobs, commutes, mode=args.commutation)

    n = len(jobs)

    if n <= args.exhaustive_limit:
        classes = exhaustive_classes(n, commutes)
        exhaustive = True
    else:
        classes = sampled_classes(
            n=n,
            commutes=commutes,
            samples=args.samples,
            seed=args.seed,
        )
        exhaustive = False

    print_class_summary(
        jobs=jobs,
        classes=classes,
        exhaustive=exhaustive,
        samples=args.samples,
    )

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
    comps = []

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

        comps.append(sorted(comp))

    return comps


def print_noncommutation_graph_summary(jobs, commutes):
    edges = noncommuting_edges(commutes)
    comps = graph_components(len(jobs), edges)

    print()
    print("=== Noncommutation graph summary ===")
    print(f"Vertices/jobs: {len(jobs)}")
    print(f"Noncommuting edges: {len(edges)}")
    print(f"Trivial upper bound on trace classes: 2^{len(edges)} = {2 ** len(edges)}")
    print()

    for k, comp in enumerate(comps):
        comp_set = set(comp)
        comp_edges = [
            (i, j)
            for i, j in edges
            if i in comp_set and j in comp_set
        ]

        labels = [jobs[i].label for i in comp]

        print(f"Component {k}:")
        print(f"  vertices: {len(comp)}")
        print(f"  edges: {len(comp_edges)}")
        print(f"  labels: {labels}")
        print()

if __name__ == "__main__":
    main()