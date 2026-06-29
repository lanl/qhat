#!/usr/bin/env python3

import argparse
import csv
import json
import pickle
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from heapq import heappop, heappush
from pathlib import Path

import numpy as np
from scipy.linalg import expm, svdvals

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
class PauliJob:
    label: str
    key: tuple
    coeff: complex


class Tee:
    """
    Write stdout/stderr both to terminal and to a log file.
    """
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def setup_logging(outdir: Path, log_file: str | None):
    outdir.mkdir(parents=True, exist_ok=True)

    if log_file is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = outdir / f"scan_h2_pauli_trotter_errors_{stamp}.log"
    else:
        log_path = Path(log_file)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = open(log_path, "w", buffering=1)

    # This is the line that saves all printed log output.
    sys.stdout = Tee(sys.__stdout__, log_fh)
    sys.stderr = Tee(sys.__stderr__, log_fh)

    return log_path, log_fh


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

        if key != () and fermion_key_has_duplicate_ladder(key):
            continue

        cleaned += FermionOperator(key, coeff)

    cleaned = normal_ordered(cleaned)
    cleaned.compress(abs_tol=tol)
    return cleaned


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

        coeff = qubit_op.terms[key]

        jobs.append(
            PauliJob(
                label=f"P{len(jobs):04d}",
                key=key,
                coeff=coeff,
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


def build_commutation_matrix(jobs):
    n = len(jobs)
    commutes = [[False for _ in range(n)] for _ in range(n)]

    for i in range(n):
        commutes[i][i] = True

    for i in range(n):
        for j in range(i + 1, n):
            val = pauli_strings_commute(jobs[i].key, jobs[j].key)
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


def canonical_order_from_raw_order(order, commutes):
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
        raise RuntimeError("Cycle detected in order-induced DAG.")

    return tuple(canon)


def lexicographic_topological_order(n, directed_edges):
    out_edges = {i: set() for i in range(n)}
    indeg = {i: 0 for i in range(n)}

    for i, j in directed_edges:
        if j not in out_edges[i]:
            out_edges[i].add(j)
            indeg[j] += 1

    heap = []
    for i in range(n):
        if indeg[i] == 0:
            heappush(heap, i)

    order = []

    while heap:
        x = heappop(heap)
        order.append(x)

        for y in sorted(out_edges[x]):
            indeg[y] -= 1
            if indeg[y] == 0:
                heappush(heap, y)

    if len(order) != n:
        return None

    return tuple(order)


def enumerate_trace_class_representatives(n, edges, max_orientation_edges):
    e = len(edges)

    if e > max_orientation_edges:
        raise ValueError(
            f"Cannot brute-force 2^{e} edge orientations. "
            f"Increase --max-orientation-edges if this is intentional."
        )

    reps = set()
    total_masks = 1 << e

    print()
    print(f"Enumerating acyclic orientations over {e} noncommuting edges.")
    print(f"Orientation masks to check: {total_masks}")

    for mask in range(total_masks):
        if mask % 10000 == 0:
            print(f"  checked {mask}/{total_masks} orientations")

        directed = []

        for bit, (i, j) in enumerate(edges):
            if (mask >> bit) & 1:
                directed.append((i, j))
            else:
                directed.append((j, i))

        topo = lexicographic_topological_order(n, directed)

        if topo is not None:
            reps.add(topo)

    reps = sorted(reps)

    print(f"Trace-class representatives found: {len(reps)}")

    return reps


def pauli_key_to_string(key):
    if key == ():
        return "I"
    return " ".join(f"{p}{q}" for q, p in key)


def kron_all(mats):
    out = mats[0]
    for mat in mats[1:]:
        out = np.kron(out, mat)
    return out


def pauli_matrix_from_key(key, n_qubits):
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    local = {
        "I": I,
        "X": X,
        "Y": Y,
        "Z": Z,
    }

    ops = dict(key)

    mats = []
    for q in range(n_qubits):
        mats.append(local[ops.get(q, "I")])

    return kron_all(mats)


def build_hamiltonian_and_term_unitaries(jobs, n_qubits, dt, tol):
    dim = 2 ** n_qubits
    I_dim = np.eye(dim, dtype=complex)

    H = np.zeros((dim, dim), dtype=complex)
    term_unitaries = []

    for job in jobs:
        P = pauli_matrix_from_key(job.key, n_qubits)
        coeff = complex(job.coeff)

        if abs(coeff.imag) > tol:
            print(f"WARNING: non-negligible imaginary coefficient for {job.label}: {coeff}")

        H += coeff * P

        if abs(coeff.imag) <= tol:
            theta = dt * coeff.real
            U_term = np.cos(theta) * I_dim - 1j * np.sin(theta) * P
        else:
            U_term = expm(-1j * dt * coeff * P)

        term_unitaries.append(U_term)

    H = 0.5 * (H + H.conj().T)
    U_exact = expm(-1j * dt * H)

    return H, U_exact, term_unitaries


def trotter_unitary(order, term_unitaries):
    dim = term_unitaries[0].shape[0]
    U = np.eye(dim, dtype=complex)

    for idx in order:
        U = U @ term_unitaries[idx]

    return U


def occupation_basis_state(n_qubits, n_electrons):
    """
    Computational basis state with the first n_electrons qubits occupied.

    With our Kronecker convention, qubit 0 is the leftmost tensor factor.
    """
    dim = 2 ** n_qubits
    bits = [1 if q < n_electrons else 0 for q in range(n_qubits)]

    basis_index = 0
    for bit in bits:
        basis_index = 2 * basis_index + bit

    psi = np.zeros(dim, dtype=complex)
    psi[basis_index] = 1.0

    bitstring = "".join(str(b) for b in bits)

    return psi, bitstring


def compute_errors(order, U_exact, term_unitaries, psi):
    U_trot = trotter_unitary(order, term_unitaries)
    diff = U_exact - U_trot

    operator_error = float(svdvals(diff)[0])

    if psi is None:
        state_error = None
    else:
        state_error = float(np.linalg.norm(diff @ psi))

    return operator_error, state_error


def safe_float_tag(x):
    return str(x).replace(".", "p").replace("-", "m")


def scan_mapping(
    mapping,
    fermion_op,
    n_qubits,
    n_electrons,
    dt,
    tol,
    outdir,
    max_orientation_edges,
):
    print()
    print("=" * 80)
    print(f"Mapping: {mapping.upper()}")
    print("=" * 80)

    if mapping == "jw":
        qubit_op = jordan_wigner(fermion_op)
    elif mapping == "bk":
        try:
            qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_qubits)
        except TypeError:
            qubit_op = bravyi_kitaev(fermion_op)
    else:
        raise ValueError(f"Unknown mapping: {mapping}")

    qubit_op = clean_qubit_operator(qubit_op, tol)

    identity_coeff = qubit_op.terms.get((), 0.0)
    jobs = qubit_operator_to_pauli_jobs(qubit_op, tol)

    print(f"Identity coefficient ignored for ordering: {identity_coeff}")
    print(f"Number of non-identity Pauli jobs: {len(jobs)}")
    print("Pauli jobs:")
    for job in jobs:
        print(f"  {job.label}: coeff={job.coeff:+.12e}, string={pauli_key_to_string(job.key)}")

    commutes = build_commutation_matrix(jobs)
    edges = noncommuting_edges(commutes)

    n = len(jobs)
    total_pairs = n * (n - 1) // 2
    commuting_pairs = total_pairs - len(edges)

    print()
    print("Commutation graph:")
    print(f"  pairwise checks: {total_pairs}")
    print(f"  commuting pairs: {commuting_pairs}")
    print(f"  noncommuting pairs: {len(edges)}")

    reps = enumerate_trace_class_representatives(
        n=n,
        edges=edges,
        max_orientation_edges=max_orientation_edges,
    )

    print()
    print("Building matrices and term exponentials.")
    H, U_exact, term_unitaries = build_hamiltonian_and_term_unitaries(
        jobs=jobs,
        n_qubits=n_qubits,
        dt=dt,
        tol=tol,
    )

    psi, bitstring = occupation_basis_state(n_qubits, n_electrons)
    print(f"State-dependent error uses occupation bitstring: {bitstring}")

    default_order = tuple(range(n))
    reverse_order = tuple(reversed(range(n)))
    coeff_desc_order = tuple(sorted(range(n), key=lambda i: (-abs(jobs[i].coeff), i)))
    coeff_asc_order = tuple(sorted(range(n), key=lambda i: (abs(jobs[i].coeff), i)))

    heuristic_reps = {
        "default": canonical_order_from_raw_order(default_order, commutes),
        "reverse": canonical_order_from_raw_order(reverse_order, commutes),
        "coeff_desc": canonical_order_from_raw_order(coeff_desc_order, commutes),
        "coeff_asc": canonical_order_from_raw_order(coeff_asc_order, commutes),
    }

    print()
    print("Heuristic class representatives:")
    for name, rep in heuristic_reps.items():
        print(f"  {name}: {rep}")

    rows = []

    print()
    print(f"Computing Trotter errors for {len(reps)} trace-class representatives.")

    start = time.time()

    for class_id, order in enumerate(reps):
        if class_id % 500 == 0:
            print(f"  computed {class_id}/{len(reps)} classes")

        operator_error, state_error = compute_errors(
            order=order,
            U_exact=U_exact,
            term_unitaries=term_unitaries,
            psi=psi,
        )

        tags = [
            name
            for name, rep in heuristic_reps.items()
            if rep == order
        ]

        rows.append(
            {
                "mapping": mapping,
                "class_id": class_id,
                "operator_error": operator_error,
                "state_error": state_error,
                "tags": ";".join(tags),
                "order_indices": " ".join(str(i) for i in order),
                "order_labels": " ".join(jobs[i].label for i in order),
            }
        )

    elapsed = time.time() - start
    print(f"Finished error scan in {elapsed:.3f} seconds.")

    rows_sorted = sorted(rows, key=lambda r: r["operator_error"])

    for rank, row in enumerate(rows_sorted, start=1):
        row["operator_error_rank"] = rank

    dt_tag = safe_float_tag(dt)
    csv_path = outdir / f"h2_pauli_trotter_errors_{mapping}_dt{dt_tag}.csv"
    json_path = outdir / f"h2_pauli_trotter_errors_{mapping}_dt{dt_tag}_summary.json"

    fieldnames = [
        "mapping",
        "class_id",
        "operator_error_rank",
        "operator_error",
        "state_error",
        "tags",
        "order_indices",
        "order_labels",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_sorted)

    op_errors = np.array([r["operator_error"] for r in rows], dtype=float)
    state_errors = np.array([r["state_error"] for r in rows], dtype=float)

    summary = {
        "mapping": mapping,
        "n_qubits": n_qubits,
        "n_electrons": n_electrons,
        "dt": dt,
        "n_pauli_jobs": n,
        "pairwise_checks": total_pairs,
        "commuting_pairs": commuting_pairs,
        "noncommuting_pairs": len(edges),
        "trace_classes": len(reps),
        "csv_path": str(csv_path),
        "operator_error": {
            "min": float(np.min(op_errors)),
            "median": float(np.median(op_errors)),
            "mean": float(np.mean(op_errors)),
            "max": float(np.max(op_errors)),
        },
        "state_error": {
            "min": float(np.min(state_errors)),
            "median": float(np.median(state_errors)),
            "mean": float(np.mean(state_errors)),
            "max": float(np.max(state_errors)),
        },
        "best_order": rows_sorted[0],
        "worst_order": rows_sorted[-1],
        "heuristic_classes": {},
    }

    for name, rep in heuristic_reps.items():
        match = None
        for row in rows:
            if tuple(int(x) for x in row["order_indices"].split()) == rep:
                match = row
                break

        if match is not None:
            summary["heuristic_classes"][name] = {
                "class_id": match["class_id"],
                "operator_error": match["operator_error"],
                "state_error": match["state_error"],
                "order_indices": match["order_indices"],
                "order_labels": match["order_labels"],
            }

    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print(f"Saved CSV results: {csv_path}")
    print(f"Saved JSON summary: {json_path}")

    print()
    print("Summary:")
    print(f"  trace classes: {len(reps)}")
    print(f"  operator error min:    {summary['operator_error']['min']:.12e}")
    print(f"  operator error median: {summary['operator_error']['median']:.12e}")
    print(f"  operator error max:    {summary['operator_error']['max']:.12e}")
    print(f"  state error min:       {summary['state_error']['min']:.12e}")
    print(f"  state error median:    {summary['state_error']['median']:.12e}")
    print(f"  state error max:       {summary['state_error']['max']:.12e}")

    print()
    print("Heuristic errors:")
    for name, data in summary["heuristic_classes"].items():
        print(
            f"  {name}: "
            f"class={data['class_id']}, "
            f"operator_error={data['operator_error']:.12e}, "
            f"state_error={data['state_error']:.12e}"
        )

    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ham2_file")
    parser.add_argument("--mapping", choices=["jw", "bk", "both"], default="both")
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--n-electrons", type=int, default=2)
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--max-orientation-edges", type=int, default=24)
    parser.add_argument(
        "--outdir",
        default="experiments/h2/trotter_error_scan",
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Optional explicit log file. If omitted, a timestamped log is created in outdir.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    log_path, log_fh = setup_logging(outdir, args.log_file)

    try:
        print(f"Log file: {log_path}")
        print(f"Started: {datetime.now().isoformat(timespec='seconds')}")
        print(f"Arguments: {vars(args)}")

        ham2 = load_ham2(Path(args.ham2_file))

        n_qubits = getattr(ham2, "n_qubits", None)
        if n_qubits is None:
            raise ValueError("Could not infer n_qubits from ham2 object.")

        fermion_op = get_fermion_operator(ham2)
        fermion_op = clean_fermion_operator(fermion_op, args.tol)

        print()
        print(f"Loaded ham2 file: {args.ham2_file}")
        print(f"n_qubits / spin orbitals: {n_qubits}")
        print(f"Cleaned FermionOperator terms including constant: {len(fermion_op.terms)}")
        print(
            "Cleaned nonconstant FermionOperator monomials: "
            f"{sum(1 for k in fermion_op.terms if k != ())}"
        )

        if args.mapping == "both":
            mappings = ["jw", "bk"]
        else:
            mappings = [args.mapping]

        summaries = []

        for mapping in mappings:
            summary = scan_mapping(
                mapping=mapping,
                fermion_op=fermion_op,
                n_qubits=n_qubits,
                n_electrons=args.n_electrons,
                dt=args.dt,
                tol=args.tol,
                outdir=outdir,
                max_orientation_edges=args.max_orientation_edges,
            )
            summaries.append(summary)

        combined_path = outdir / f"h2_pauli_trotter_errors_dt{safe_float_tag(args.dt)}_combined_summary.json"

        with open(combined_path, "w") as f:
            json.dump(summaries, f, indent=2)

        print()
        print(f"Saved combined summary: {combined_path}")
        print(f"Finished: {datetime.now().isoformat(timespec='seconds')}")

    finally:
        log_fh.flush()
        log_fh.close()


if __name__ == "__main__":
    main()
