#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.linalg import svdvals
from scipy.stats import spearmanr, kendalltau


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def setup_logging(outdir, log_file=None):
    outdir.mkdir(parents=True, exist_ok=True)

    if log_file is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = outdir / f"evaluate_commutator_surrogate_{stamp}.log"
    else:
        log_path = Path(log_file)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_fh = open(log_path, "w", buffering=1)

    orig_stdout = sys.stdout
    orig_stderr = sys.stderr

    sys.stdout = Tee(sys.stdout, log_fh)
    sys.stderr = Tee(sys.stderr, log_fh)

    return log_path, log_fh, orig_stdout, orig_stderr


def import_scan_module(path):
    path = Path(path)

    spec = importlib.util.spec_from_file_location("scan_h2_pauli_trotter_errors_imported", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import scan module from {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    return module


def safe_float_tag(x):
    return str(x).replace(".", "p").replace("-", "m")


def read_rows(csv_path):
    rows = []

    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            row["class_id"] = int(row["class_id"])
            row["operator_error_rank"] = int(row["operator_error_rank"])
            row["operator_error"] = float(row["operator_error"])
            row["state_error"] = float(row["state_error"])
            row["order_tuple"] = tuple(int(x) for x in row["order_indices"].split())
            rows.append(row)

    return rows


def write_rows(csv_path, rows):
    fieldnames = [
        "mapping",
        "class_id",
        "operator_error_rank",
        "operator_error",
        "state_error",
        "commutator_surrogate",
        "predicted_leading_error",
        "surrogate_rank",
        "tags",
        "order_indices",
        "order_labels",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            out = {k: row[k] for k in fieldnames}
            writer.writerow(out)


def build_pauli_jobs(scan, ham2_file, mapping, tol):
    ham2 = scan.load_ham2(Path(ham2_file))

    n_qubits = getattr(ham2, "n_qubits", None)
    if n_qubits is None:
        raise ValueError("Could not infer n_qubits from ham2 object.")

    fermion_op = scan.get_fermion_operator(ham2)
    fermion_op = scan.clean_fermion_operator(fermion_op, tol)

    if mapping == "jw":
        qubit_op = scan.jordan_wigner(fermion_op)
    elif mapping == "bk":
        try:
            qubit_op = scan.bravyi_kitaev(fermion_op, n_qubits=n_qubits)
        except TypeError:
            qubit_op = scan.bravyi_kitaev(fermion_op)
    else:
        raise ValueError(f"Unknown mapping: {mapping}")

    qubit_op = scan.clean_qubit_operator(qubit_op, tol)
    jobs = scan.qubit_operator_to_pauli_jobs(qubit_op, tol)

    return jobs, n_qubits


def build_term_matrices(scan, jobs, n_qubits):
    term_mats = []

    for job in jobs:
        P = scan.pauli_matrix_from_key(job.key, n_qubits)
        H_j = complex(job.coeff) * P
        H_j = 0.5 * (H_j + H_j.conj().T)
        term_mats.append(H_j)

    return term_mats


def build_pair_commutators(term_mats, tol):
    n = len(term_mats)
    comm = {}

    dim = term_mats[0].shape[0]
    zero = np.zeros((dim, dim), dtype=complex)

    nonzero_pairs = 0

    for i in range(n):
        for j in range(i + 1, n):
            C = term_mats[i] @ term_mats[j] - term_mats[j] @ term_mats[i]
            norm_C = svdvals(C)[0]

            if norm_C > tol:
                comm[(i, j)] = C
                comm[(j, i)] = -C
                nonzero_pairs += 1

    return comm, zero, nonzero_pairs


def surrogate_for_order(order, comm, zero):
    Csum = zero.copy()

    for a, i in enumerate(order):
        for j in order[a + 1:]:
            Csum += comm.get((i, j), zero)

    return float(svdvals(Csum)[0])


def evaluate_one(mapping, dt, args, scan, outdir):
    print()
    print("=" * 80)
    print(f"Mapping: {mapping.upper()}, dt={dt}")
    print("=" * 80)

    dt_tag = safe_float_tag(dt)

    input_csv = Path(args.csv_dir) / f"h2_pauli_trotter_errors_{mapping}_dt{dt_tag}.csv"
    if not input_csv.exists():
        raise FileNotFoundError(f"Missing input CSV: {input_csv}")

    rows = read_rows(input_csv)

    print(f"Loaded exhaustive error CSV: {input_csv}")
    print(f"Rows / trace classes: {len(rows)}")

    jobs, n_qubits = build_pauli_jobs(
        scan=scan,
        ham2_file=args.ham2_file,
        mapping=mapping,
        tol=args.tol,
    )

    print(f"Rebuilt Pauli jobs: {len(jobs)}")
    print(f"n_qubits: {n_qubits}")

    term_mats = build_term_matrices(scan, jobs, n_qubits)
    comm, zero, nonzero_pairs = build_pair_commutators(term_mats, args.tol)

    print(f"Nonzero commutator pairs: {nonzero_pairs}")

    print()
    print("Computing commutator surrogate for every trace class.")

    start = time.time()

    for k, row in enumerate(rows):
        if k % 500 == 0:
            print(f"  processed {k}/{len(rows)} rows")

        s = surrogate_for_order(row["order_tuple"], comm, zero)

        row["commutator_surrogate"] = s
        row["predicted_leading_error"] = 0.5 * (dt ** 2) * s

    elapsed = time.time() - start
    print(f"Finished surrogate computation in {elapsed:.3f} seconds.")

    rows_by_actual = sorted(rows, key=lambda r: r["operator_error"])
    rows_by_surrogate = sorted(rows, key=lambda r: r["commutator_surrogate"])

    for rank, row in enumerate(rows_by_surrogate, start=1):
        row["surrogate_rank"] = rank

    actual = np.array([r["operator_error"] for r in rows], dtype=float)
    surrogate = np.array([r["commutator_surrogate"] for r in rows], dtype=float)
    predicted = np.array([r["predicted_leading_error"] for r in rows], dtype=float)

    spearman = spearmanr(surrogate, actual)
    kendall = kendalltau(surrogate, actual)

    best_actual = rows_by_actual[0]
    best_surrogate = rows_by_surrogate[0]

    topk_summary = {}
    for k in [1, 5, 10, 25, 50, 100, 250, 500]:
        candidates = rows_by_surrogate[:k]
        best_in_topk = min(candidates, key=lambda r: r["operator_error"])

        topk_summary[f"top_{k}"] = {
            "best_actual_error_in_topk": best_in_topk["operator_error"],
            "ratio_to_global_best": best_in_topk["operator_error"] / best_actual["operator_error"],
            "class_id": best_in_topk["class_id"],
            "actual_rank": best_in_topk["operator_error_rank"],
            "surrogate_rank": best_in_topk["surrogate_rank"],
            "order_indices": best_in_topk["order_indices"],
        }

    enriched_csv = outdir / f"h2_commutator_surrogate_{mapping}_dt{dt_tag}.csv"
    summary_json = outdir / f"h2_commutator_surrogate_{mapping}_dt{dt_tag}_summary.json"

    write_rows(enriched_csv, rows_by_surrogate)

    summary = {
        "mapping": mapping,
        "dt": dt,
        "input_csv": str(input_csv),
        "enriched_csv": str(enriched_csv),
        "n_trace_classes": len(rows),
        "n_pauli_jobs": len(jobs),
        "n_qubits": n_qubits,
        "nonzero_commutator_pairs": nonzero_pairs,
        "spearman_correlation": {
            "rho": float(spearman.statistic),
            "pvalue": float(spearman.pvalue),
        },
        "kendall_correlation": {
            "tau": float(kendall.statistic),
            "pvalue": float(kendall.pvalue),
        },
        "operator_error": {
            "min": float(np.min(actual)),
            "median": float(np.median(actual)),
            "mean": float(np.mean(actual)),
            "max": float(np.max(actual)),
        },
        "commutator_surrogate": {
            "min": float(np.min(surrogate)),
            "median": float(np.median(surrogate)),
            "mean": float(np.mean(surrogate)),
            "max": float(np.max(surrogate)),
        },
        "predicted_leading_error": {
            "min": float(np.min(predicted)),
            "median": float(np.median(predicted)),
            "mean": float(np.mean(predicted)),
            "max": float(np.max(predicted)),
        },
        "global_best_actual": {
            "class_id": best_actual["class_id"],
            "operator_error": best_actual["operator_error"],
            "operator_error_rank": best_actual["operator_error_rank"],
            "surrogate_rank": best_actual["surrogate_rank"],
            "commutator_surrogate": best_actual["commutator_surrogate"],
            "predicted_leading_error": best_actual["predicted_leading_error"],
            "order_indices": best_actual["order_indices"],
            "order_labels": best_actual["order_labels"],
        },
        "surrogate_best": {
            "class_id": best_surrogate["class_id"],
            "operator_error": best_surrogate["operator_error"],
            "operator_error_rank": best_surrogate["operator_error_rank"],
            "surrogate_rank": best_surrogate["surrogate_rank"],
            "ratio_to_global_best": best_surrogate["operator_error"] / best_actual["operator_error"],
            "commutator_surrogate": best_surrogate["commutator_surrogate"],
            "predicted_leading_error": best_surrogate["predicted_leading_error"],
            "order_indices": best_surrogate["order_indices"],
            "order_labels": best_surrogate["order_labels"],
        },
        "topk_summary": topk_summary,
    }

    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print(f"Saved enriched CSV: {enriched_csv}")
    print(f"Saved summary JSON: {summary_json}")

    print()
    print("Summary:")
    print(f"  Spearman rho: {summary['spearman_correlation']['rho']:.6f}")
    print(f"  Kendall tau:  {summary['kendall_correlation']['tau']:.6f}")
    print(f"  Global best actual error: {best_actual['operator_error']:.12e}")
    print(f"  Surrogate-best actual error: {best_surrogate['operator_error']:.12e}")
    print(f"  Surrogate-best / global-best ratio: {summary['surrogate_best']['ratio_to_global_best']:.6f}")
    print(f"  Surrogate-best actual rank: {best_surrogate['operator_error_rank']}")

    for key, val in topk_summary.items():
        print(
            f"  {key}: best_actual={val['best_actual_error_in_topk']:.12e}, "
            f"ratio={val['ratio_to_global_best']:.6f}, "
            f"actual_rank={val['actual_rank']}, "
            f"surrogate_rank={val['surrogate_rank']}"
        )

    return summary


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "ham2_file",
        help="Path to H2 ham2 pickle.",
    )
    parser.add_argument(
        "--scan-script",
        default="experiments/h2/scan_h2_pauli_trotter_errors.py",
        help="Path to the previous exhaustive scan script. This script imports helper functions from it.",
    )
    parser.add_argument(
        "--csv-dir",
        default="experiments/h2/trotter_error_scan_dt_sweep",
        help="Directory containing h2_pauli_trotter_errors_{mapping}_dt*.csv files.",
    )
    parser.add_argument(
        "--mapping",
        choices=["jw", "bk", "both"],
        default="both",
    )
    parser.add_argument(
        "--dt",
        type=float,
        required=True,
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-10,
    )
    parser.add_argument(
        "--outdir",
        default="experiments/h2/commutator_surrogate",
    )
    parser.add_argument(
        "--log-file",
        default=None,
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    log_path, log_fh, orig_stdout, orig_stderr = setup_logging(outdir, args.log_file)

    try:
        print(f"Log file: {log_path}")
        print(f"Started: {datetime.now().isoformat(timespec='seconds')}")
        print(f"Arguments: {vars(args)}")

        scan = import_scan_module(args.scan_script)

        if args.mapping == "both":
            mappings = ["jw", "bk"]
        else:
            mappings = [args.mapping]

        summaries = []

        for mapping in mappings:
            summaries.append(
                evaluate_one(
                    mapping=mapping,
                    dt=args.dt,
                    args=args,
                    scan=scan,
                    outdir=outdir,
                )
            )

        combined_path = outdir / f"h2_commutator_surrogate_dt{safe_float_tag(args.dt)}_combined_summary.json"

        with open(combined_path, "w") as f:
            json.dump(summaries, f, indent=2)

        print()
        print(f"Saved combined summary: {combined_path}")
        print(f"Finished: {datetime.now().isoformat(timespec='seconds')}")

    finally:
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        log_fh.flush()
        log_fh.close()


if __name__ == "__main__":
    main()
