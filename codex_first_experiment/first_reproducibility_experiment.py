#!/usr/bin/env python3
"""Read-only QHAT baseline audit. No arguments: preflight only; --execute: run.

Uses existing QHAT physics routines without modifying them. Each repetition
runs in a fresh process against a frozen tensor. Results never overwrite a
previous campaign. This tests repeatability, not historical provenance or
novel ordering performance. See FIRST_EXPERIMENT.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile


CASES = {
    "f2": "hamiltonian_generator/hgbs5_diatomic_basis_transfer/night/F-F/1.28/hgbs-5/F-F_1.28_hgbs-5_as-014-002.tensors.npz",
    "nh3": "hamiltonian_generator/polyatomic_library/NH3/s-1.00/hgbs-5/NH3_s-1.00_hgbs-5_as-006-010.tensors.npz",
    "li2": "hamiltonian_generator/hgbs5_diatomic_basis_transfer/lower/Li-Li/2.66/hgbs-5/Li-Li_2.66_hgbs-5_as-006-006.tensors.npz",
    "smoke": "hamiltonian_generator/hgbs5_diatomic_basis_transfer/lower/B-B/1.70/hgbs-5/B-B_1.70_hgbs-5_as-002-002.tensors.npz",
}
ORDERINGS = (
    "fermionic_signed_coefficient_lexicographic",
    "signed_coefficient_lexicographic",
    "jw_magnitude_descending_lexicographic",
)
HISTORICAL_SCHEDULES = dict(zip(
    ("fermionic_signed_reference", "jw_signed_baseline", "jw_magnitude_baseline"),
    ORDERINGS,
))
HISTORY = (
    "analysis/cancellation_hypothesis_validation/full_ablation_results.csv",
    "analysis/fermionic_body_rank_ablation_20case.csv",
)
METRICS = ("one_minus_overlap", "bch2_hf_state_norm", "number_of_pauli_terms",
           "number_of_fermionic_terms")


def validate_tensor_archive(path: Path) -> None:
    with zipfile.ZipFile(path) as archive:
        missing = {"constant.npy", "one_body.npy", "two_body.npy"} - set(archive.namelist())
        if missing:
            raise ValueError(f"Tensor is incompatible with the current QHAT loader: {path}; missing {sorted(missing)}")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def save(path: Path, obj: object) -> None:
    # Exclusive creation: a prior result is never overwritten.
    with path.open("x", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write("\n")


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def source_fingerprints(repo: Path) -> dict[str, str]:
    paths = git(repo, "ls-files", "--", "*.py", "pyproject.toml").splitlines()
    return {name: digest(repo / name) for name in paths}


def assert_sources(repo: Path, manifest: dict) -> None:
    if source_fingerprints(repo) != manifest["source_sha256"]:
        raise RuntimeError("QHAT sources changed during the campaign; start a new campaign.")
    if digest(Path(__file__)) != manifest["driver_sha256"]:
        raise RuntimeError("Audit script changed during the campaign.")


def copy_verified(source: Path, target: Path) -> str:
    before = digest(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, target.open("xb") as dst:
        shutil.copyfileobj(src, dst)
    if before != digest(source) or before != digest(target):
        raise RuntimeError(f"Input changed while copying: {source}")
    return before


def compare_metric(name: str, fresh: float, reference: float, *, repeat: bool) -> str:
    if not all(math.isfinite(v) for v in (fresh, reference)):
        return "nonfinite"
    if name.startswith("number_of_"):
        return "match" if fresh == reference else "different"
    # This is an audit threshold, not a claim about scientific significance.
    atol = 1e-13 if name == "one_minus_overlap" else 1e-10
    rtol = 1e-5 if repeat else 1e-3
    if name == "one_minus_overlap" and min(fresh, reference) <= 1e-12:
        return "below_floor_match" if math.isclose(fresh, reference, abs_tol=atol,
                                                   rel_tol=rtol) else "below_floor_different"
    return "match" if math.isclose(fresh, reference, abs_tol=atol, rel_tol=rtol) else "different"


def history_rows(campaign: Path, manifest: dict) -> list[dict]:
    rows = []
    selected = {v["case_id"] for v in manifest["cases"].values()}
    for name in manifest["history_sha256"]:
        with (campaign / "history" / name).open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("status") != "success" or row.get("case_id") not in selected:
                    continue
                ordering = HISTORICAL_SCHEDULES.get(row.get("schedule"))
                if not ordering:
                    continue
                settings = manifest["settings"]
                if any(float(row[col]) != settings[key] for col, key in (
                    ("trotter_steps", "steps"), ("evolution_time", "evolution_time"),
                    ("coefficient_tolerance", "tolerance"))):
                    continue
                rows.append({**row, "ordering": ordering, "source": name})
    return rows


def summarize(campaign: Path, manifest: dict, runs: list[dict], errors: list[str]) -> bool:
    comparisons = []
    repeat_failures = []
    for case in manifest["cases"]:
        group = [r for r in runs if r["case"] == case]
        if len(group) != manifest["repeats"]:
            repeat_failures.append(f"{case}: incomplete independent repetitions")
        if not group:
            continue
        first = group[0]
        for other in group[1:]:
            if other["fingerprints"] != first["fingerprints"]:
                repeat_failures.append(f"{case}: tensor/Hamiltonian/parent/order hashes changed")
            if other["versions"] != first["versions"]:
                repeat_failures.append(f"{case}: dependency versions changed")
            for a, b in zip(first["rows"], other["rows"], strict=True):
                for metric in METRICS:
                    result = compare_metric(metric, float(b[metric]), float(a[metric]), repeat=True)
                    comparisons.append(dict(kind="repeat", case=case, ordering=a["ordering"],
                                            metric=metric, first=a[metric], fresh=b[metric], result=result))
                    if result not in ("match", "below_floor_match"):
                        repeat_failures.append(f"{case}/{a['ordering']}/{metric}: {result}")
    for old in history_rows(campaign, manifest):
        for run in runs:
            if run["repetition"] != 1 or run["rows"][0]["case_id"] != old["case_id"]:
                continue
            fresh = next(r for r in run["rows"] if r["ordering"] == old["ordering"])
            for metric in METRICS:
                comparisons.append(dict(kind="historical", case=run["case"], source=old["source"],
                    ordering=old["ordering"], metric=metric, historical=float(old[metric]),
                    fresh=float(fresh[metric]), result=compare_metric(
                        metric, float(fresh[metric]), float(old[metric]), repeat=False)))
    save(campaign / "comparisons.json", comparisons)
    stable = not errors and not repeat_failures
    changed = [r for r in comparisons if r["kind"] == "historical" and r["result"]
               not in ("match", "below_floor_match")]
    lines = ["# First reproducibility experiment", "",
             f"Independent-repeat checks: {'PASS' if stable else 'FAIL / INCOMPLETE'}",
             f"Completed worker jobs: {len(runs)}/{len(manifest['cases']) * manifest['repeats']}",
             f"Historical metric comparisons flagged: {len(changed)}", "",
             "Primary physical error: 1 - abs(overlap). Raw vector error is NOT used:",
             "the existing QHAT exact/Trotter routines use different scalar energy offsets.", "",
             "Repeat PASS establishes only reproducibility with the CURRENT frozen input/code.",
             "Historical agreement does not establish identical historical tensors or code.",
             "Historical disagreement does not identify its cause. Never replace old results.", "",
             "## Historical differences", ""]
    for row in changed:
        lines.append(f"- {row['case']} / {row['source']} / {row['ordering']} / {row['metric']}: "
                     f"historical={row['historical']:.8g}, current={row['fresh']:.8g} ({row['result']})")
    if not changed:
        lines.append("No flagged historical comparisons (see comparisons.json for coverage).")
    lines += ["", "## Execution / repeat failures", ""] + [f"- {e}" for e in errors + repeat_failures]
    lines += ["", "## Next decision", "",
              "If repeats fail: investigate environment, numerical tolerances, and code first.",
              "If repeats pass but history differs: trace tensor/orbital/code revisions before pooling data.",
              "Do not start a large new sweep until input lineage and scalar-phase conventions are resolved."]
    with (campaign / "SUMMARY.md").open("x", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return stable


def worker(campaign: Path, case: str, repetition: int) -> None:
    # Parent sets thread counts/cache locations before this process imports NumPy.
    manifest = json.loads((campaign / "manifest.json").read_text())
    repo = Path(manifest["repo"])
    assert_sources(repo, manifest)
    sys.path[:0] = [str(repo.parent), str(repo)]
    import numpy as np
    from openfermion import get_fermion_operator, jordan_wigner
    from analysis import benchmark_b2_signed_coefficient_baseline as baseline
    for fn in (baseline.benchmark_case, baseline.exact_reference_state,
               baseline.build_hermitian_fermion_terms):
        module_file = Path(sys.modules[fn.__module__].__file__).resolve()
        if not module_file.is_relative_to(repo):
            raise RuntimeError(f"Imported QHAT code outside selected repo: {module_file}")
    item = manifest["cases"][case]
    tensor = campaign / "inputs" / item["relative_path"]
    if digest(tensor) != item["sha256"]:
        raise RuntimeError("Frozen tensor checksum mismatch")
    args = argparse.Namespace(**manifest["settings"])
    interaction, n = baseline.load_interaction_operator(tensor)
    if n > manifest["max_qubits"]:
        raise RuntimeError(f"{n} qubits exceeds explicit cap {manifest['max_qubits']}")
    fermion = baseline.clean_fermion_operator(get_fermion_operator(interaction), args.tolerance)
    jw = jordan_wigner(fermion)
    jw.compress(abs_tol=args.tolerance)
    coefficients = {k: v for k, v in jw.terms.items() if k and abs(v) > args.tolerance}
    parents = baseline.build_hermitian_fermion_terms(fermion, args.tolerance)
    orders = baseline.build_deterministic_orderings(
        fermion, list(coefficients), coefficients, n, args.tolerance)
    def canonical_terms(terms):
        return [[key, float(complex(value).real).hex(), float(complex(value).imag).hex()]
                for key, value in sorted(terms.items())]
    parent_terms = sorted([canonical_terms(p.operator.terms) for p in parents], key=repr)
    fingerprints = {
        "tensor_sha256": digest(tensor),
        "identity_free_jw_sha256": json_digest(canonical_terms(coefficients)),
        "hermitian_parents_sha256": json_digest(parent_terms),
        "physical_order_sha256": {name: json_digest([
            [key, complex(coefficients[key]).real.hex(), complex(coefficients[key]).imag.hex()]
            for key in orders[name]]) for name in ORDERINGS},
    }
    save(campaign / f"{case}.repeat{repetition}.operators.json", {
        "pauli_coefficients": canonical_terms(coefficients), "parents": parent_terms,
        "orders": {name: orders[name] for name in ORDERINGS}})
    baseline.warm_up_numba()
    started = time.monotonic()
    rows = baseline.benchmark_case(tensor, args, ordering_names=ORDERINGS)
    for row in rows:
        row["number_of_fermionic_terms"] = len(parents)
    if [r["ordering"] for r in rows] != list(ORDERINGS):
        raise RuntimeError("Missing or reordered benchmark outputs")
    if any(r["status"] != "success" or not all(math.isfinite(float(r[m])) for m in METRICS)
           or not 0 <= float(r["one_minus_overlap"]) <= 1 for r in rows):
        raise RuntimeError("Failed or invalid numerical result")
    assert_sources(repo, manifest)
    if digest(tensor) != item["sha256"]:
        raise RuntimeError("Frozen tensor changed during calculation")
    versions = {p: importlib.metadata.version(p) for p in
                ("numpy", "scipy", "openfermion", "numba", "pandas")}
    versions["python"] = sys.version
    save(campaign / f"{case}.repeat{repetition}.json", {
        "case": case, "repetition": repetition, "fingerprints": fingerprints,
        "versions": versions, "number_of_fermionic_terms": len(parents),
        "exact_minus_trotter_scalar": float((jw.terms.get((), 0) - fermion.terms.get((), 0)).real),
        "benchmark_wall_seconds": time.monotonic() - started, "rows": rows})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("/Users/albertlee0125/Repos/qhat"))
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parent / "repro_runs")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=["f2", "nh3", "li2"])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout-minutes", type=float, default=30)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--time", type=float, default=1.0, dest="evolution_time")
    parser.add_argument("--tolerance", type=float, default=1e-12)
    parser.add_argument("--max-qubits", type=int, default=16)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--worker", nargs=3, metavar=("CAMPAIGN", "CASE", "REPETITION"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(Path(args.worker[0]), args.worker[1], int(args.worker[2]))
        return 0
    if (args.repeats < 2 or args.threads < 1 or args.steps < 1 or args.max_qubits < 1
        or len(args.cases) != len(set(args.cases))
        or any(not math.isfinite(v) or v <= 0 for v in
               (args.timeout_minutes, args.evolution_time, args.tolerance))):
        parser.error("Use >=2 repeats, positive finite settings, and unique cases.")
    repo = args.repo.resolve()
    commit = git(repo, "rev-parse", "HEAD")
    for case in args.cases:
        tensor = repo / CASES[case]
        if not tensor.is_file():
            raise FileNotFoundError(f"Missing input; no regeneration will be attempted: {tensor}")
        validate_tensor_archive(tensor)
        counts = tensor.name.split("_as-")[1].split(".")[0].split("-")
        n = sum(map(int, counts))
        if n > args.max_qubits:
            raise ValueError(f"{case}: {n} qubits exceeds --max-qubits")
        print(f"{case}: {n}q, {tensor.stat().st_size} bytes, SHA256={digest(tensor)}")
    print(f"QHAT commit: {commit}\nPlan: {len(args.cases)} cases x {len(ORDERINGS)} orders x "
          f"{args.repeats} independent repeats = {len(args.cases)*len(ORDERINGS)*args.repeats} rows.")
    print(f"First order, HF state, T={args.evolution_time}, r={args.steps}; "
          f"serial jobs, {args.threads} thread(s), {args.timeout_minutes} min timeout per job.")
    if not args.execute:
        print("Preflight only. Nothing written. Add --execute using your QHAT Python environment.")
        return 0
    args.output_root.mkdir(parents=True, exist_ok=True)
    campaign = Path(tempfile.mkdtemp(prefix="first-repro-", dir=args.output_root.resolve()))
    print(f"Campaign: {campaign}", flush=True)
    manifest = {
        "repo": str(repo), "commit": commit, "branch": git(repo, "branch", "--show-current"),
        "tracked_worktree_status": git(repo, "status", "--porcelain", "--untracked-files=no"),
        "source_sha256": source_fingerprints(repo), "driver_sha256": digest(Path(__file__)),
        "python_executable": sys.executable, "platform": platform.platform(),
        "repeats": args.repeats, "threads": args.threads, "max_qubits": args.max_qubits,
        "settings": dict(steps=args.steps, evolution_time=args.evolution_time,
                         tolerance=args.tolerance, parallel_threshold=2**16, no_spin_sector=False),
        "cases": {}, "history_sha256": {},
    }
    for case in args.cases:
        rel = CASES[case]
        # Preserve the molecule/geometry/basis directory layout used by QHAT metadata.
        target_rel = rel
        sha = copy_verified(repo / rel, campaign / "inputs" / target_rel)
        manifest["cases"][case] = dict(relative_path=target_rel, source=str(repo / rel),
                                     case_id=Path(rel).name.removesuffix(".tensors.npz"), sha256=sha)
    for name in HISTORY:
        manifest["history_sha256"][name] = copy_verified(repo / name, campaign / "history" / name)
    save(campaign / "manifest.json", manifest)
    environment = os.environ.copy()
    environment.update({key: str(args.threads) for key in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                         "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS")})
    environment.update(PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0",
                       MPLCONFIGDIR=str(campaign / "cache" / "mpl"),
                       NUMBA_CACHE_DIR=str(campaign / "cache" / "numba"))
    runs, errors = [], []
    for case in args.cases:
        for repetition in range(1, args.repeats + 1):
            print(f"Running {case}, repeat {repetition}/{args.repeats} ...", flush=True)
            try:
                assert_sources(repo, manifest)
                with (campaign / f"{case}.repeat{repetition}.log").open("x") as log:
                    subprocess.run([sys.executable, "-u", str(Path(__file__).resolve()), "--worker",
                                    str(campaign), case, str(repetition)], cwd=repo, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, check=True,
                                   timeout=60 * args.timeout_minutes)
                runs.append(json.loads((campaign / f"{case}.repeat{repetition}.json").read_text()))
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                errors.append(f"{case} repeat {repetition}: {exc}; inspect the corresponding log")
                break
        if errors:
            break  # A failed prerequisite should not launch more expensive jobs.
    stable = summarize(campaign, manifest, runs, errors)
    print(f"{'PASS' if stable else 'FAIL / INCOMPLETE'}: {campaign / 'SUMMARY.md'}")
    return 0 if stable else 1


if __name__ == "__main__":
    raise SystemExit(main())
