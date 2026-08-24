#!/usr/bin/env python3
"""Memory-efficient lowest-energy solver for QHAT/OpenFermion InteractionOperator pickles.

This program avoids constructing the 2**n_qubits by 2**n_qubits sparse matrix.
Instead, it converts the spin-orbital tensors in an OpenFermion
InteractionOperator back to spatial-orbital integrals and calls PySCF's direct
FCI Davidson solver in a fixed (N_alpha, N_beta) electron sector.

The method is CPU based. It intentionally leaves the GPU mostly idle because
building the full OpenFermion sparse matrix is the operation that exhausted
host RAM. PySCF contracts H with CI vectors directly and never stores the full
Hamiltonian matrix.

Typical use:
  python qhat_low_energies_fci.py hamiltonian.pickle -o q24_fci \
      --num-electrons 10 --spin 0 -k 10 --threads 32

Here --spin means N_alpha - N_beta = 2*S_z. For an even-electron singlet use 0.
For a doublet use 1, for a triplet M_s=1 sector use 2, and so on.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np


def human_bytes(n: int | float) -> str:
    x = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if x < 1024.0 or unit == "PiB":
            return f"{x:.3f} {unit}"
        x /= 1024.0
    return f"{x:.3f} PiB"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compute the lowest energies of a QHAT/OpenFermion "
            "InteractionOperator without constructing the full sparse matrix."
        ),
    )
    p.add_argument("input", type=Path, help="QHAT/OpenFermion pickle")
    p.add_argument("-o", "--output-dir", type=Path, default=Path("fci_output"))
    p.add_argument("-k", "--num-roots", type=int, default=5,
                   help="Number of lowest energies/eigenstates")
    p.add_argument("--num-electrons", type=int, required=True,
                   help="Total number of active electrons")
    p.add_argument("--spin", type=int, default=0,
                   help="N_alpha - N_beta (twice S_z), not spin multiplicity")
    p.add_argument("--threads", type=int, default=0,
                   help="OpenMP/BLAS threads; 0 leaves environment unchanged")
    p.add_argument("--tol", type=float, default=1e-9,
                   help="Davidson convergence tolerance")
    p.add_argument("--max-cycle", type=int, default=100)
    p.add_argument("--max-space", type=int, default=20,
                   help="Maximum Davidson subspace size before restart")
    p.add_argument("--lindep", type=float, default=1e-14)
    p.add_argument("--level-shift", type=float, default=1e-3)
    p.add_argument("--pspace-size", type=int, default=400,
                   help="Small determinant-space preconditioner size")
    p.add_argument("--max-memory-mb", type=int, default=0,
                   help="PySCF memory limit in MB; 0 uses 80%% of physical RAM")
    p.add_argument("--verbose", type=int, default=5,
                   help="PySCF verbosity; 4-5 shows Davidson progress")
    p.add_argument("--save-vectors", action="store_true",
                   help="Save CI eigenvectors; energies are always saved")
    p.add_argument("--inspect-only", action="store_true")
    p.add_argument("--imag-tol", type=float, default=1e-10,
                   help="Maximum allowed imaginary integral magnitude")
    p.add_argument("--spin-block-tol", type=float, default=1e-8,
                   help="Tolerance used to validate even/odd spin blocks")
    p.add_argument("--skip-spin-block-check", action="store_true",
                   help="Skip validation that the pickle uses OpenFermion even/odd spin ordering")
    p.add_argument("--validate-small", action="store_true",
                   help="For small cases, compare PySCF roots with an OpenFermion sector sparse solve")
    return p.parse_args()


def physical_memory_bytes() -> int:
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return int(pages * page_size)
    except (ValueError, OSError, AttributeError):
        return 0


def configure_threads(nthreads: int) -> None:
    if nthreads <= 0:
        return
    text = str(nthreads)
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = text


def load_interaction_operator(path: Path) -> Any:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("rb") as f:
        obj = pickle.load(f)
    required = ("constant", "one_body_tensor", "two_body_tensor")
    missing = [name for name in required if not hasattr(obj, name)]
    if missing:
        raise TypeError(
            f"Input type {type(obj)!r} is not an InteractionOperator-like object; "
            f"missing attributes: {missing}"
        )
    return obj


def split_electrons(nelec: int, spin: int, norb: int) -> tuple[int, int]:
    if nelec < 0 or nelec > 2 * norb:
        raise ValueError(f"--num-electrons must be between 0 and {2*norb}")
    if abs(spin) > nelec:
        raise ValueError("abs(--spin) cannot exceed --num-electrons")
    if (nelec + spin) % 2 != 0:
        raise ValueError(
            "--num-electrons and --spin must have the same parity so that "
            "N_alpha=(N+spin)/2 is an integer"
        )
    nalpha = (nelec + spin) // 2
    nbeta = nelec - nalpha
    if nalpha > norb or nbeta > norb or nalpha < 0 or nbeta < 0:
        raise ValueError(
            f"Invalid sector for {norb} spatial orbitals: "
            f"N_alpha={nalpha}, N_beta={nbeta}"
        )
    return nalpha, nbeta


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if a.size else 0.0


def validate_spin_blocks(one: np.ndarray, two: np.ndarray, tol: float) -> dict[str, float]:
    """Validate OpenFermion's default even/odd spin-orbital ordering.

    MolecularData.get_molecular_hamiltonian uses indices 2*p and 2*p+1 for
    alpha/beta versions of spatial orbital p. The one-body alpha and beta
    blocks should agree and spin-flip blocks should vanish. The same-spin
    samples of the two-body tensor should also agree.
    """
    even = slice(0, None, 2)
    odd = slice(1, None, 2)
    metrics = {
        "one_alpha_beta_difference": max_abs(one[even, even] - one[odd, odd]),
        "one_spin_flip_even_odd": max_abs(one[even, odd]),
        "one_spin_flip_odd_even": max_abs(one[odd, even]),
        "two_eeee_oooo_difference": max_abs(
            two[even, even, even, even] - two[odd, odd, odd, odd]
        ),
    }
    worst = max(metrics.values())
    if worst > tol:
        details = ", ".join(f"{k}={v:.3e}" for k, v in metrics.items())
        raise ValueError(
            "The tensors do not look like a standard spin-independent molecular "
            "InteractionOperator in even/odd spin ordering. " + details + ". "
            "Do not bypass this check unless you know the tensor convention."
        )
    return metrics


def interaction_to_pyscf(obj: Any, imag_tol: float, spin_tol: float,
                         skip_check: bool) -> tuple[float, np.ndarray, np.ndarray, dict[str, float]]:
    """Convert OpenFermion InteractionOperator tensors to PySCF spatial integrals.

    OpenFermion MolecularData creates the InteractionOperator as
      one_body_tensor = spinorb_from_spatial(h1)
      two_body_tensor = 0.5 * spinorb_from_spatial(h2_openfermion)

    Therefore the spatial OpenFermion tensor is recovered from an all-alpha
    block by multiplying by two. OpenFermion's spatial two-electron tensor is
    converted to chemist ordering for PySCF with transpose(0,3,1,2).
    """
    one = np.asarray(obj.one_body_tensor)
    two = np.asarray(obj.two_body_tensor)
    if one.ndim != 2 or one.shape[0] != one.shape[1]:
        raise ValueError(f"Unexpected one_body_tensor shape {one.shape}")
    nspin = one.shape[0]
    if nspin % 2:
        raise ValueError(f"Expected an even number of spin orbitals; got {nspin}")
    if two.shape != (nspin, nspin, nspin, nspin):
        raise ValueError(f"Unexpected two_body_tensor shape {two.shape}")

    imag = max(max_abs(np.asarray(obj.constant).imag), max_abs(one.imag), max_abs(two.imag))
    if imag > imag_tol:
        raise ValueError(
            f"Complex integrals detected (maximum imaginary magnitude {imag:.3e}). "
            "This script currently targets the real molecular Hamiltonians handled "
            "efficiently by PySCF direct_spin1."
        )
    one = np.asarray(one.real, dtype=np.float64, order="C")
    two = np.asarray(two.real, dtype=np.float64, order="C")

    metrics: dict[str, float] = {}
    if not skip_check:
        metrics = validate_spin_blocks(one, two, spin_tol)

    # Default OpenFermion ordering is alpha,beta,alpha,beta,...
    h1 = np.asarray(one[0::2, 0::2], dtype=np.float64, order="C")
    h2_of = np.asarray(2.0 * two[0::2, 0::2, 0::2, 0::2], dtype=np.float64)
    # OpenFermion spatial ordering -> chemist ordering (pq|rs).
    eri = np.asarray(h2_of.transpose(0, 3, 1, 2), dtype=np.float64, order="C")
    ecore = float(np.real(obj.constant))
    return ecore, h1, eri, metrics


def sector_size(norb: int, nalpha: int, nbeta: int) -> tuple[int, int, int]:
    na = math.comb(norb, nalpha)
    nb = math.comb(norb, nbeta)
    return na, nb, na * nb


def normalize_roots(energies: Any, vectors: Any, k: int) -> tuple[np.ndarray, list[np.ndarray]]:
    e = np.atleast_1d(np.asarray(energies, dtype=np.float64))
    if k == 1:
        vecs = [np.asarray(vectors)]
    else:
        vecs = [np.asarray(v) for v in vectors]
    order = np.argsort(e)
    return e[order], [vecs[int(i)] for i in order]


def save_results(args: argparse.Namespace, energies: np.ndarray,
                 vectors: list[np.ndarray], metadata: dict[str, Any]) -> None:
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "eigenvalues.npy", energies)
    np.savetxt(
        out / "lowest_energies.txt",
        np.column_stack((np.arange(energies.size), energies)),
        header="root_index energy_hartree",
        fmt=["%d", "%.16e"],
    )
    # This is the returned low-energy diagonal only, not a full-spectrum matrix.
    np.save(out / "diagonal_values.npy", energies)
    if args.save_vectors:
        for i, vec in enumerate(vectors):
            np.save(out / f"ci_vector_{i:04d}.npy", vec)
    with (out / "metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def validate_against_openfermion(obj: Any, energies: np.ndarray, nspin: int,
                                 nelec: int, spin: int, k: int) -> dict[str, Any]:
    """Optional small-problem cross-check against OpenFermion sector matrix."""
    if nspin > 18:
        raise ValueError("--validate-small is restricted to at most 18 spin orbitals")
    try:
        from openfermion.transforms import get_fermion_operator
        from openfermion.linalg import get_number_preserving_sparse_operator
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:
        raise RuntimeError("Validation requires OpenFermion and SciPy") from exc

    if spin != 0:
        raise ValueError(
            "The automatic validation currently supports --spin 0 only because "
            "OpenFermion's helper is being used with spin_preserving=True."
        )
    fermion = get_fermion_operator(obj)
    reference = np.zeros(nspin, dtype=bool)
    nalpha = nelec // 2
    # OpenFermion even indices are alpha and odd are beta.
    reference[2 * np.arange(nalpha)] = True
    reference[2 * np.arange(nalpha) + 1] = True
    h = get_number_preserving_sparse_operator(
        fermion,
        num_qubits=nspin,
        num_electrons=nelec,
        spin_preserving=True,
        reference_determinant=reference,
    ).tocsr()
    kval = min(k, h.shape[0] - 1)
    ref = np.sort(eigsh(h, k=kval, which="SA", return_eigenvectors=False))
    err = float(np.max(np.abs(ref - energies[:kval])))
    return {
        "validation_sector_shape": list(h.shape),
        "validation_max_abs_energy_difference": err,
        "validation_reference_energies": ref.tolist(),
    }


def main() -> int:
    args = parse_args()
    configure_threads(args.threads)

    try:
        from pyscf import fci, lib
    except ImportError as exc:
        raise RuntimeError(
            "PySCF is required. Install it with: python -m pip install pyscf"
        ) from exc

    if args.threads > 0:
        lib.num_threads(args.threads)

    t0 = time.perf_counter()
    obj = load_interaction_operator(args.input)
    load_seconds = time.perf_counter() - t0
    ecore, h1, eri, spin_metrics = interaction_to_pyscf(
        obj, args.imag_tol, args.spin_block_tol, args.skip_spin_block_check
    )
    nspin = int(np.asarray(obj.one_body_tensor).shape[0])
    del obj
    norb = h1.shape[0]
    nalpha, nbeta = split_electrons(args.num_electrons, args.spin, norb)
    na, nb, dim = sector_size(norb, nalpha, nbeta)

    print(f"Spin orbitals/qubits: {nspin}", flush=True)
    print(f"Spatial orbitals: {norb}", flush=True)
    print(f"Electron sector: N_alpha={nalpha}, N_beta={nbeta}, total={nalpha+nbeta}", flush=True)
    print(f"CI dimension: C({norb},{nalpha}) * C({norb},{nbeta}) = {na:,} * {nb:,} = {dim:,}", flush=True)
    print(f"One float64 CI vector: {human_bytes(dim * 8)}", flush=True)
    print(f"Requested roots: {args.num_roots}", flush=True)
    print(f"Saved-root storage if requested: {human_bytes(dim * 8 * args.num_roots)}", flush=True)
    print(f"Integral storage: h1={human_bytes(h1.nbytes)}, eri={human_bytes(eri.nbytes)}", flush=True)

    if args.num_roots < 1 or args.num_roots > dim:
        raise ValueError(f"-k/--num-roots must be in [1, {dim}]")

    if args.max_memory_mb <= 0:
        phys = physical_memory_bytes()
        max_memory_mb = int(0.80 * phys / (1024**2)) if phys else 120000
    else:
        max_memory_mb = args.max_memory_mb
    print(f"PySCF max_memory: {max_memory_mb:,} MB", flush=True)
    print(f"PySCF/OpenMP threads: {lib.num_threads()}", flush=True)

    if args.inspect_only:
        return 0

    solver = fci.direct_spin1.FCI()
    solver.verbose = args.verbose
    solver.conv_tol = args.tol
    solver.max_cycle = args.max_cycle
    solver.max_space = args.max_space
    solver.lindep = args.lindep
    solver.level_shift = args.level_shift
    solver.pspace_size = args.pspace_size
    solver.max_memory = max_memory_mb
    solver.nroots = args.num_roots
    solver.davidson_only = True

    print("Starting direct FCI Davidson solve; no full Hamiltonian matrix is constructed.", flush=True)
    solve_start = time.perf_counter()
    energies_raw, vectors_raw = solver.kernel(
        h1,
        eri,
        norb,
        (nalpha, nbeta),
        ecore=ecore,
        nroots=args.num_roots,
    )
    solve_seconds = time.perf_counter() - solve_start
    energies, vectors = normalize_roots(energies_raw, vectors_raw, args.num_roots)

    converged_raw = solver.converged
    if isinstance(converged_raw, (list, tuple, np.ndarray)):
        converged = [bool(x) for x in converged_raw]
    else:
        converged = [bool(converged_raw)]

    metadata: dict[str, Any] = {
        "input": str(args.input),
        "method": "PySCF direct_spin1 FCI Davidson, matrix-free Hamiltonian contraction",
        "spin_orbitals": nspin,
        "spatial_orbitals": norb,
        "num_electrons": args.num_electrons,
        "spin_nalpha_minus_nbeta": args.spin,
        "nalpha": nalpha,
        "nbeta": nbeta,
        "ci_dimension": dim,
        "num_roots_requested": args.num_roots,
        "num_roots_returned": int(energies.size),
        "converged": converged,
        "conv_tol": args.tol,
        "max_cycle": args.max_cycle,
        "max_space": args.max_space,
        "threads": int(lib.num_threads()),
        "max_memory_mb": max_memory_mb,
        "constant_energy": ecore,
        "spin_block_validation": spin_metrics,
        "load_and_conversion_seconds": load_seconds,
        "solve_seconds": solve_seconds,
        "total_seconds": time.perf_counter() - t0,
        "vectors_saved": bool(args.save_vectors),
        "diagonal_meaning": "Only the returned lowest-energy diagonal entries, not the complete spectrum",
    }

    if args.validate_small:
        metadata.update(validate_against_openfermion(
            load_interaction_operator(args.input), energies, nspin,
            args.num_electrons, args.spin, args.num_roots
        ))

    save_results(args, energies, vectors, metadata)

    print("\nLowest energies (Hartree):", flush=True)
    for i, energy in enumerate(energies):
        flag = "converged" if i < len(converged) and converged[i] else "NOT CONVERGED"
        print(f"  {i:4d}  {energy: .16e}  {flag}", flush=True)
    print(f"Solve time: {solve_seconds:.3f} s", flush=True)
    print(f"Output: {args.output_dir.resolve()}", flush=True)

    if not all(converged):
        print(
            "WARNING: one or more roots did not meet the requested convergence tolerance. "
            "Increase --max-cycle and/or --max-space, or request fewer roots.",
            file=sys.stderr,
            flush=True,
        )
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr, flush=True)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        traceback.print_exc()
        raise SystemExit(1)
