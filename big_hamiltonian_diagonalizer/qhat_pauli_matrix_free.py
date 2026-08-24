#!/usr/bin/env python3
"""Matrix-free GPU low-energy solver for QHAT/OpenFermion Hamiltonians.

The program reads the same QHAT pickle used by qhat_gpu_diagonalize_v2.py,
converts an InteractionOperator/FermionOperator to an OpenFermion
QubitOperator, and applies the Pauli sum directly to vectors on an NVIDIA GPU.
It NEVER constructs the exponentially large Hamiltonian as a dense or CSR
matrix.

Two basis modes are supported:

1. Full qubit space (default): dimension 2**n_qubits.
2. Fixed-particle-number space (--num-electrons N): dimension C(n_qubits, N).
   This is strongly recommended for QHAT molecular Hamiltonians.  The Pauli
   action is projected into that sector, which is exact when the Hamiltonian
   conserves particle number.

Outputs (same core layout as qhat_gpu_diagonalize_v2.py):
  eigenvalues.npy
  eigenvectors.npy             unless --values-only
  diagonal_sparse.npz          diagonal matrix of returned eigenvalues
  residuals.npy
  metadata.json
  timing.json
  pauli_terms.npz              with --save-hamiltonian/--save-pauli-terms

Important limitations
---------------------
* This computes only the requested lowest k eigenpairs.  It does not produce
  the complete spectrum unless k is essentially the full dimension, which is
  not practical for large systems.
* Matrix-free removes the huge CSR allocation, but every Krylov vector still
  has one value per basis state.  Therefore ncv and k still determine GPU
  memory use.
* A Pauli-sum matvec can be compute intensive.  The implementation groups terms
  with the same X/Y flip mask so that basis-index mapping is performed once per
  group rather than once per Pauli term.

Tested target environment:
  Python 3.11, CuPy 14.x, CUDA 13.x, NVIDIA Blackwell.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import sys
import threading
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np


# ------------------------------- utilities ---------------------------------


def human_bytes(value: int | float) -> str:
    x = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if x < 1024.0 or unit == "PiB":
            return f"{x:.3f} {unit}"
        x /= 1024.0
    return f"{x:.3f} PiB"


@dataclass
class Timings:
    start: float = field(default_factory=time.perf_counter)
    stages: dict[str, float] = field(default_factory=dict)
    _last: float = field(init=False)

    def __post_init__(self) -> None:
        self._last = self.start

    def mark(self, name: str) -> None:
        now = time.perf_counter()
        self.stages[name] = now - self._last
        self._last = now
        print(f"[{now-self.start:10.3f} s total] {name}", flush=True)

    def as_dict(self) -> dict[str, Any]:
        return {
            "total_seconds": time.perf_counter() - self.start,
            "stage_seconds": self.stages,
        }


class SolverHeartbeat:
    def __init__(self, operator: "PauliMatrixFreeOperator", interval: float) -> None:
        self.operator = operator
        self.interval = max(1.0, interval)
        self.stop_event = threading.Event()
        self.thread: Optional[threading.Thread] = None
        self.start = 0.0

    def __enter__(self) -> "SolverHeartbeat":
        self.start = time.perf_counter()

        def worker() -> None:
            while not self.stop_event.wait(self.interval):
                elapsed = time.perf_counter() - self.start
                print(
                    f"[working {elapsed:10.1f} s] eigsh; "
                    f"matrix-free Hx calls={self.operator.matvec_calls}, "
                    f"last Hx={self.operator.last_matvec_seconds:.3f} s",
                    flush=True,
                )

        self.thread = threading.Thread(target=worker, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.stop_event.set()
        if self.thread is not None:
            self.thread.join(timeout=2.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Find low energies of a QHAT/OpenFermion Hamiltonian with a "
            "matrix-free GPU Pauli-string operator."
        ),
    )
    p.add_argument("input", type=Path, help="QHAT/OpenFermion pickle")
    p.add_argument("-o", "--output-dir", type=Path, default=Path("diag_output"))
    p.add_argument("--gpu", type=int, default=0)
    p.add_argument("--inspect-only", action="store_true")
    p.add_argument("--n-qubits", type=int, default=None)

    # Keep compatible names with the earlier program.
    p.add_argument(
        "--solver",
        choices=("auto", "sparse", "matrix-free"),
        default="auto",
        help="All choices select the matrix-free eigsh implementation.",
    )
    p.add_argument(
        "-k", "--num-eigenpairs", type=int, default=8,
        help="Number of lowest energies/eigenpairs to return.",
    )
    p.add_argument("--values-only", action="store_true")
    p.add_argument("--dtype", choices=("complex128", "complex64"), default="complex128")
    p.add_argument("--tol", type=float, default=1e-9)
    p.add_argument("--maxiter", type=int, default=None)
    p.add_argument(
        "--ncv", type=int, default=None,
        help="Krylov subspace size. Must be greater than k."
    )
    p.add_argument("--heartbeat", type=float, default=15.0)
    p.add_argument(
        "--progress-every", type=int, default=1,
        help="Synchronize/timestamp every N matrix-free Hx calls; 0 disables."
    )

    p.add_argument(
        "--num-electrons", "--num_electrons", type=int, default=None,
        help=(
            "Restrict to this fixed Hamming-weight/electron-number sector. "
            "For QHAT use the active occupied single-occupancy count from the .dat file."
        ),
    )
    p.add_argument(
        "--dat-file", type=Path, default=None,
        help="Read --num-electrons automatically from the matching QHAT .dat header."
    )
    p.add_argument(
        "--full-space", action="store_true",
        help="Explicitly use all 2**n basis states even if --dat-file is supplied."
    )

    p.add_argument("--check-hermitian", action="store_true")
    p.add_argument("--residual-checks", type=int, default=5)
    p.add_argument(
        "--save-hamiltonian", action="store_true",
        help="Compatibility flag: saves compact Pauli terms, not a CSR matrix."
    )
    p.add_argument("--save-pauli-terms", action="store_true")
    p.add_argument(
        "--validate-small", action="store_true",
        help="For dimension <= 65536, compare random Hx products with OpenFermion CSR."
    )
    p.add_argument(
        "--memory-fraction", type=float, default=0.80,
        help="Refuse an estimated eigsh workspace above this fraction of free VRAM."
    )
    p.add_argument("--force-memory", action="store_true")
    p.add_argument(
        "--threads-per-block", type=int, default=128,
        help="CUDA kernel threads per block."
    )
    return p.parse_args()


def import_dependencies() -> tuple[Any, Any, Any, Any]:
    try:
        import scipy.sparse as sp
    except ImportError as exc:
        raise RuntimeError("Missing SciPy: python -m pip install scipy") from exc

    try:
        import cupy as cp
        import cupyx.scipy.sparse.linalg as cpx_linalg
    except ImportError as exc:
        raise RuntimeError(
            "Missing CuPy CUDA 13 wheel: python -m pip install cupy-cuda13x"
        ) from exc

    try:
        from openfermion import count_qubits
        from openfermion.ops import QubitOperator
        from openfermion.transforms import jordan_wigner
    except ImportError as exc:
        raise RuntimeError(
            "Missing OpenFermion: python -m pip install openfermion"
        ) from exc

    return sp, cp, cpx_linalg, (count_qubits, QubitOperator, jordan_wigner)


def read_num_electrons_from_dat(path: Path) -> int:
    needle = "number of active, occupied, single-occupancy orbitals"
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if needle in line:
                try:
                    return int(line.rsplit("=", 1)[1].strip())
                except (IndexError, ValueError) as exc:
                    raise ValueError(f"Could not parse electron count from: {line.rstrip()}") from exc
    raise ValueError(f"Did not find '{needle}' in {path}")


def load_pickle(path: Path) -> Any:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("rb") as f:
        return pickle.load(f)


def infer_qubits(obj: Any, count_qubits: Any) -> int:
    for attr in ("n_qubits", "num_qubits"):
        value = getattr(obj, attr, None)
        if isinstance(value, (int, np.integer)):
            return int(value)
    return int(count_qubits(obj))


# -------------------------- Pauli term preparation --------------------------


@dataclass
class PreparedPauliTerms:
    n_qubits: int
    identity: complex
    group_x: np.ndarray       # uint64, one per unique X/Y flip mask
    group_offsets: np.ndarray # int32, len(groups)+1
    term_z: np.ndarray        # uint64, grouped by X mask
    term_coeff: np.ndarray    # complex128; includes i**nY phase
    raw_term_count: int

    @property
    def group_count(self) -> int:
        return int(self.group_x.size)

    @property
    def nonidentity_count(self) -> int:
        return int(self.term_z.size)


def prepare_pauli_terms(
    qubit_operator: Any,
    n_qubits: int,
    hermitian_tolerance: float = 1e-10,
) -> PreparedPauliTerms:
    """Convert OpenFermion terms to grouped bit masks.

    For a basis bitstring s, a Pauli term maps it to s xor xmask with phase
      (i ** number_of_Y) * (-1) ** popcount(s & zmask),
    where zmask contains both Y and Z locations.
    """
    identity = 0.0 + 0.0j
    grouped: dict[int, list[tuple[int, complex]]] = {}
    raw_count = len(qubit_operator.terms)

    max_supported = 63
    if n_qubits > max_supported:
        raise ValueError(
            f"This implementation currently supports <= {max_supported} qubits "
            "because basis states are encoded in uint64 masks."
        )

    for term, coefficient in qubit_operator.terms.items():
        coefficient = complex(coefficient)
        if not term:
            identity += coefficient
            continue

        xmask = 0
        zmask = 0
        ny = 0
        for qubit, action in term:
            q = int(qubit)
            if q < 0 or q >= n_qubits:
                raise ValueError(f"Term references qubit {q}, outside n_qubits={n_qubits}")
            bit = 1 << q
            if action == "X":
                xmask |= bit
            elif action == "Y":
                xmask |= bit
                zmask |= bit
                ny += 1
            elif action == "Z":
                zmask |= bit
            else:
                raise ValueError(f"Unsupported Pauli action: {action!r}")

        effective = coefficient * (1j ** ny)
        grouped.setdefault(xmask, []).append((zmask, effective))

    # Place diagonal xmask=0 first; otherwise sort by mask for reproducibility.
    masks = sorted(grouped, key=lambda m: (m != 0, m))
    group_x: list[int] = []
    offsets = [0]
    zs: list[int] = []
    coeffs: list[complex] = []

    for xmask in masks:
        group_x.append(xmask)
        # Combine any accidental duplicate z masks.
        merged: dict[int, complex] = {}
        for zmask, coeff in grouped[xmask]:
            merged[zmask] = merged.get(zmask, 0.0j) + coeff
        for zmask in sorted(merged):
            coeff = merged[zmask]
            if abs(coeff) > 1e-15:
                zs.append(zmask)
                coeffs.append(coeff)
        offsets.append(len(zs))

    # For a Hermitian QubitOperator the original Pauli coefficients should be
    # real.  effective coefficients can be imaginary because Y action phases
    # are folded into them.
    max_imag = max((abs(complex(c).imag) for c in qubit_operator.terms.values()), default=0.0)
    if max_imag > hermitian_tolerance:
        print(
            f"WARNING: maximum imaginary part of original Pauli coefficient is {max_imag:.3e}; "
            "verify Hermiticity.",
            flush=True,
        )

    return PreparedPauliTerms(
        n_qubits=n_qubits,
        identity=identity,
        group_x=np.asarray(group_x, dtype=np.uint64),
        group_offsets=np.asarray(offsets, dtype=np.int32),
        term_z=np.asarray(zs, dtype=np.uint64),
        term_coeff=np.asarray(coeffs, dtype=np.complex128),
        raw_term_count=raw_count,
    )


def save_compact_terms(path: Path, terms: PreparedPauliTerms) -> None:
    np.savez_compressed(
        path,
        n_qubits=np.asarray(terms.n_qubits, dtype=np.int32),
        identity=np.asarray(terms.identity, dtype=np.complex128),
        group_x=terms.group_x,
        group_offsets=terms.group_offsets,
        term_z=terms.term_z,
        term_coeff=terms.term_coeff,
    )


# ------------------------------- CUDA kernel --------------------------------


def cuda_source(dtype_name: str, fixed_weight: bool) -> str:
    if dtype_name == "complex128":
        ctype = "cuDoubleComplex"
        make = "make_cuDoubleComplex"
        cadd = "cuCadd"
        cmul = "cuCmul"
    else:
        ctype = "cuFloatComplex"
        make = "make_cuFloatComplex"
        cadd = "cuCaddf"
        cmul = "cuCmulf"

    common = f"""
#include <cuComplex.h>
extern "C" {{

__device__ __forceinline__ int parity64(unsigned long long x) {{
    return __popcll(x) & 1;
}}

__device__ __forceinline__ {ctype} negate_complex({ctype} z) {{
    return {make}(-z.x, -z.y);
}}
"""

    if not fixed_weight:
        body = f"""
__global__ void pauli_matvec_full(
    const {ctype}* __restrict__ x,
    {ctype}* __restrict__ y,
    const unsigned long long dimension,
    const unsigned long long* __restrict__ group_x,
    const int* __restrict__ group_offsets,
    const unsigned long long* __restrict__ term_z,
    const {ctype}* __restrict__ term_coeff,
    const int group_count,
    const {ctype} identity)
{{
    unsigned long long out =
        (unsigned long long)blockDim.x * blockIdx.x + threadIdx.x;
    if (out >= dimension) return;

    {ctype} acc = {cmul}(identity, x[out]);
    for (int g = 0; g < group_count; ++g) {{
        unsigned long long in_state = out ^ group_x[g];
        {ctype} amplitude = {make}(0.0, 0.0);
        int begin = group_offsets[g];
        int end = group_offsets[g + 1];
        for (int t = begin; t < end; ++t) {{
            {ctype} c = term_coeff[t];
            if (parity64(in_state & term_z[t])) c = negate_complex(c);
            amplitude = {cadd}(amplitude, c);
        }}
        acc = {cadd}(acc, {cmul}(amplitude, x[in_state]));
    }}
    y[out] = acc;
}}
}} // extern C
"""
        return common + body

    body = f"""
__device__ __forceinline__ unsigned long long choose_value(
    const unsigned long long* __restrict__ table,
    const int stride,
    const int n,
    const int k)
{{
    if (n < 0 || k < 0 || k > n) return 0ULL;
    return table[n * stride + k];
}}

// Colexicographic unranking: rank -> n-bit word with exactly weight bits.
__device__ __forceinline__ unsigned long long unrank_colex(
    unsigned long long rank,
    const int n,
    const int weight,
    const unsigned long long* __restrict__ binom,
    const int stride)
{{
    unsigned long long state = 0ULL;
    int x = n - 1;
    for (int j = weight; j >= 1; --j) {{
        while (x >= j && choose_value(binom, stride, x, j) > rank) --x;
        state |= (1ULL << x);
        rank -= choose_value(binom, stride, x, j);
        --x;
    }}
    return state;
}}

// n-bit word with fixed Hamming weight -> colexicographic rank.
__device__ __forceinline__ unsigned long long rank_colex(
    unsigned long long state,
    const int n,
    const unsigned long long* __restrict__ binom,
    const int stride)
{{
    unsigned long long rank = 0ULL;
    int j = 1;
    while (state) {{
        int pos = __ffsll((long long)state) - 1;
        rank += choose_value(binom, stride, pos, j);
        state &= (state - 1ULL);
        ++j;
    }}
    return rank;
}}

__global__ void pauli_matvec_sector(
    const {ctype}* __restrict__ x,
    {ctype}* __restrict__ y,
    const unsigned long long dimension,
    const int n_qubits,
    const int weight,
    const unsigned long long* __restrict__ binom,
    const int binom_stride,
    const unsigned long long* __restrict__ group_x,
    const int* __restrict__ group_offsets,
    const unsigned long long* __restrict__ term_z,
    const {ctype}* __restrict__ term_coeff,
    const int group_count,
    const {ctype} identity)
{{
    unsigned long long out_rank =
        (unsigned long long)blockDim.x * blockIdx.x + threadIdx.x;
    if (out_rank >= dimension) return;

    unsigned long long out_state =
        unrank_colex(out_rank, n_qubits, weight, binom, binom_stride);
    {ctype} acc = {cmul}(identity, x[out_rank]);

    for (int g = 0; g < group_count; ++g) {{
        unsigned long long in_state = out_state ^ group_x[g];
        if (__popcll(in_state) != weight) continue;

        {ctype} amplitude = {make}(0.0, 0.0);
        int begin = group_offsets[g];
        int end = group_offsets[g + 1];
        for (int t = begin; t < end; ++t) {{
            {ctype} c = term_coeff[t];
            if (parity64(in_state & term_z[t])) c = negate_complex(c);
            amplitude = {cadd}(amplitude, c);
        }}
        unsigned long long in_rank =
            rank_colex(in_state, n_qubits, binom, binom_stride);
        acc = {cadd}(acc, {cmul}(amplitude, x[in_rank]));
    }}
    y[out_rank] = acc;
}}
}} // extern C
"""
    return common + body


def make_binomial_table(n: int, k: int) -> np.ndarray:
    table = np.zeros((n + 1, k + 1), dtype=np.uint64)
    for i in range(n + 1):
        table[i, 0] = 1
        for j in range(1, min(i, k) + 1):
            table[i, j] = math.comb(i, j)
    return table


class PauliMatrixFreeOperator:
    def __init__(
        self,
        cp: Any,
        terms: PreparedPauliTerms,
        dtype_name: str,
        num_electrons: Optional[int],
        threads_per_block: int,
        progress_every: int,
    ) -> None:
        self.cp = cp
        self.terms = terms
        self.dtype_name = dtype_name
        self.dtype = cp.complex128 if dtype_name == "complex128" else cp.complex64
        self.num_electrons = num_electrons
        self.threads = threads_per_block
        self.progress_every = max(0, progress_every)
        self.matvec_calls = 0
        self.last_matvec_seconds = 0.0

        n = terms.n_qubits
        if num_electrons is None:
            self.dimension = 1 << n
            source = cuda_source(dtype_name, fixed_weight=False)
            self.kernel = cp.RawKernel(
                source,
                "pauli_matvec_full",
                options=("--std=c++17",),
            )
            self.binom_gpu = None
            self.binom_stride = 0
        else:
            if not 0 <= num_electrons <= n:
                raise ValueError("--num-electrons must satisfy 0 <= N <= n_qubits")
            self.dimension = math.comb(n, num_electrons)
            table = make_binomial_table(n, num_electrons)
            self.binom_gpu = cp.asarray(table)
            self.binom_stride = table.shape[1]
            source = cuda_source(dtype_name, fixed_weight=True)
            self.kernel = cp.RawKernel(
                source,
                "pauli_matvec_sector",
                options=("--std=c++17",),
            )

        coeff_dtype = np.complex128 if dtype_name == "complex128" else np.complex64
        self.group_x_gpu = cp.asarray(terms.group_x)
        self.group_offsets_gpu = cp.asarray(terms.group_offsets)
        self.term_z_gpu = cp.asarray(terms.term_z)
        self.term_coeff_gpu = cp.asarray(terms.term_coeff.astype(coeff_dtype, copy=False))
        self.identity_gpu_scalar = coeff_dtype(terms.identity)
        self.blocks = (self.dimension + self.threads - 1) // self.threads

    def matvec(self, x: Any) -> Any:
        cp = self.cp
        x = cp.asarray(x, dtype=self.dtype).reshape(-1)
        if x.size != self.dimension:
            raise ValueError(f"matvec expected length {self.dimension}, got {x.size}")
        y = cp.empty_like(x)

        timed = self.progress_every > 0 and (self.matvec_calls % self.progress_every == 0)
        if timed:
            cp.cuda.Stream.null.synchronize()
            t0 = time.perf_counter()

        if self.num_electrons is None:
            self.kernel(
                (self.blocks,),
                (self.threads,),
                (
                    x,
                    y,
                    np.uint64(self.dimension),
                    self.group_x_gpu,
                    self.group_offsets_gpu,
                    self.term_z_gpu,
                    self.term_coeff_gpu,
                    np.int32(self.terms.group_count),
                    self.identity_gpu_scalar,
                ),
            )
        else:
            self.kernel(
                (self.blocks,),
                (self.threads,),
                (
                    x,
                    y,
                    np.uint64(self.dimension),
                    np.int32(self.terms.n_qubits),
                    np.int32(self.num_electrons),
                    self.binom_gpu,
                    np.int32(self.binom_stride),
                    self.group_x_gpu,
                    self.group_offsets_gpu,
                    self.term_z_gpu,
                    self.term_coeff_gpu,
                    np.int32(self.terms.group_count),
                    self.identity_gpu_scalar,
                ),
            )

        self.matvec_calls += 1
        if timed:
            cp.cuda.Stream.null.synchronize()
            self.last_matvec_seconds = time.perf_counter() - t0
            print(
                f"[Hx {self.matvec_calls:6d}] {self.last_matvec_seconds:.3f} s",
                flush=True,
            )
        return y

    def matmat(self, x: Any) -> Any:
        cp = self.cp
        x = cp.asarray(x, dtype=self.dtype)
        if x.ndim == 1:
            return self.matvec(x)
        columns = [self.matvec(x[:, i]) for i in range(x.shape[1])]
        return cp.stack(columns, axis=1)


# ------------------------- validation and output ----------------------------


def validate_small_operator(
    obj: Any,
    qubit_operator: Any,
    operator: PauliMatrixFreeOperator,
    cp: Any,
    num_electrons: Optional[int],
) -> float:
    """Compare matrix-free Hx to OpenFermion's explicit sparse matrix.

    Fixed-sector validation builds the full matrix only for small full-space
    dimensions and extracts the selected Hamming-weight rows/columns.
    """
    if (1 << operator.terms.n_qubits) > 65536:
        raise ValueError("--validate-small is limited to full dimension <= 65536")

    from openfermion.linalg import get_sparse_operator

    h_full = get_sparse_operator(qubit_operator, n_qubits=operator.terms.n_qubits).tocsr()
    rng = np.random.default_rng(12345)

    if num_electrons is None:
        h_ref = h_full
    else:
        states = np.asarray(
            [i for i in range(1 << operator.terms.n_qubits) if i.bit_count() == num_electrons],
            dtype=np.int64,
        )
        # Reorder to the colex rank used by the CUDA kernel.
        def colex_rank(state: int) -> int:
            rank = 0
            j = 1
            for pos in range(operator.terms.n_qubits):
                if (state >> pos) & 1:
                    rank += math.comb(pos, j)
                    j += 1
            return rank
        order = np.argsort([colex_rank(int(s)) for s in states])
        states = states[order]
        h_ref = h_full[states][:, states]

    x = rng.normal(size=operator.dimension) + 1j * rng.normal(size=operator.dimension)
    x = x.astype(np.complex128 if operator.dtype_name == "complex128" else np.complex64)
    y_ref = h_ref @ x
    y_gpu = cp.asnumpy(operator.matvec(cp.asarray(x)))
    relative = float(np.linalg.norm(y_ref - y_gpu) / max(np.linalg.norm(y_ref), 1e-300))
    print(f"Matrix-free validation relative Hx error: {relative:.6e}", flush=True)
    return relative


def estimate_workspace_bytes(
    dimension: int,
    itemsize: int,
    k: int,
    ncv: int,
    values_only: bool,
) -> int:
    # ARPACK/CuPy implementation details can change.  This deliberately uses a
    # conservative count rather than claiming exact allocation behavior.
    vector_count = ncv + 10 + (0 if values_only else k)
    return dimension * itemsize * vector_count


def save_vectors_incrementally(path: Path, vectors_gpu: Any, cp: Any) -> None:
    shape = tuple(int(x) for x in vectors_gpu.shape)
    mmap = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.dtype(vectors_gpu.dtype.name), shape=shape
    )
    for i in range(shape[1]):
        mmap[:, i] = cp.asnumpy(vectors_gpu[:, i])
        mmap.flush()
        print(f"Saved eigenvector {i+1}/{shape[1]}", flush=True)
    del mmap


def compute_residuals(
    operator: PauliMatrixFreeOperator,
    eigenvalues_gpu: Any,
    eigenvectors_gpu: Any,
    checks: int,
    cp: Any,
) -> np.ndarray:
    if eigenvectors_gpu is None or checks <= 0:
        return np.empty(0, dtype=np.float64)
    count = min(checks, int(eigenvalues_gpu.size))
    residuals = np.empty(count, dtype=np.float64)
    for i in range(count):
        v = eigenvectors_gpu[:, i]
        hv = operator.matvec(v)
        numerator = cp.linalg.norm(hv - eigenvalues_gpu[i] * v)
        denominator = cp.maximum(cp.linalg.norm(hv), cp.asarray(1e-300))
        residuals[i] = float((numerator / denominator).get())
        print(
            f"Residual eigenpair {i}: {residuals[i]:.3e}; "
            f"energy={float(cp.real(eigenvalues_gpu[i]).get()):.16g}",
            flush=True,
        )
    return residuals


def main() -> int:
    args = parse_args()
    timings = Timings()
    sp, cp, cpx_linalg, of_deps = import_dependencies()
    count_qubits, QubitOperator, jordan_wigner = of_deps

    if args.dat_file is not None and args.num_electrons is None and not args.full_space:
        args.num_electrons = read_num_electrons_from_dat(args.dat_file)
        print(f"Read --num-electrons {args.num_electrons} from {args.dat_file}", flush=True)
    if args.full_space:
        args.num_electrons = None

    cp.cuda.Device(args.gpu).use()
    props = cp.cuda.runtime.getDeviceProperties(args.gpu)
    gpu_name = props["name"]
    if isinstance(gpu_name, bytes):
        gpu_name = gpu_name.decode(errors="replace")
    free_vram, total_vram = cp.cuda.runtime.memGetInfo()
    print(f"GPU {args.gpu}: {gpu_name}", flush=True)
    print(f"CuPy: {cp.__version__}", flush=True)
    print(f"CUDA runtime seen by CuPy: {cp.cuda.runtime.runtimeGetVersion()}", flush=True)
    print(f"VRAM: {human_bytes(free_vram)} free / {human_bytes(total_vram)} total", flush=True)

    obj = load_pickle(args.input)
    timings.mark("Input loaded")
    print(f"Input object type: {type(obj)!r}", flush=True)

    n_qubits = args.n_qubits if args.n_qubits is not None else infer_qubits(obj, count_qubits)
    print(f"Number of qubits: {n_qubits}", flush=True)

    if isinstance(obj, QubitOperator):
        qop = obj
    else:
        print("Applying OpenFermion Jordan-Wigner transform (compact Pauli representation)", flush=True)
        qop = jordan_wigner(obj)
    del obj
    timings.mark("Pauli representation ready")

    terms = prepare_pauli_terms(qop, n_qubits)
    print(f"Raw Pauli terms including identity: {terms.raw_term_count:,}", flush=True)
    print(f"Nonidentity Pauli terms after merging: {terms.nonidentity_count:,}", flush=True)
    print(f"Unique X/Y flip-mask groups: {terms.group_count:,}", flush=True)
    print(f"Identity coefficient: {terms.identity!r}", flush=True)

    if args.check_hermitian:
        max_imag = max((abs(complex(c).imag) for c in qop.terms.values()), default=0.0)
        print(f"Maximum imaginary part of Pauli coefficients: {max_imag:.6e}", flush=True)
        if max_imag > 1e-9:
            raise ValueError("Pauli coefficients are not real within Hermitian tolerance")
    del qop

    if args.num_electrons is None:
        dimension = 1 << n_qubits
        basis_label = "full_qubit_space"
    else:
        dimension = math.comb(n_qubits, args.num_electrons)
        basis_label = "fixed_particle_number"

    itemsize = np.dtype(args.dtype).itemsize
    print(f"Basis mode: {basis_label}", flush=True)
    if args.num_electrons is not None:
        print(f"Active electrons/Hamming weight: {args.num_electrons}", flush=True)
    print(f"Operator dimension: {dimension:,}", flush=True)
    print(f"One {args.dtype} vector: {human_bytes(dimension * itemsize)}", flush=True)

    k = args.num_eigenpairs
    if not 1 <= k < dimension:
        raise ValueError(f"-k must satisfy 1 <= k < dimension ({dimension:,})")
    ncv = args.ncv if args.ncv is not None else min(dimension, max(2 * k + 1, 20))
    if not k < ncv <= dimension:
        raise ValueError(f"--ncv must satisfy k < ncv <= dimension; got k={k}, ncv={ncv}")

    workspace = estimate_workspace_bytes(dimension, itemsize, k, ncv, args.values_only)
    print(
        f"Conservative eigsh vector-workspace estimate: {human_bytes(workspace)} "
        f"(k={k}, ncv={ncv})",
        flush=True,
    )
    if workspace > args.memory_fraction * free_vram and not args.force_memory:
        raise MemoryError(
            f"Estimated matrix-free eigsh workspace {human_bytes(workspace)} exceeds "
            f"{args.memory_fraction:.0%} of free VRAM ({human_bytes(free_vram)}). "
            "Reduce k/ncv, use --num-electrons, use --values-only, or override with "
            "--force-memory at your own risk."
        )

    if args.inspect_only:
        print("Inspection complete; no eigensolver was run.", flush=True)
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.save_hamiltonian or args.save_pauli_terms:
        save_compact_terms(args.output_dir / "pauli_terms.npz", terms)
        timings.mark("Compact Pauli terms saved")

    operator = PauliMatrixFreeOperator(
        cp=cp,
        terms=terms,
        dtype_name=args.dtype,
        num_electrons=args.num_electrons,
        threads_per_block=args.threads_per_block,
        progress_every=args.progress_every,
    )
    timings.mark("GPU matrix-free operator initialized")

    # Force JIT compilation and benchmark one Hx before entering eigsh.
    rng = cp.random.RandomState(12345)
    trial = rng.standard_normal(operator.dimension, dtype=cp.float64)
    trial = trial.astype(operator.dtype)
    trial /= cp.linalg.norm(trial)
    print("Compiling/benchmarking one matrix-free Hx operation...", flush=True)
    _ = operator.matvec(trial)
    cp.cuda.Stream.null.synchronize()
    del trial, _
    cp.get_default_memory_pool().free_all_blocks()
    timings.mark("First matrix-free Hx completed")

    validation_error = None
    if args.validate_small:
        # Reload because the original object was deliberately freed.
        validation_obj = load_pickle(args.input)
        validation_qop = validation_obj if isinstance(validation_obj, QubitOperator) else jordan_wigner(validation_obj)
        validation_error = validate_small_operator(
            validation_obj, validation_qop, operator, cp, args.num_electrons
        )
        del validation_obj, validation_qop
        timings.mark("Small-system validation completed")

    linear_operator = cpx_linalg.LinearOperator(
        shape=(dimension, dimension),
        matvec=operator.matvec,
        rmatvec=operator.matvec,
        matmat=operator.matmat,
        dtype=operator.dtype,
    )

    print(
        f"Running matrix-free GPU eigsh: k={k}, ncv={ncv}, which='SA', "
        f"tol={args.tol:g}, maxiter={args.maxiter}",
        flush=True,
    )
    v0 = cp.ones(dimension, dtype=operator.dtype)
    v0 /= cp.linalg.norm(v0)

    with SolverHeartbeat(operator, args.heartbeat):
        result = cpx_linalg.eigsh(
            linear_operator,
            k=k,
            which="SA",
            v0=v0,
            ncv=ncv,
            maxiter=args.maxiter,
            tol=args.tol,
            return_eigenvectors=not args.values_only,
        )
        cp.cuda.Stream.null.synchronize()
    del v0
    timings.mark("Matrix-free eigsh completed")

    if args.values_only:
        eigenvalues_gpu = result
        eigenvectors_gpu = None
    else:
        eigenvalues_gpu, eigenvectors_gpu = result

    order = cp.argsort(cp.real(eigenvalues_gpu))
    eigenvalues_gpu = eigenvalues_gpu[order]
    if eigenvectors_gpu is not None:
        eigenvectors_gpu = eigenvectors_gpu[:, order]

    residuals = compute_residuals(
        operator, eigenvalues_gpu, eigenvectors_gpu, args.residual_checks, cp
    )
    timings.mark("Residual checks completed")

    eigenvalues = cp.asnumpy(eigenvalues_gpu)
    np.save(args.output_dir / "eigenvalues.npy", eigenvalues)
    np.savetxt(
        args.output_dir / "lowest_energies.txt",
        np.column_stack((np.arange(eigenvalues.size), np.real(eigenvalues))),
        header="state_index energy",
        fmt=["%d", "%.17e"],
    )
    sp.save_npz(
        args.output_dir / "diagonal_sparse.npz",
        sp.diags(eigenvalues, offsets=0, format="csr"),
    )
    np.save(args.output_dir / "residuals.npy", residuals)

    if eigenvectors_gpu is not None:
        save_vectors_incrementally(
            args.output_dir / "eigenvectors.npy", eigenvectors_gpu, cp
        )
    timings.mark("Numerical outputs saved")

    metadata = {
        "input": str(args.input),
        "dat_file": None if args.dat_file is None else str(args.dat_file),
        "gpu_index": args.gpu,
        "gpu_name": str(gpu_name),
        "cupy_version": cp.__version__,
        "cuda_runtime": int(cp.cuda.runtime.runtimeGetVersion()),
        "n_qubits": n_qubits,
        "num_electrons": args.num_electrons,
        "basis_mode": basis_label,
        "dimension": int(dimension),
        "raw_pauli_terms": terms.raw_term_count,
        "nonidentity_pauli_terms": terms.nonidentity_count,
        "flip_mask_groups": terms.group_count,
        "dtype": args.dtype,
        "k": k,
        "ncv": ncv,
        "tol": args.tol,
        "maxiter": args.maxiter,
        "values_only": args.values_only,
        "matvec_calls": operator.matvec_calls,
        "last_matvec_seconds": operator.last_matvec_seconds,
        "workspace_estimate_bytes": int(workspace),
        "validation_relative_error": validation_error,
        "lowest_energy": float(np.real(eigenvalues[0])),
        "highest_returned_energy": float(np.real(eigenvalues[-1])),
    }
    with (args.output_dir / "metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    with (args.output_dir / "timing.json").open("w", encoding="utf-8") as f:
        json.dump(timings.as_dict(), f, indent=2)

    timings.mark("Metadata saved")
    print("", flush=True)
    print(f"Returned {eigenvalues.size} lowest energy value(s)", flush=True)
    for i, value in enumerate(eigenvalues):
        print(f"E[{i:4d}] = {value:.17g}", flush=True)
    print(f"Output directory: {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr, flush=True)
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr, flush=True)
        traceback.print_exc()
        raise SystemExit(1)
