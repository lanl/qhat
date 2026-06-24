#!/usr/bin/env python
# coding: utf-8

import math
import h5py
from typing import Optional
import numpy as np






class ResourceEstimator:
    """Estimate Toffoli gate count and logical qubit count for the Bloch--UPAW block encoding.

    Implements the resource formulas from Appendix E of the paper:
        - Toffoli count per block-encoding step (Eq. 98)
        - Total logical qubit count (Eq. 99)

    with supporting bit-width definitions from Eqs. (90)--(91). The total
    number of LCU labels is L = 2·(Npw + P)·Nk, where P = Σ_a na·(na+1)/2
    counts the upper-triangular partial-wave pairs across all atoms.

    The ``qroam='optimum'`` mode (default) minimises the QROAM batch-size
    parameters analytically, yielding the scaling in Eq. (41).

    Typical usage::

        eps_chem = 1.6e-3  # chemical accuracy in Hartree (1 kcal/mol)
        eps_qpe = eps_chem / 5  # QPE share of 5-way error budget

        est = ResourceEstimator.from_hdf5("integrals.h5")
        toffolis = est.toffoli_count_per_be(Rl=R_avg, R0=R0)
        qubits = est.total_qubits(Rl=R_avg, R0=R0,
                                   lam=lam, eps_qpe=eps_qpe)

    Precision parameters br, N1, N2, B have sensible defaults (see
    ``toffoli_count_per_be`` docstring for guidance on when to change them).

    Args:
        Nk: Number of k-points.
        Npw: Number of plane waves (G ≠ 0).
        P: Total number of upper-triangular partial-wave pairs.
        Nb: Number of bands per k-point (optional, inferred from HDF5).
        eta: QROAM state-preparation parameter — the smallest integer such
            that 2^eta divides L.  Conventionally hardcoded to 10 following
            Ivanov & Rubin (default 10).
    """

    # ---------------- instance fields loaded from HDF5 ----------------
    def __init__(self, Nk: int, Npw: int, P: int, Nb: Optional[int] = None, eta: int = 10):
        if Nk <= 0:
            raise ValueError("Nk must be positive.")
        if Npw < 0 or P < 0:
            raise ValueError("Npw and P must be nonnegative.")
        self.Nk   = int(Nk)
        self.Npw  = int(Npw)
        self.P    = int(P)
        self.L    = int(2 * (self.Npw + self.P) * self.Nk)
        self.Nb   = (int(Nb) if Nb is not None else None)
        self.eta  = int(eta)

    # ---------------- builder from HDF5 written by export_hdf5 ----------------
    @staticmethod
    def _infer_Nb_from_one_body(h5: "h5py.File") -> Optional[int]:
        """Try to deduce Nb from one-body dataset shapes."""
        if "one_body" not in h5:
            return None
        g = h5["one_body"]
        # Preferred dataset name used by export_hdf5 in earlier code:
        candidates = []
        if "H_kikj" in g:
            candidates.append(g["H_kikj"])
        # Otherwise, try any dataset under one_body:
        for name, ds in g.items():
            if isinstance(ds, h5py.Dataset) and ds not in candidates:
                candidates.append(ds)

        for ds in candidates:
            shp = ds.shape
            # Common cases:
            # (Nk, Nb, Nk, Nb)
            if len(shp) == 4 and shp[0] == shp[2] and shp[1] == shp[3]:
                return int(shp[1])
            # (Nk, Nb, Nb)
            if len(shp) == 3 and shp[1] == shp[2]:
                return int(shp[1])
            # (Nb, Nb, Nk)
            if len(shp) == 3 and shp[0] == shp[1]:
                return int(shp[0])
            # (Nb, Nb) for single-k export
            if len(shp) == 2 and shp[0] == shp[1]:
                return int(shp[0])
        return None

    @classmethod
    def from_hdf5(cls, filepath: str, *, kmesh_region: str = "bz") -> "ResourceEstimator":
        if kmesh_region not in {"bz", "ibz"}:
            raise ValueError("kmesh_region must be 'bz' or 'ibz'.")

        with h5py.File(filepath, "r") as h5:
            # --- Nk ---
            Nk = None
            try:
                Nk = int(h5["kmesh"][kmesh_region]["reduced"].shape[0])
            except Exception:
                try:
                    Nk = int(h5["kmesh"][kmesh_region]["cart_1_perA"].shape[0])
                except Exception:
                    other = "bz" if kmesh_region == "ibz" else "ibz"
                    for name in ("reduced", "cart_1_perA"):
                        try:
                            Nk = int(h5["kmesh"][other][name].shape[0]); break
                        except Exception:
                            pass
            if Nk is None:
                raise KeyError("Could not determine Nk from HDF5 kmesh.")

            # --- Npw ---
            try:
                Npw = int(h5["Npw"][()])
            except Exception as e:
                raise KeyError("Dataset 'Npw' not found in HDF5.") from e

            # --- P from C_tensor (preferred) or D_tensor fallback ---
            P = 0
            if "C_tensor" in h5:
                for _, ds in h5["C_tensor"].items():
                    if not isinstance(ds, h5py.Dataset):
                        continue
                    na = int(ds.shape[0])
                    if ds.shape != (na, na, na, na):
                        raise ValueError(f"C_tensor dataset has unexpected shape {ds.shape}.")
                    P += na * (na + 1) // 2
            elif "D_tensor" in h5:
                seen = False
                for _, ds in h5["D_tensor"].items():
                    if not isinstance(ds, h5py.Dataset):
                        continue
                    if len(ds.shape) < 2 or ds.shape[-1] != ds.shape[-2]:
                        raise ValueError(f"D_tensor dataset has unexpected shape {ds.shape}.")
                    na = int(ds.shape[-1])
                    P += na * (na + 1) // 2
                    seen = True
                if not seen:
                    raise KeyError("D_tensor group present but empty.")
            else:
                raise KeyError("Neither 'C_tensor' nor 'D_tensor' present to infer P.")

            # eta is a QROAM state-preparation parameter (2-adic valuation
            # of L), conventionally hardcoded to 10 following Ivanov & Rubin.
            # It is NOT the electron count (n_electrons is separate metadata).

            # --- Nb from one-body (robust to a few layout variants)
            Nb = cls._infer_Nb_from_one_body(h5)

        return cls(Nk=Nk, Npw=Npw, P=P, Nb=Nb)

    # ---------------- bit-width helpers (Eqs. 90–91) ----------------
    @staticmethod
    def _bits(x: float) -> int:
        if x <= 1:
            return 1
        return math.ceil(math.log2(x))

    @staticmethod
    def _ceil_log2_nonneg(x: float) -> int:
        if x <= 1:
            return 0
        return math.ceil(math.log2(x))

    @staticmethod
    def _nL(L: int) -> int:
        return ResourceEstimator._bits(L)

    @staticmethod
    def _nk(Nk: int) -> int:
        return ResourceEstimator._bits(Nk)

    @staticmethod
    def _nR(Rl: int) -> int:
        return ResourceEstimator._bits(Rl)

    @staticmethod
    def _nLR(L: int, Rl: int, R0: int, Nk: int) -> int:
        val = L * Rl + (Nk * R0) 
        return ResourceEstimator._bits(val)

    @staticmethod
    def _bo(nk: int, nR: int, nLR: int, br: int) -> int:
        return nk + nR + nLR + br + 2

    @staticmethod
    def _bp1(nL: int, N1: int) -> int:
        return nL + N1

    @staticmethod
    def _bp2(nL: int, N2: int) -> int:
        return nL + N2

    @staticmethod
    def _ceil_ratio(numer: float, denom: float) -> int:
        return math.ceil(numer / float(denom))

    # ---------------- internal cores with optional bo override ----------------
    @classmethod
    def _toffoli_core(cls,
                      L, kp1, kp1_prime, ko, ko_prime, kp2, kp2_prime, kr, kr_prime,
                      br, Rl, R0, Nb, Nk, N1, N2, eta, B,
                      bo_override: float | None = None) -> int:
        """
        Eq. 98: Toffoli count per block-encoding step.

        QROAM batch-size parameters (each minimizes the cost of one QROAM oracle;
        the 'prime' variant is the cost of the uncompute/inverse call):
          kp1, kp1_prime  — PREPARE oracle reading bp1-bit phase angles over the
                            L+1 label table (one-body + select amplitudes)
          ko,  ko_prime   — OUTPUT oracle reading bo-bit coefficients over the
                            L+1 label table (rotation-angle register load)
          kp2, kp2_prime  — PREPARE oracle reading bp2-bit phase angles over the
                            LR-sized two-body amplitude table (rank-weighted)
          kr,  kr_prime   — rotation oracle reading (4·Nb·B + nk)-bit Givens
                            angles over the LR-sized table

        Bit-width shorthands (Eqs. 90–91):
          nL  = ceil(log2(L))
          nk  = ceil(log2(Nk))
          nR  = ceil(log2(Rl))
          nLR = ceil(log2(L·Rl + Nk·R0/2))
          bo  = nk + nR + nLR + br + 2   (output register width)
          bp1 = nL + N1                  (PREPARE phase-angle width, one-body)
          bp2 = nL + N2                  (PREPARE phase-angle width, two-body)
        """
        nL   = cls._nL(L)
        nk   = cls._nk(Nk)
        nR   = cls._nR(Rl)
        nLR  = cls._nLR(L, Rl, R0, Nk)
        bo   = (bo_override if bo_override is not None else cls._bo(nk, nR, nLR, br))
        bp1  = cls._bp1(nL, N1)
        bp2  = cls._bp2(nL, N2)

        # LR = L·Rl is the dominant table size for two-body QROAM.
        # half_term = Nk·R0/2 is the one-body contribution to that table.
        LR        = L * Rl
        half_term = (Nk * R0) 

        # ---- QROAM read costs over the L+1 label table (Eq. 98, first block) ----
        # Each term ceil((L+1)/k) is the Toffoli cost of one QROAM call with
        # batch size k over a table of size L+1.
        t  = cls._ceil_ratio(L + 1, kp1)         # PREPARE (forward)
        t += cls._ceil_ratio(L + 1, kp1_prime)   # PREPARE (uncompute)
        t += cls._ceil_ratio(L + 1, ko)           # OUTPUT  (forward)
        t += cls._ceil_ratio(L + 1, ko_prime)     # OUTPUT  (uncompute)

        # ---- QROAM read costs over the LR-sized two-body table (Eq. 98, second block) ----
        # Each oracle appears twice: once for the full table (LR + half_term)
        # and once for the rank-only portion (LR), giving the two ceil(·/k) terms.
        t += cls._ceil_ratio(LR + half_term, kp2)        # PREPARE two-body (forward, full)
        t += cls._ceil_ratio(LR,             kp2)        # PREPARE two-body (forward, rank)

        t += cls._ceil_ratio(LR + half_term, kr)         # rotation oracle (forward, full)
        t += cls._ceil_ratio(LR,             kr)         # rotation oracle (forward, rank)

        t += cls._ceil_ratio(LR + half_term, kr_prime)   # rotation oracle (uncompute, full)
        t += cls._ceil_ratio(LR,             kr_prime)   # rotation oracle (uncompute, rank)

        t += cls._ceil_ratio(LR + half_term, kp2_prime)  # PREPARE two-body (uncompute, full)
        t += cls._ceil_ratio(LR,             kp2_prime)  # PREPARE two-body (uncompute, rank)

        # ---- QROAM ancilla costs (Eq. 98, third block) ----
        # Each QROAM oracle with batch size k stores k output copies; the ancillae
        # for the forward calls are cleaned up at cost (bitwidth)·(k-1) Toffolis.
        # The uncompute calls cost only k Toffolis regardless of bitwidth.
        t += bp1 * (kp1 - 1)               # PREPARE ancillae (forward)
        t += bo  * (ko  - 1)               # OUTPUT  ancillae (forward)
        t += 2 * bp2 * (kp2 - 1)           # PREPARE two-body ancillae (forward, factor 2 for both halves)
        t += (4 * Nb * B + nk) * (kr - 1)  # rotation ancillae (forward, angle bitwidth = 4·Nb·B + nk)
        t += kp1_prime + ko_prime + 2 * kr_prime + 2 * kp2_prime  # all uncompute ancillae

        # ---- Arithmetic costs (Eq. 98, fourth block) ----
        # Fixed-point arithmetic for Hamiltonian coefficient loading, Givens
        # rotations, and phase kickback: depends on bit-widths and system size.
        t += (9 * nL + 34 * nR + 8 * nLR + 3 * N1 + 6 * N2
              + 12 * br - 6 * eta
              + 16 * Nb * B - 32 * Nb + 6 * R0 * Nk
              + 12 * nk - 43)
        return int(t)

    @classmethod
    def _total_qubits_core(cls,
                           L, kr,
                           br, Rl, R0, Nb, Nk, N1, N2, eta, B,
                           I: Optional[int] = None, lam: Optional[float] = None, eps_qpe: Optional[float] = None,
                           bo_override: float | None = None) -> int:
        """
        Eq. 99: total logical qubit count.

        Register layout (each term in the sum below):
          Nk·R0 + Nb          — system register: k-point (Nk·R0 qubits) + band index (Nb)
          2·nL + nR            — label registers: two L-index registers + one rank register
          2·N1 + N2 + B + bo + bp2 — coefficient/rotation-angle registers
          2·kr·Nb·B            — QROAM rotation ancillae (kr output copies × Nb·B angle bits × 2 halves)
          2·ceil(log2(m))      — QROAM output ancillae for the rotation oracle (m = table size / kr)
          2·ceil(log2(I+1))    — phase estimation register (I = ceil(π·λ/ε_QPE) QPE iterations)
          9                    — miscellaneous ancillae (flag qubits, comparison bits, etc.)
        """
        nL   = cls._nL(L)
        nk   = cls._nk(Nk)
        nR   = cls._nR(Rl)
        nLR  = cls._nLR(L, Rl, R0, Nk)
        bo   = (bo_override if bo_override is not None else cls._bo(nk, nR, nLR, br))
        bp2  = cls._bp2(nL, N2)
        LR        = L * Rl
        half_term = (Nk * R0) 

        # Number of QPE iterations: I = ceil(π·λ/ε_QPE); if not provided, omit QPE register.
        if I is None:
            if lam is not None and eps_qpe is not None and eps_qpe > 0:
                I = math.ceil(math.pi * lam / eps_qpe)
            else:
                I = 0

        # m = ceil((LR + half_term) / kr) is the effective table size per QROAM output copy.
        m = cls._ceil_ratio(LR + half_term, kr)

        qubits = (
            # System register
            2 * Nk * R0 + Nb
            # Label and coefficient registers
            + 2 * nL + nR + 2 * N1 + N2 + B + bo + bp2
            # QROAM rotation ancillae
            + 2 * kr * Nb * B
            # QROAM output ancillae (log2 of QROAM output table depth)
            + 2 * cls._bits(m)
            # Phase estimation register
            + 2 * cls._ceil_log2_nonneg(I + 1)
            # Miscellaneous ancillae
            + 9
        )
        return int(qubits)

    # ---------------- instance APIs with qroam ----------------
    def _opt_params(self, *, br, Rl, R0, Nb: Optional[int], N1, N2, eta, B):
        """
        Compute the 'optimum' parameters (floats allowed).
        Defaults Nb to self.Nb when not provided.
        """
        if Nb is None:
            if self.Nb is None:
                raise ValueError("Nb not provided and not found in HDF5.")
            Nb = self.Nb

        L   = float(self.L)
        Nk  = float(self.Nk)
        Rl  = float(Rl)
        Nb  = float(Nb)
        brf = float(br)

        nL  = self._nL(self.L)
        bp1 = nL + int(N1)
        bp2 = nL + int(N2)

        bo_opt = (np.log2(Nk) + np.log2(Rl) + np.log2(L * Rl + (Nk * R0) / 2.0) + brf + 2.0)

        kp1        = np.sqrt((L + 1.0) / bp1)
        kp1_prime  = np.sqrt(L + 1.0)
        kp2        = np.sqrt((2.0 * L * Rl + (R0 * Nk) ) / (2.0 * bp2))
        kp2_prime  = np.sqrt((2.0 * L * Rl + (R0 * Nk) ) / 2.0)
        ko         = np.sqrt((L + 1.0) / bo_opt)
        ko_prime   = np.sqrt(L + 1.0)
        kr         = np.sqrt((2.0 * L * Rl + (R0 * Nk) ) / (4.0 * Nb * B + np.log2(Nk)))
        kr_prime   = np.sqrt((2.0 * L * Rl + (R0 * Nk) ) / 2.0)

        return dict(bo=bo_opt, bp1=bp1, bp2=bp2,
                    kp1=kp1, kp1_prime=kp1_prime,
                    ko=ko, ko_prime=ko_prime,
                    kp2=kp2, kp2_prime=kp2_prime,
                    kr=kr, kr_prime=kr_prime)

    def toffoli_count_per_be(
        self,
        *,
        Rl, R0, eta=None,
        br: int = 20,
        N1: int = 10,
        N2: int = 10,
        B: int = 20,
        Nb: Optional[int] = None,
        qroam: str = "optimum",
        kp1: float | None = None, kp1_prime: float | None = None,
        ko: float | None = None,  ko_prime: float | None = None,
        kp2: float | None = None, kp2_prime: float | None = None,
        kr: float | None = None,  kr_prime: float | None = None
    ) -> int:
        """Compute the Toffoli gate count per block-encoding step (Eq. 98).

        System-specific parameters:
            Rl: Average two-body rank R^(ℓ≠0), from OneNormCalculator output.
            R0: One-body rank R^(0), from OneNormCalculator output.
            eta: QROAM parameter (2-adic valuation of L, default 10).

        Precision parameters (have defaults):
            br: Bit precision for the rotation-angle register (default 20).
                Controls the accuracy of Hamiltonian coefficient loading into
                quantum registers.  20 bits gives ~6 decimal digits of
                precision, which is well below chemical accuracy (1.6 mHa
                = 1 kcal/mol) for typical systems.  Increase to 25--30
                for sub-μHa precision studies; decrease to 10--15 for
                quick estimates.
            N1: Phase-angle bit precision for the one-body PREPARE oracle
                (default 10).  Determines how accurately the LCU
                amplitudes for kinetic + local-potential terms are encoded.
                10 bits (~3 decimal digits) is sufficient for the one-body
                terms, which are typically smaller than two-body terms.
                Rarely needs adjustment.
            N2: Phase-angle bit precision for the two-body PREPARE oracle
                (default 10).  Same role as N1 but for the Coulomb
                integrals.  10 bits is sufficient for typical accuracy
                targets (ε_chem ≈ 1.6 mHa).  Increase for high-precision studies.
            B: Givens rotation bit precision (default 20).  Controls the
                accuracy of the Givens-rotation network that maps between
                the plane-wave and partial-wave bases.  20 bits matches the
                rotation-angle register precision (br=20).  Should generally
                be kept equal to br; change both together.
            Nb: Number of bands per k-point (defaults to self.Nb from HDF5).
            qroam: 'optimum' (default) or 'non-optimum'.
            kp1, kp1_prime, ...: QROAM batch sizes (required for non-optimum).

        Returns:
            Integer Toffoli count.
        """
        if eta is None:
            eta = self.eta

        if qroam.lower() == "optimum":
            p = self._opt_params(br=br, Rl=Rl, R0=R0, Nb=Nb, N1=N1, N2=N2, eta=eta, B=B)
            Nb_eff = self.Nb if Nb is None else Nb
            return self._toffoli_core(
                self.L,
                p["kp1"], p["kp1_prime"], p["ko"], p["ko_prime"],
                p["kp2"], p["kp2_prime"], p["kr"], p["kr_prime"],
                br, Rl, R0, Nb_eff, self.Nk, N1, N2, eta, B,
                bo_override=p["bo"]
            )
        elif qroam.lower() == "non-optimum":
            if Nb is None:
                if self.Nb is None:
                    raise ValueError("Provide Nb or ensure it is stored from HDF5.")
                Nb = self.Nb
            needed = dict(kp1=kp1, kp1_prime=kp1_prime, ko=ko, ko_prime=ko_prime,
                          kp2=kp2, kp2_prime=kp2_prime, kr=kr, kr_prime=kr_prime)
            missing = [k for k, v in needed.items() if v is None]
            if missing:
                raise ValueError(f"qroam='non-optimum' requires parameters: {', '.join(missing)}")
            return self._toffoli_core(
                self.L, kp1, kp1_prime, ko, ko_prime, kp2, kp2_prime, kr, kr_prime,
                br, Rl, R0, Nb, self.Nk, N1, N2, eta, B, bo_override=None
            )
        else:
            raise ValueError("qroam must be 'optimum' or 'non-optimum'.")

    def total_qubits(
        self,
        *,
        Rl, R0, eta=None,
        lam: Optional[float] = None, eps_qpe: Optional[float] = None,
        br: int = 20,
        N1: int = 10,
        N2: int = 10,
        B: int = 20,
        Nb: Optional[int] = None,
        qroam: str = "optimum",
        kr: Optional[float] = None,
        I: Optional[int] = None,
    ) -> int:
        """Compute the total logical qubit count (Eq. 99).

        System-specific parameters:
            Rl: Average two-body rank R^(ℓ≠0), from OneNormCalculator output.
            R0: One-body rank R^(0), from OneNormCalculator output.
            eta: QROAM parameter (2-adic valuation of L, default 10).

        QPE parameters (needed for the phase-estimation register size):
            lam: One-norm λ from OneNormCalculator.
            eps_qpe: Target QPE precision.  Chemical accuracy is 1.6e-3 Ha
                (1 kcal/mol); with a 5-way error budget (QPE, LCU
                truncation, coefficient loading, finite basis, finite k-mesh) the QPE
                share is eps_chem/5 ≈ 3.2e-4.  The number of QPE
                iterations is I = ceil(π·λ/ε_QPE).  If neither lam/eps_qpe
                nor I is given, the QPE register is omitted from the count.
            I: Number of QPE iterations (overrides lam/eps_qpe if given).

        Precision parameters (have defaults — see toffoli_count_per_be for
        detailed guidance on when to change these):
            br: Rotation-angle register bit precision (default 20).
            N1: One-body PREPARE phase-angle bits (default 10).
            N2: Two-body PREPARE phase-angle bits (default 10).
            B: Givens rotation bit precision (default 20).
            Nb: Number of bands per k-point (defaults to self.Nb from HDF5).
            qroam: 'optimum' (default) or 'non-optimum'.
            kr: QROAM batch size for rotation oracle (required for non-optimum).

        Returns:
            Integer qubit count.
        """
        if eta is None:
            eta = self.eta

        if Nb is None:
            if self.Nb is None:
                raise ValueError("Provide Nb or ensure it is stored from HDF5.")
            Nb = self.Nb

        if qroam.lower() == "optimum":
            p = self._opt_params(br=br, Rl=Rl, R0=R0, Nb=Nb, N1=N1, N2=N2, eta=eta, B=B)
            kr_eff = p["kr"]
            return self._total_qubits_core(
                self.L, kr_eff,
                br, Rl, R0, Nb, self.Nk, N1, N2, eta, B,
                I=I, lam=lam, eps_qpe=eps_qpe,
                bo_override=p["bo"]
            )
        elif qroam.lower() == "non-optimum":
            if kr is None:
                raise ValueError("qroam='non-optimum' requires kr.")
            return self._total_qubits_core(
                self.L, kr,
                br, Rl, R0, Nb, self.Nk, N1, N2, eta, B,
                I=I, lam=lam, eps_qpe=eps_qpe,
                bo_override=None
            )
        else:
            raise ValueError("qroam must be 'optimum' or 'non-optimum'.")

