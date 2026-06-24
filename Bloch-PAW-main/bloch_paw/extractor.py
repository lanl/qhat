#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import numpy as np
from gpaw import GPAW
from gpaw.transformers import Transformer
import h5py
import json
import datetime
from functools import partial


# ---- small helper: pack upper-triangular Hermitian matrix (row-major, i<=j) ----
def pack_hermitian(H: np.ndarray) -> np.ndarray:
    """Pack an (n,n) Hermitian matrix into length n(n+1)/2 vector (i<=j)."""
    n = H.shape[0]
    try:
        cache = pack_hermitian._triu_cache  # type: ignore[attr-defined]
    except AttributeError:
        cache = pack_hermitian._triu_cache = {}  # type: ignore[attr-defined]
    ij = cache.get(n)
    if ij is None:
        ij = np.triu_indices(n)
        cache[n] = ij
    i, j = ij
    return H[i, j]


class PawExtractor:
    """Extract PAW integrals from a GPAW calculator for quantum resource estimation.

    This version applies Bloch/k-point normalization consistently for *all* modes
    (PW, FD, LCAO) whenever Nk > 0:
      - D^a /= sqrt(Nk)
      - rho_tilde /= sqrt(Nk)
      - h /= Nk
      - kappa_soft, kappa_paw /= Nk

    The compensation-charge build (Z̃) is made consistent with version 1 by:
      - forcing real LFC coefficients (QL.real) when writing to W_aL[a]
      - ensuring W_aL[a] has a complex dtype if needed (dtype guard)
      - zeroing atom entry after finishing the atom to avoid cross-atom contamination

    Energy units:
      - GPAW internal operators are typically in Hartree, while some public APIs
        (notably get_eigenvalues) often return eV.
      - This class includes a robust unit-detection helper that determines whether
        the "raw" Hamiltonian / energy-like quantities it constructs are in Ha or eV,
        and converts consistently to the requested output unit (default eV).
    """

    def __init__(
        self,
        calc,
        nbands: Optional[int] = None,
        thr_rho: float = 0.0,
        thr_D: float = 0.0,
        thr_C: float = 0.0,
        thr_h: float = 0.0,
        thr_kappa: float = 0.0,
    ):
        self.calc = calc
        self.wfs = calc.wfs
        self.default_nbands = nbands
        self._thr_rho = float(thr_rho)
        self._thr_D = float(thr_D)
        self._thr_C = float(thr_C)
        self._thr_h = float(thr_h)
        self._thr_kappa = float(thr_kappa)

        self._kd = getattr(self.wfs, "kd", None)
        self._pd_ghat = getattr(getattr(self.calc.density, "ghat", None), "pd", None)
        self._lfc = getattr(self.calc.density, "ghat", None)
        self._gd_work = self.wfs.gd
        self._gd_lfc = getattr(self._pd_ghat, "gd", None)
        self._T_Z2W = (
            None
            if (self._gd_lfc is None or tuple(self._gd_lfc.N_c) == tuple(self._gd_work.N_c))
            else Transformer(self._gd_lfc, self._gd_work)
        )
        self._B2pi = 2.0 * np.pi * self.calc.atoms.cell.reciprocal()
        self._Bno2pi = np.linalg.inv(self.calc.atoms.cell.array).T
        self._P_cache: Dict[tuple, Dict[int, np.ndarray]] = {}

        # caches for repeated small index maps / temporary work buffers
        self._triu_idx_cache: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
        self._W_aL_cache: Dict[Tuple[int, int], Dict[int, np.ndarray]] = {}
        self._buf_xG_cache: Dict[int, np.ndarray] = {}

        # cached energy-unit detection: scale from "raw energy-like quantity" -> eV
        self._energy_scale_to_ev: Optional[float] = None

    # -------------------------------------------------------------------------
    # Mode detection helpers (kept, but normalization no longer depends on mode)
    # -------------------------------------------------------------------------
    def _is_pw_mode(self) -> bool:
        """Return True if calculation is in PW mode."""
        mode = getattr(self.calc, "mode", None)
        name = getattr(getattr(mode, "__class__", type(None)), "__name__", "").lower()
        if "pw" in name:
            return True
        try:
            kpt0 = self.wfs.kpt_u[0]
            return hasattr(kpt0, "psit_nG") and (kpt0.psit_nG is not None)
        except (AttributeError, IndexError):
            return False

    def _projections_array(self, kpt):
        """
        Return PAW projector overlaps as a raw 2D numpy array (no casting/transpose).
        Works in PW, FD, and LCAO.
        """
        P = getattr(kpt, "projections", None)
        if P is None:
            return None

        arr = getattr(P, "array", None)
        if arr is None:
            if hasattr(P, "toarray"):
                arr = P.toarray()
            elif hasattr(P, "asarray"):
                arr = P.asarray()
            else:
                arr = np.array(P)

        A = np.asarray(arr)  # IMPORTANT: no dtype / no copy=False (NumPy 2.0 safe)
        if A.ndim != 2:
            raise ValueError(f"Unexpected projections array shape {A.shape}; expected 2-D.")
        return A

    # -------------------------------------------------------------------------
    # Guards / utilities
    # -------------------------------------------------------------------------
    def _ensure_initialized(self) -> None:
        """Ensure GPAW has built k-points/wavefunctions."""
        if getattr(self, "wfs", None) is None or getattr(self.wfs, "kpt_u", None) is None:
            atoms = getattr(self.calc, "atoms", None)
            if atoms is not None:
                try:
                    self.calc.initialize(atoms)
                    self.wfs = self.calc.wfs
                except (AttributeError, RuntimeError):
                    pass
        if self.wfs is None or getattr(self.wfs, "kpt_u", None) is None:
            raise RuntimeError("GPAW not initialized; call atoms.get_potential_energy() first.")

    @staticmethod
    def _threshold_inplace(arr: np.ndarray, thr: float) -> np.ndarray:
        if thr is None or thr <= 0.0:
            return arr
        if not arr.flags.writeable:
            arr = arr.copy()
        m = np.abs(arr) < thr
        if m.any():
            arr[m] = 0.0
        return arr

    def _triu_indices_cached(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        ij = self._triu_idx_cache.get(int(n))
        if ij is None:
            ij = np.triu_indices(int(n))
            self._triu_idx_cache[int(n)] = ij
        return ij

    # -------------------------------------------------------------------------
    # Energy scale helpers (raw -> eV), inspired by the reference snippet
    # -------------------------------------------------------------------------
    @staticmethod
    def _hartree_ev_value() -> float:
        """Return 1 Hartree in eV."""
        try:
            from ase.units import Hartree as _Ha
            return float(_Ha)
        except Exception:
            return 27.211386245988

    def _detect_energy_scale_to_ev(self, spin: Optional[int] = None) -> float:
        """
        Detect whether the 'raw' energy-like quantities we construct (e.g. via
        apply_pseudo_hamiltonian matrix elements) are in Hartree or already in eV.

        Returns:
            scale_to_ev: multiply raw energies by this factor to get eV.
                         Typically either 1.0 or Hartree_in_eV.
        """
        if self._energy_scale_to_ev is not None:
            return float(self._energy_scale_to_ev)

        self._ensure_initialized()
        Ha_ev = self._hartree_ev_value()

        iks = self._kblocks_for_spin(spin)
        if not iks:
            self._energy_scale_to_ev = 1.0
            return 1.0

        ik0 = iks[0]
        kpt = self.wfs.kpt_u[ik0]

        # Reference energies from kpt or get_eigenvalues (usually eV in GPAW/ASE interfaces)
        eps_ref = getattr(kpt, "eps_n", None)
        if eps_ref is None:
            try:
                eps_ref = self.calc.get_eigenvalues(kpt=ik0)
            except Exception:
                eps_ref = None

        # If we can't get a reference, fall back to a heuristic:
        # if apply_pseudo_hamiltonian exists, assume raw is Hartree; else assume eV.
        if eps_ref is None:
            have_apply = hasattr(self.wfs, "apply_pseudo_hamiltonian")
            self._energy_scale_to_ev = Ha_ev if have_apply else 1.0
            return float(self._energy_scale_to_ev)

        eps_ref = np.asarray(eps_ref, dtype=float)
        eps_ref = eps_ref[np.isfinite(eps_ref)]
        if eps_ref.size == 0:
            have_apply = hasattr(self.wfs, "apply_pseudo_hamiltonian")
            self._energy_scale_to_ev = Ha_ev if have_apply else 1.0
            return float(self._energy_scale_to_ev)

        have_apply = hasattr(self.wfs, "apply_pseudo_hamiltonian")
        psit = getattr(kpt, "psit", None)

        # If we can't form operator matrix elements, we assume the only available energies are eV.
        if not have_apply or psit is None:
            self._energy_scale_to_ev = 1.0
            return 1.0

        try:
            ham = self.calc.hamiltonian
            Ht = partial(self.wfs.apply_pseudo_hamiltonian, kpt, ham)

            try:
                tmp = psit.new()
            except Exception:
                tmp = None

            Hloc = psit.matrix_elements(operator=Ht, result=tmp, out=None, symmetric=True, cc=True)
            Harr = Hloc.array if hasattr(Hloc, "array") else np.asarray(Hloc)

            ntest = int(min(8, Harr.shape[0], Harr.shape[1], eps_ref.size))
            if ntest <= 0:
                self._energy_scale_to_ev = Ha_ev
                return float(self._energy_scale_to_ev)

            diag_raw = np.diag(np.asarray(Harr[:ntest, :ntest], dtype=np.complex128)).real
            diag_raw = diag_raw[np.isfinite(diag_raw)]
            eps_t = eps_ref[:diag_raw.size]
            eps_t = eps_t[np.isfinite(eps_t)]

            if diag_raw.size == 0 or eps_t.size == 0:
                self._energy_scale_to_ev = Ha_ev
                return float(self._energy_scale_to_ev)

            a = float(np.median(np.abs(diag_raw)))
            b = float(np.median(np.abs(eps_t)))
            if a <= 0.0 or b <= 0.0:
                # Heuristic fallback: small raw magnitudes often imply Hartree scale
                self._energy_scale_to_ev = Ha_ev if float(np.median(np.abs(diag_raw))) < 10.0 else 1.0
                return float(self._energy_scale_to_ev)

            ratio = b / a
            self._energy_scale_to_ev = Ha_ev if abs(ratio - Ha_ev) < abs(ratio - 1.0) else 1.0
            return float(self._energy_scale_to_ev)

        except Exception:
            # Conservative: matrix-elements path likely in Hartree
            self._energy_scale_to_ev = Ha_ev
            return float(self._energy_scale_to_ev)

    def _convert_energy_matrix(self, M: np.ndarray, *, from_scale_to_ev: float, to_unit: str) -> np.ndarray:
        """
        Convert an energy-like matrix M from 'raw' units to desired unit.

        Args:
            from_scale_to_ev: multiplier to convert raw -> eV.
                              (1.0 if raw already eV, Hartree_in_eV if raw is Hartree)
            to_unit: 'eV' or 'Ha'
        """
        to = to_unit.lower()
        Ha_ev = self._hartree_ev_value()
        if to == "ev":
            return M * float(from_scale_to_ev)
        if to == "ha":
            return M * (float(from_scale_to_ev) / Ha_ev)
        raise ValueError("to_unit must be 'eV' or 'Ha'")

    # -------------------------------------------------------------------------
    # K-mesh helpers
    # -------------------------------------------------------------------------
    def kmesh_vectors(self, where: str = "ibz", units: str = "cart") -> np.ndarray:
        """Return k-point vectors from the GPAW k-point descriptor."""
        self._ensure_initialized()
        kd = getattr(self.wfs, "kd", None)
        if kd is None:
            raise RuntimeError("No k-point descriptor (kd) on wavefunctions.")
        wl = where.lower()
        if wl == "ibz":
            k_red = kd.ibzk_kc
        elif wl == "bz":
            k_red = kd.bzk_kc
        else:
            raise ValueError("where must be 'ibz' or 'bz'")

        u = units.lower()
        if u == "reduced":
            return np.array(k_red, dtype=np.float64)
        if u == "cart":
            return (k_red @ self._B2pi).astype(np.float64, copy=False)
        if u == "cart_no2pi":
            return (k_red @ self._Bno2pi).astype(np.float64, copy=False)
        raise ValueError("units must be 'reduced' or 'cart' or 'cart_no2pi'")

    def kpoint_weights(self, where: str = "ibz"):
        self._ensure_initialized()
        kd = getattr(self.wfs, "kd", None)
        wl = where.lower()
        if wl == "ibz":
            return np.array(kd.ibzk_kc, dtype=np.float64), np.array(kd.weight_k, dtype=np.float64)
        if wl == "bz":
            kpts = np.array(kd.bzk_kc, dtype=np.float64)
            Nk = kpts.shape[0]
            return kpts, np.full(Nk, 1.0 / Nk, dtype=np.float64)
        raise ValueError("where must be 'ibz' or 'bz'")

    def plane_waves_count(self) -> int:
        """Return the average number of plane waves per k-point, ⟨N_pw(k)⟩ = (1/Nk) Σ_k N_pw(k)."""
        pd = getattr(self.wfs, "pd", None)
        if pd is not None and hasattr(pd, "ng_q"):
            try:
                ng_q = np.asarray(pd.ng_q, dtype=np.int64)
                # ng_q is the number of plane waves for each k-point (q-point) on this rank.
                # Return the average over k-points.
                if ng_q.size == 0:
                    return 0
                return int(np.rint(float(np.sum(ng_q)) / float(ng_q.size)))
            except Exception:
                pass

        # Fallback: estimate PW basis size from the real-space grid spacing h by
        # converting to an effective plane-wave cutoff and counting reciprocal vectors.
        #
        # For an FFT grid with spacing h (in Bohr), a consistent PW cutoff is
        #   E_cut (Ha) = (π / h)^2 / 4,
        # which implies |G| <= G_cut = sqrt(2 E_cut).
        #
        # NOTE: this counts reciprocal lattice vectors inside a sphere (PW basis),
        # not all Fourier modes on the FFT lattice.
        gd = getattr(self.wfs, "gd", getattr(self.calc.density, "gd", None))
        if gd is None or not hasattr(gd, "N_c"):
            return 0

        try:
            from ase.units import Bohr as _Bohr
            bohr_in_A = float(_Bohr)  # Bohr in Å
        except Exception:
            bohr_in_A = 0.529177210903

        Nx, Ny, Nz = map(int, gd.N_c)
        if Nx <= 0 or Ny <= 0 or Nz <= 0:
            return 0

        # Estimate per-axis real-space spacing from the cell vectors.
        # (Works best for near-orthorhombic cells; for skewed cells this is an approximation.)
        A = np.asarray(self.calc.atoms.cell.array, dtype=np.float64)  # rows are cell vectors in Å
        hA = np.array(
            [np.linalg.norm(A[0]) / Nx, np.linalg.norm(A[1]) / Ny, np.linalg.norm(A[2]) / Nz],
            dtype=np.float64,
        )
        h_bohr = float(np.min(hA) / bohr_in_A)  # choose the finest spacing (conservative PW cutoff)

        if h_bohr <= 0.0:
            return 0

        Ecut_Ha = (np.pi / h_bohr) ** 2 / 4.0
        Gcut = np.sqrt(2.0 * Ecut_Ha)  # in Bohr^-1

        # Reciprocal basis (2π convention) in Å^-1 -> convert to Bohr^-1
        B2pi_Ainv = self._B2pi  # rows are reciprocal vectors in Å^-1
        B2pi_Bohrinv = np.asarray(B2pi_Ainv, dtype=np.float64) * bohr_in_A

        # Build a bounding box in integer coefficient space, then count vectors in the sphere.
        bnorm = np.linalg.norm(B2pi_Bohrinv, axis=1)
        bmin = float(np.min(bnorm[bnorm > 0.0])) if np.any(bnorm > 0.0) else 0.0
        if bmin <= 0.0:
            return 0

        nmax = int(np.ceil(Gcut / bmin)) + 1
        n = np.arange(-nmax, nmax + 1, dtype=np.int32)

        n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
        nn = np.stack([n1, n2, n3], axis=-1).reshape(-1, 3)  # (N,3)

        G_cart = nn @ B2pi_Bohrinv  # (N,3) in Bohr^-1
        G2 = np.einsum("ij,ij->i", G_cart, G_cart, optimize=True)

        # PW basis: |G|^2/2 < Ecut, excluding G=0
        mask = (G2 > 0.0) & (0.5 * G2 < Ecut_Ha)
        return int(np.count_nonzero(mask))

    # -------------------------------------------------------------------------
    # k-blocks (spin) and band counting (PW & LCAO safe)
    # -------------------------------------------------------------------------
    def _kblocks_for_spin(self, spin: Optional[int]) -> List[int]:
        self._ensure_initialized()
        if spin is None:
            spin = self.wfs.kpt_u[0].s
        return [ik for ik, kpt in enumerate(self.wfs.kpt_u) if kpt.s == spin]

    def _nbands_at_k(self, ik: int) -> int:
        kpt = self.wfs.kpt_u[ik]
        for attr in ("eps_n", "f_n", "psit_nG", "C_nM"):
            arr = getattr(kpt, attr, None)
            if arr is not None:
                return len(arr)
        raise RuntimeError(f"Cannot determine nbands at k={ik}")

    def _choose_nb(self, iks: List[int], nbands: Optional[int]) -> int:
        self._ensure_initialized()
        nb_all = min(self._nbands_at_k(ik) for ik in iks)
        N = nbands if nbands is not None else (self.default_nbands or nb_all)
        if N > nb_all:
            raise ValueError(f"Requested nbands={N} > available min per k={nb_all}.")
        return int(N)

    def _get_q_index_map(self, iks: List[int]) -> Dict[int, int]:
        pd = self._pd_ghat
        if pd is None:
            return {ik: 0 for ik in iks}
        if hasattr(pd, "ng_q"):
            nq_pd = len(pd.ng_q)
        elif hasattr(pd, "myng_q"):
            nq_pd = len(pd.myng_q)
        else:
            nq_pd = 1
        return {ik: (getattr(self.wfs.kpt_u[ik], "q", 0) % max(nq_pd, 1)) for ik in iks}

    # -------------------------------------------------------------------------
    # PAW overlaps from GRID wavefunctions (PW & LCAO)
    # -------------------------------------------------------------------------
    def _P_overlaps_all(self, iks: List[int], Nb: int, use_cache: bool = True) -> Dict[int, np.ndarray]:
        """Build per-atom projector overlaps P_by_atom[a][ik_idx, n, i] from kpt.projections."""
        key = (tuple(iks), int(Nb))
        if use_cache and key in self._P_cache:
            return self._P_cache[key]

        ni_a = [getattr(setup, "ni", None) for setup in self.wfs.setups]
        if any(n is None for n in ni_a):
            raise AttributeError("Could not determine projector channel counts (setup.ni).")
        offsets = np.cumsum([0] + ni_a[:-1])
        nI_total = int(sum(ni_a))

        P_by_atom: Dict[int, np.ndarray] = {}
        for ik_idx, ik in enumerate(iks):
            kpt = self.wfs.kpt_u[ik]
            P_full = self._projections_array(kpt)
            if P_full is None:
                raise AttributeError("kpt.projections not available; cannot form PAW overlaps.")

            # Decide orientation by matching known projector width nI_total
            if P_full.shape[1] == nI_total:
                pass  # (Nb_full, nI_total)
            elif P_full.shape[0] == nI_total:
                P_full = P_full.T
            else:
                raise ValueError(f"Projections shape {P_full.shape} incompatible with total projectors {nI_total}.")

            P_full = np.asarray(P_full, dtype=np.complex128)
            Nb_full = P_full.shape[0]
            if Nb > Nb_full:
                raise ValueError(f"Requested nbands={Nb} but projections only have {Nb_full} bands.")

            P_nb = P_full[:Nb, :]

            for a, (off, ni) in enumerate(zip(offsets, ni_a)):
                sl = slice(off, off + ni)
                Pi = P_nb[:, sl]
                if a not in P_by_atom:
                    P_by_atom[a] = np.empty((len(iks), Nb, ni), dtype=np.complex128)
                P_by_atom[a][ik_idx, :, :] = Pi

        self._P_cache[key] = P_by_atom
        return P_by_atom

    # -------------------------------------------------------------------------
    # D and rho_tilde (normalized for all modes when Nk>0)
    # -------------------------------------------------------------------------
    def D_matrix(
        self,
        nbands: Optional[int] = None,
        spin: Optional[int] = None,
        thr: Optional[float] = None,
        use_cache: bool = True,
    ) -> Dict[int, np.ndarray]:
        """Compute the projector density matrix D^a for each atom, with D /= sqrt(Nk) if Nk>0."""
        iks = self._kblocks_for_spin(spin)
        Nk = len(iks)
        N = self._choose_nb(iks, nbands)
        P_all = self._P_overlaps_all(iks, N, use_cache=use_cache)
        tol = self._thr_D if (thr is None) else float(thr)

        out: Dict[int, np.ndarray] = {}
        for a, P in P_all.items():
            D = np.einsum("kpi,lqj->kplqij", P.conj(), P, optimize=True)
            if Nk > 0:
                D = D / np.sqrt(float(Nk))
            out[a] = self._threshold_inplace(D, tol)
        return out

    def _build_psi_stack(self, iks: List[int], N: int, spin: Optional[int], broadcast: bool) -> np.ndarray:
        psi_blocks = []
        for ik in iks:
            kpt = self.wfs.kpt_u[ik]
            q, s = getattr(kpt, "q", 0), kpt.s
            psi_k = [self.calc.get_pseudo_wave_function(n, q, s, broadcast) for n in range(N)]
            psi_blocks.append(np.stack([np.asarray(ps, dtype=np.complex128) for ps in psi_k], axis=0))
        psi = np.stack(psi_blocks, axis=0)
        return psi

    def _rho_rs_from_psi(self, psi: np.ndarray, r_lin: np.ndarray, s_lin: np.ndarray, Nb: int) -> np.ndarray:
        Nk = psi.shape[0]
        M = Nk * Nb
        r_lin = np.asarray(r_lin, dtype=np.int64)
        s_lin = np.asarray(s_lin, dtype=np.int64)
        if (
            (r_lin.min(initial=0) < 0)
            or (s_lin.min(initial=0) < 0)
            or (r_lin.max(initial=-1) >= M)
            or (s_lin.max(initial=-1) >= M)
        ):
            raise ValueError("r_lin/s_lin out of bounds for (Nk*Nb).")
        kr = r_lin // Nb
        ur = r_lin % Nb
        ks = s_lin // Nb
        vs = s_lin % Nb
        return psi[kr, ur].conj() * psi[ks, vs]

    # ------------------------ UPDATED HERE ------------------------
    def _add_compensation_charges(
        self,
        iks: List[int],
        Nk: int,
        N: int,
        P_all: Dict[int, np.ndarray],
        pd_ghat,
        lfc,
        rho_shape: tuple,
        early_Z_threshold: float,
    ) -> np.ndarray:
        """Compute the compensation charge contribution Z̃ to the pair density.

        Changes vs original v1 to match v2:
          1) Force real LFC coefficients: W_aL[a] <- QL[p,q].real
          2) Dtype-guard W_aL[a] (ensure complex128 buffer if needed)
          3) Zero W_aL[a] after finishing each atom to avoid cross-atom contamination
        """
        ik_to_qref = self._get_q_index_map(iks)
        Z_work = np.zeros(rho_shape, dtype=np.complex128)

        for a, P in P_all.items():
            setup = self.wfs.setups[a]
            Delta_pL = setup.Delta_pL
            nL = Delta_pL.shape[1]

            ni = int(P.shape[-1])
            ti, tj = self._triu_indices_cached(ni)
            QL = np.empty((N, N, nL), dtype=np.complex128)

            # We will track which qref dicts we touched for this atom, and clean them at end.
            touched_WaL: List[Tuple[int, Dict[int, np.ndarray]]] = []

            for ikp_idx in range(Nk):
                qref = ik_to_qref[iks[ikp_idx]]
                if qref not in self._buf_xG_cache:
                    self._buf_xG_cache[qref] = pd_ghat.zeros(1, q=qref)
                buf_xG = self._buf_xG_cache[qref]

                w_key = (int(qref), int(nL))
                W_aL = self._W_aL_cache.get(w_key)
                if W_aL is None:
                    W_aL = lfc.dict()
                    self._W_aL_cache[w_key] = W_aL

                # ---- dtype guard (like v2 minimal patch) ----
                Wa = W_aL.get(a, None)
                if Wa is not None and Wa.dtype != np.complex128:
                    W_aL[a] = Wa.astype(np.complex128, copy=False)
                elif Wa is None:
                    # Some GPAW builds may omit atom key until first use; ensure it exists.
                    # lfc.dict() typically contains all atoms, but be safe.
                    Wtmp = lfc.dict()
                    if a in Wtmp:
                        W_aL[a] = np.asarray(Wtmp[a]).astype(np.complex128, copy=False)

                touched_WaL.append((qref, W_aL))

                for ikq_idx in range(ikp_idx, Nk):
                    Pp = P[ikp_idx]
                    Pq = P[ikq_idx]

                    Dpq = np.einsum("pi,qj->pqij", Pp.conj(), Pq, optimize=True)
                    Hpq = 0.5 * (Dpq + Dpq.swapaxes(-1, -2).conj())
                    di = np.arange(Hpq.shape[-1])
                    Hpq[..., di, di] = Hpq[..., di, di].real

                    dp_all = Hpq[..., ti, tj]
                    QL[...] = dp_all @ Delta_pL

                    for p in range(N):
                        for q in range(N):
                            if early_Z_threshold and early_Z_threshold > 0.0:
                                if np.max(np.abs(QL[p, q])) < early_Z_threshold:
                                    continue

                            # Only touch this atom's entry (do not clear entire dict).
                            W_aL[a][...] = 0.0
                            # ---- enforce real coefficients (match v2) ----
                            W_aL[a][...] = np.asarray(QL[p, q].real, dtype=W_aL[a].dtype)

                            buf_xG.fill(0.0)
                            lfc.add(buf_xG, W_aL, q=qref)
                            z_r = pd_ghat.ifft(buf_xG[0], q=qref)

                            if self._T_Z2W is not None:
                                z_w = self._gd_work.empty(dtype=z_r.dtype)
                                self._T_Z2W.apply(z_r, z_w)
                                z_r = z_w

                            Z_work[ikp_idx, p, ikq_idx, q] += z_r.astype(np.complex128, copy=False)
                            if ikq_idx != ikp_idx:
                                Z_work[ikq_idx, q, ikp_idx, p] += z_r.conj().astype(np.complex128, copy=False)

            # ---- end-of-atom cleanup to avoid contamination ----
            for _qref, W_aL in touched_WaL:
                if a in W_aL:
                    W_aL[a][...] = 0.0

        return Z_work
    # ------------------------------------------------------------

    def pair_density_matrix(
        self,
        nbands: Optional[int] = None,
        spin: Optional[int] = None,
        broadcast: bool = True,
        thr: Optional[float] = None,
        include_compensation: bool = True,
        use_cache: bool = True,
        early_Z_threshold: float = 0.0,
    ) -> np.ndarray:
        """Compute rho_tilde with rho_tilde /= sqrt(Nk) if Nk>0 (all modes)."""
        iks = self._kblocks_for_spin(spin)
        Nk = len(iks)
        N = self._choose_nb(iks, nbands)

        psi = self._build_psi_stack(iks, N, spin, broadcast)
        rho_val = psi.conj()[:, :, None, None, ...] * psi[None, None, :, :, ...]
        del psi

        if not include_compensation:
            rho = rho_val
        else:
            lfc = self._lfc
            pd_ghat = self._pd_ghat
            if pd_ghat is None or lfc is None:
                rho = rho_val
            else:
                P_all = self._P_overlaps_all(iks, N, use_cache=use_cache)
                Z_work = self._add_compensation_charges(
                    iks, Nk, N, P_all, pd_ghat, lfc, rho_val.shape, early_Z_threshold
                )
                rho = rho_val + Z_work

        if Nk > 0:
            rho = rho / float(np.sqrt(Nk))

        return self._threshold_inplace(rho, self._thr_rho if (thr is None) else float(thr))

    # -------------------------------------------------------------------------
    # C tensor
    # -------------------------------------------------------------------------
    def C_tensor(
        self,
        a: Optional[int] = None,
        unpack: bool = True,
        ensure_fresh: bool = False,
        thr: Optional[float] = None,
    ):
        def _get_Mpp(setup):
            if ensure_fresh and hasattr(setup, "calculate_coulomb_corrections"):
                try:
                    _Mp, Mpp = setup.calculate_coulomb_corrections(None, None, None, None, None)
                    return np.asarray(Mpp, dtype=np.float64)
                except Exception:
                    pass
            if not hasattr(setup, "M_pp"):
                raise AttributeError("Setup has no attribute 'M_pp'")
            return np.asarray(setup.M_pp, dtype=np.float64)

        def _unpack(setup, Mpp):
            ni = getattr(setup, "ni", None)
            if ni is None:
                l_j = getattr(setup, "l_j", None)
                ni = len(l_j) if l_j is not None else None
            if ni is None:
                raise AttributeError("Could not determine number of projector channels (ni).")
            pairs = [(i1, i2) for i1 in range(ni) for i2 in range(i1, ni)]
            if Mpp.shape != (len(pairs), len(pairs)):
                raise ValueError(f"Unexpected M_pp shape {Mpp.shape} for ni={ni}")
            C = np.zeros((ni, ni, ni, ni), dtype=Mpp.dtype)
            for p_idx, (i1, i2) in enumerate(pairs):
                for q_idx, (i3, i4) in enumerate(pairs):
                    v = Mpp[p_idx, q_idx]
                    C[i1, i2, i3, i4] = v
                    C[i2, i1, i3, i4] = v
                    C[i1, i2, i4, i3] = v
                    C[i2, i1, i4, i3] = v
            return C

        tol = self._thr_C if (thr is None) else float(thr)

        if isinstance(a, int):
            setup = self.wfs.setups[a]
            if unpack:
                C = _unpack(setup, _get_Mpp(setup))
                return self._threshold_inplace(C, tol)
            else:
                Mpp = _get_Mpp(setup)
                return self._threshold_inplace(Mpp, tol)

        out: Dict[int, np.ndarray] = {}
        for aa, setup in enumerate(self.wfs.setups):
            if unpack:
                C = _unpack(setup, _get_Mpp(setup))
                out[aa] = self._threshold_inplace(C, tol)
            else:
                Mpp = _get_Mpp(setup)
                out[aa] = self._threshold_inplace(Mpp, tol)
        return out

    # -------------------------------------------------------------------------
    # One-body k-diagonal H(k)
    # -------------------------------------------------------------------------
    def one_body_Hk_list(
        self,
        nbands: Optional[int] = None,
        spin: Optional[int] = None,
        *,
        unit: str = "eV",
        orthonormalize: bool = False,
    ):
        ham = self.calc.hamiltonian
        iks = self._kblocks_for_spin(spin)
        Nb = self._choose_nb(iks, nbands)

        unit_l = unit.lower()
        if unit_l not in ("ev", "ha"):
            raise ValueError("unit must be 'eV' or 'Ha'")

        # Detect raw->eV scaling (1.0 if raw already eV; Ha_ev if raw is Hartree)
        scale_to_ev = self._detect_energy_scale_to_ev(spin=spin)

        Hk_list = []
        have_apply = hasattr(self.wfs, "apply_pseudo_hamiltonian")

        for ik in iks:
            kpt = self.wfs.kpt_u[ik]

            if orthonormalize and not getattr(self.wfs, "orthonormalized", True):
                try:
                    self.wfs.orthonormalize(kpt)
                except Exception:
                    pass

            if have_apply:
                try:
                    Ht = partial(self.wfs.apply_pseudo_hamiltonian, kpt, ham)

                    psit = kpt.psit
                    try:
                        tmp = psit.new()
                    except Exception:
                        tmp = None

                    Hloc = psit.matrix_elements(operator=Ht, result=tmp, out=None, symmetric=True, cc=True)
                    Harr = Hloc.array if hasattr(Hloc, "array") else np.asarray(Hloc)
                    Hkk_raw = np.array(Harr[:Nb, :Nb], dtype=np.complex128, copy=True)

                    # PAW augmentation: <ψ| dH |ψ>
                    try:
                        P_raw = self._projections_array(kpt)
                        P_raw = np.asarray(P_raw, dtype=np.complex128)
                        P_nb = P_raw[:Nb, :]

                        P2 = kpt.projections.new()
                        ham.dH(kpt.projections, out=P2)
                        P2_arr = getattr(P2, "array", None)
                        if P2_arr is None:
                            if hasattr(P2, "toarray"):
                                P2_arr = P2.toarray()
                            elif hasattr(P2, "asarray"):
                                P2_arr = P2.asarray()
                            else:
                                P2_arr = np.array(P2)
                        P2_arr = np.asarray(P2_arr, dtype=np.complex128)

                        if P2_arr.shape != P_raw.shape:
                            if P2_arr.T.shape == P_raw.shape:
                                P2_arr = P2_arr.T
                            else:
                                raise ValueError(
                                    f"dH·P shape {P2_arr.shape} incompatible with projections {P_raw.shape}"
                                )

                        P2_nb = P2_arr[:Nb, :]
                        Hkk_raw += P_nb.conj() @ P2_nb.T
                    except Exception:
                        pass

                    try:
                        ham.xc.correct_hamiltonian_matrix(kpt, Hkk_raw)
                    except Exception:
                        pass

                    Hkk_raw = 0.5 * (Hkk_raw + Hkk_raw.conj().T)

                    # Convert from raw to requested unit consistently
                    Hkk = self._convert_energy_matrix(Hkk_raw, from_scale_to_ev=scale_to_ev, to_unit=unit)
                    Hk_list.append(Hkk)
                    continue
                except Exception:
                    pass

            # Fallback: get_eigenvalues is treated as eV (ASE/GPAW user-facing convention)
            evals_ev = np.asarray(self.calc.get_eigenvalues(kpt=ik), dtype=float)[:Nb]
            if unit_l == "ev":
                evals = evals_ev
            else:
                evals = evals_ev / self._hartree_ev_value()
            Hkk = np.diag(evals.astype(np.complex128))
            Hk_list.append(Hkk)

        return Hk_list

    def one_body_integrals_kdiag(
        self,
        nbands: Optional[int] = None,
        spin: Optional[int] = None,
        *,
        unit: str = "eV",
    ) -> np.ndarray:
        """Assemble full one-body tensor, and normalize by Nk if Nk>0 (all modes)."""
        Hk_list = self.one_body_Hk_list(nbands=nbands, spin=spin, unit=unit)
        Nk = len(Hk_list)
        Nb = Hk_list[0].shape[0]
        h = np.zeros((Nk, Nb, Nk, Nb), dtype=np.complex128)

        thr = self._thr_h
        for k, Hkk in enumerate(Hk_list):
            if thr > 0.0:
                M = np.array(Hkk, dtype=np.complex128, copy=True)
                m = np.abs(M) < thr
                if m.any():
                    M[m] = 0.0
                M = 0.5 * (M + M.conj().T)
                Hkk = M
            h[k, :, k, :] = Hkk

        return h

    # -------------------------------------------------------------------------
    # Coulomb kernel
    # -------------------------------------------------------------------------
    def _coulomb_kernel_vG(self, g0: str = "zero") -> np.ndarray:
        Nx, Ny, Nz = map(int, self._gd_work.N_c)
        fx = np.fft.fftfreq(Nx, d=1.0 / Nx)
        fy = np.fft.fftfreq(Ny, d=1.0 / Ny)
        fz = np.fft.fftfreq(Nz, d=1.0 / Nz)
        nn = np.stack(np.meshgrid(fx, fy, fz, indexing="ij"), axis=-1)
        G_cart = np.tensordot(nn, self._B2pi, axes=(3, 0))
        G2 = np.einsum("...i,...i->...", G_cart, G_cart, optimize=True)
        V = float(abs(np.linalg.det(self.calc.atoms.cell.array)))
        vG = np.zeros_like(G2, dtype=np.float64)
        mask = G2 > 0.0
        vG[mask] = 4.0 * np.pi / (V * G2[mask])
        vG[~mask] = 0.0
        return vG

    # -------------------------------------------------------------------------
    # Two-body integrals (position space) – FULL rank-8
    # -------------------------------------------------------------------------
    @staticmethod
    def _ut_index_to_rs(idx_vec: np.ndarray, M: int) -> Tuple[np.ndarray, np.ndarray]:
        idx = idx_vec.astype(np.int64)
        M2 = 2 * M
        r = np.floor((M2 + 1 - np.sqrt((M2 + 1) ** 2 - 8 * idx)) / 2).astype(np.int64)
        f_r = (r * (M2 - r + 1)) // 2
        too_big = f_r > idx
        if np.any(too_big):
            r[too_big] -= 1
            f_r = (r * (M2 - r + 1)) // 2
        off = idx - f_r
        s = r + off
        return r, s

    # -------------------------------------------------------------------------
    # Two-body integrals (position space) – REDUCED block-diagonal components only
    # -------------------------------------------------------------------------
    def _compute_soft_kappa_blockdiag(self, Nk, Nb, rho_all, psi, vG) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute ONLY the two required soft (grid) blocks:

          K_A[k,i,j,kp,l] = kappa[k,i,kp,l,k,j,kp,l]
          K_B[k,i,j,kp,l] = kappa[k,i,k,j,kp,l,kp,l]

        Returns:
            K_A_soft, K_B_soft with shape (Nk, Nb, Nb, Nk, Nb)
        """
        # Grid shape + volume element
        if rho_all is not None:
            Nx, Ny, Nz = map(int, rho_all.shape[-3:])
        else:
            Nx, Ny, Nz = map(int, psi.shape[-3:])

        Vcell = float(abs(np.linalg.det(self.calc.atoms.cell.array)))
        dV = Vcell / float(Nx * Ny * Nz)
        Ng = int(Nx * Ny * Nz)

        K_A = np.zeros((Nk, Nb, Nb, Nk, Nb), dtype=np.complex128)
        K_B = np.zeros((Nk, Nb, Nb, Nk, Nb), dtype=np.complex128)

        # --- accessors for rho blocks ---
        def rho_kikpl(k: int, i: int, kp: int, l: int) -> np.ndarray:
            if rho_all is not None:
                return rho_all[k, i, kp, l]
            # include_compensation==False path: use psi directly
            return psi[k, i].conj() * psi[kp, l]

        def rho_kikj(k: int, i: int, j: int) -> np.ndarray:
            if rho_all is not None:
                return rho_all[k, i, k, j]
            return psi[k, i].conj() * psi[k, j]

        def rho_kplkpl(kp: int, l: int) -> np.ndarray:
            if rho_all is not None:
                return rho_all[kp, l, kp, l]
            return psi[kp, l].conj() * psi[kp, l]

        # --- main loops ---
        # For each fixed q=(kp,l):
        #   Type A: Gram matrix over {rho[(k,i),q]} for each k separately
        #   Type B: contraction of rho[(k,i),(k,j)] against phi_q = v * rho[q,q]
        for kp in range(Nk):
            for l in range(Nb):
                # phi_qq = v * rho[q,q]
                rho_qq = rho_kplkpl(kp, l)
                phi_qq = np.fft.ifftn(
                    np.fft.fftn(rho_qq, axes=(-3, -2, -1)) * vG,
                    axes=(-3, -2, -1),
                )
                phi_qq_flat = phi_qq.reshape(Ng)

                # --- Type A ---
                # For each k, build R_i = rho[(k,i), (kp,l)] vectors (i=0..Nb-1)
                for k in range(Nk):
                    R = np.empty((Nb, Ng), dtype=np.complex128)
                    for i in range(Nb):
                        R[i, :] = rho_kikpl(k, i, kp, l).reshape(Ng)

                    Phi = np.fft.ifftn(
                        np.fft.fftn(R.reshape(Nb, Nx, Ny, Nz), axes=(-3, -2, -1)) * vG[None, ...],
                        axes=(-3, -2, -1),
                    ).reshape(Nb, Ng)

                    # K_A[k,i,j,kp,l] = dV * sum_g conj(R_i[g]) * Phi_j[g]
                    K_A[k, :, :, kp, l] = dV * (R.conj() @ Phi.T)

                    # --- Type B ---
                    # K_B[k,i,j,kp,l] = dV * sum_g conj(rho[(k,i),(k,j)](g)) * phi_qq(g)
                    # Build rho_k (Nb,Nb,Ng) and contract once.
                    rho_k = np.empty((Nb, Nb, Ng), dtype=np.complex128)
                    for i in range(Nb):
                        for j in range(Nb):
                            rho_k[i, j, :] = rho_kikj(k, i, j).reshape(Ng)
                    K_B[k, :, :, kp, l] = dV * np.tensordot(rho_k.conj(), phi_qq_flat, axes=([2], [0]))

        # Hermitian symmetrize the (i,j) blocks (numerically safer)
        K_A = 0.5 * (K_A + np.swapaxes(K_A.conj(), 1, 2))
        K_B = 0.5 * (K_B + np.swapaxes(K_B.conj(), 1, 2))
        return K_A, K_B

    def _compute_paw_kappa_blockdiag(self, Nk, Nb, spin) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute ONLY the PAW on-site corrections for the two required blocks:

          K_A_paw[k,i,j,kp,l] = kappa_paw[k,i,kp,l,k,j,kp,l]
          K_B_paw[k,i,j,kp,l] = kappa_paw[k,i,k,j,kp,l,kp,l]

        Returns:
            K_A_paw, K_B_paw with shape (Nk, Nb, Nb, Nk, Nb)
        """
        K_A = np.zeros((Nk, Nb, Nb, Nk, Nb), dtype=np.complex128)
        K_B = np.zeros((Nk, Nb, Nb, Nk, Nb), dtype=np.complex128)

        C_all = self.C_tensor(a=None, unpack=True)
        D_all = self.D_matrix(nbands=Nb, spin=spin)

        for a, C in C_all.items():
            D = D_all[a]  # (Nk,Nb,Nk,Nb,na,na)
            na = C.shape[0]
            Pna = na * na
            Cflat = C.reshape(Pna, Pna)

            # For each q=(kp,l), precompute V_q := Cflat @ vec(D[q,q])
            for kp in range(Nk):
                for l in range(Nb):
                    Dqq = D[kp, l, kp, l, :, :].reshape(Pna)  # vec
                    Vq = Cflat @ Dqq  # (Pna,)

                    # --- Type A (Gram under Cflat) ---
                    # For each k: D_i := vec(D[(k,i), q]) for i=0..Nb-1
                    for k in range(Nk):
                        Dkq = D[k, :, kp, l, :, :].reshape(Nb, Pna)  # (Nb,Pna)

                        # K_A[k,i,j,kp,l] += D_i^* · (Cflat @ D_j)
                        Y = (Cflat @ Dkq.T)  # (Pna,Nb)
                        K_A[k, :, :, kp, l] += Dkq.conj() @ Y  # (Nb,Nb)

                        # --- Type B ---
                        # Need D[(k,i),(k,j)] contracted with Vq
                        Dkk = D[k, :, k, :, :, :].reshape(Nb, Nb, Pna)  # (Nb,Nb,Pna)
                        K_B[k, :, :, kp, l] += np.tensordot(Dkk.conj(), Vq, axes=([2], [0]))

        # Hermitian symmetrize (i,j)
        K_A = 0.5 * (K_A + np.swapaxes(K_A.conj(), 1, 2))
        K_B = 0.5 * (K_B + np.swapaxes(K_B.conj(), 1, 2))
        return K_A, K_B

    def two_body_integrals(
        self,
        nbands: Optional[int] = None,
        *,
        spin: Optional[int] = None,
        include_compensation: bool = True,
        g0: str = "zero",
        unit: str = "eV",
        thr: float = 0.0,
        block_pairs: int = 64,
        reuse_rho: Optional[np.ndarray] = None,
        reuse_vG: Optional[np.ndarray] = None,
        return_split: bool = False,
    ):
        self._ensure_initialized()
        iks = self._kblocks_for_spin(spin)
        Nk = len(iks)
        Nb = self._choose_nb(iks, nbands)

        rho_all = reuse_rho
        psi = None
        if rho_all is None:
            if include_compensation:
                rho_all = self.pair_density_matrix(nbands=Nb, spin=spin, include_compensation=include_compensation)
            else:
                psi = self._build_psi_stack(iks, Nb, spin, broadcast=True)

        vG = reuse_vG if reuse_vG is not None else self._coulomb_kernel_vG(g0=g0)

        # ---- reduced block-diagonal soft/paw blocks ----
        K_A_soft, K_B_soft = self._compute_soft_kappa_blockdiag(Nk, Nb, rho_all, psi, vG)
        K_A_paw, K_B_paw = self._compute_paw_kappa_blockdiag(Nk, Nb, spin)

        # ---- UPDATED UNITS: detect + convert consistently ----
        unit_l = unit.lower()
        if unit_l not in ("ev", "ha"):
            raise ValueError("unit must be 'Ha' or 'eV'")

        scale_to_ev = self._detect_energy_scale_to_ev(spin=spin)  # raw -> eV (1 or Ha_ev)

        if unit_l == "ev":
            fac = float(scale_to_ev)
        else:
            Ha_ev = self._hartree_ev_value()
            fac = float(scale_to_ev) / Ha_ev

        K_A_soft *= fac
        K_B_soft *= fac
        K_A_paw *= fac
        K_B_paw *= fac

        K_A = K_A_soft + K_A_paw
        K_B = K_B_soft + K_B_paw

        thr_eff = float(thr) if (thr and thr > 0.0) else self._thr_kappa
        if thr_eff > 0.0:
            mA = np.abs(K_A) < thr_eff
            if mA.any():
                K_A = K_A.copy()
                K_A[mA] = 0.0
            mB = np.abs(K_B) < thr_eff
            if mB.any():
                K_B = K_B.copy()
                K_B[mB] = 0.0

            msA = np.abs(K_A_soft) < thr_eff
            if msA.any():
                K_A_soft = K_A_soft.copy()
                K_A_soft[msA] = 0.0
            mpA = np.abs(K_A_paw) < thr_eff
            if mpA.any():
                K_A_paw = K_A_paw.copy()
                K_A_paw[mpA] = 0.0

            msB = np.abs(K_B_soft) < thr_eff
            if msB.any():
                K_B_soft = K_B_soft.copy()
                K_B_soft[msB] = 0.0
            mpB = np.abs(K_B_paw) < thr_eff
            if mpB.any():
                K_B_paw = K_B_paw.copy()
                K_B_paw[mpB] = 0.0

        # NOTE: return shape changed by design:
        #   (K_A, K_B) instead of full rank-8 kappa
        return (K_A, K_B, K_A_soft, K_A_paw, K_B_soft, K_B_paw) if return_split else (K_A, K_B)

    # -------------------------------------------------------------------------
    # HDF5 export
    # -------------------------------------------------------------------------
    def export_hdf5(
        self,
        filepath: str = "paw_ingredients.h5",
        *,
        spin: Optional[int] = None,
        include_compensation: bool = True,
        supercell_size: Optional[Tuple[int, int, int]] = None,
        write_one_body: bool = True,
        one_body_unit: str = "eV",
        write_two_body: bool = True,
        two_body_unit: str = "eV",
        two_body_g0: str = "zero",
        two_body_thr: float = 0.0,
        write_two_body_split: bool = False,
        compression: str = "gzip",
        compression_opts: int = 4,
        chunks: bool = True,
    ) -> str:
        nbands = self.default_nbands

        mode_obj = getattr(self.calc, "mode", None)
        if isinstance(mode_obj, str):
            mode_name = mode_obj
            ecut_ev = None
        else:
            mode_name = getattr(mode_obj, "__class__", type(mode_obj)).__name__
            ecut_ev = getattr(mode_obj, "ecut", None)

        n_electrons = int(self.calc.get_number_of_electrons())
        Npw = self.plane_waves_count()
        k_ibz_red, w_ibz = self.kpoint_weights("ibz")
        k_ibz_cart = self.kmesh_vectors("ibz", "cart")
        k_bz_red = self.kmesh_vectors("bz", "reduced")
        k_bz_cart = self.kmesh_vectors("bz", "cart")
        A = np.array(self.calc.atoms.cell.array, dtype=np.float64)
        B2pi = 2.0 * np.pi * self.calc.atoms.cell.reciprocal()
        L_size = np.array(supercell_size if supercell_size is not None else (1, 1, 1), dtype=np.int64)

        rho_tilde = self.pair_density_matrix(nbands=nbands, spin=spin, include_compensation=include_compensation)
        Ca_dict = self.C_tensor(a=None, unpack=True)
        Da_dict = self.D_matrix(nbands=nbands, spin=spin)

        h_pq = None
        if write_one_body:
            h_pq = self.one_body_integrals_kdiag(nbands=nbands, spin=spin, unit=one_body_unit)

        # ---- reduced two-body outputs ----
        kappa_A = kappa_B = None
        kappa_A_soft = kappa_A_paw = None
        kappa_B_soft = kappa_B_paw = None

        if write_two_body:
            vG_kernel = self._coulomb_kernel_vG(g0=two_body_g0)
            if write_two_body_split:
                (kappa_A, kappa_B, kappa_A_soft, kappa_A_paw, kappa_B_soft, kappa_B_paw) = self.two_body_integrals(
                    nbands=nbands,
                    spin=spin,
                    include_compensation=include_compensation,
                    g0=two_body_g0,
                    unit=two_body_unit,
                    thr=0.0,
                    reuse_rho=rho_tilde,
                    reuse_vG=vG_kernel,
                    return_split=True,
                )
            else:
                kappa_A, kappa_B = self.two_body_integrals(
                    nbands=nbands,
                    spin=spin,
                    include_compensation=include_compensation,
                    g0=two_body_g0,
                    unit=two_body_unit,
                    thr=0.0,
                    reuse_rho=rho_tilde,
                    reuse_vG=vG_kernel,
                    return_split=False,
                )

            if two_body_thr and two_body_thr > 0.0:
                mA = np.abs(kappa_A) < two_body_thr
                if mA.any():
                    kappa_A = kappa_A.copy()
                    kappa_A[mA] = 0.0
                mB = np.abs(kappa_B) < two_body_thr
                if mB.any():
                    kappa_B = kappa_B.copy()
                    kappa_B[mB] = 0.0

                if write_two_body_split:
                    msA = np.abs(kappa_A_soft) < two_body_thr
                    if msA.any():
                        kappa_A_soft = kappa_A_soft.copy()
                        kappa_A_soft[msA] = 0.0
                    mpA = np.abs(kappa_A_paw) < two_body_thr
                    if mpA.any():
                        kappa_A_paw = kappa_A_paw.copy()
                        kappa_A_paw[mpA] = 0.0

                    msB = np.abs(kappa_B_soft) < two_body_thr
                    if msB.any():
                        kappa_B_soft = kappa_B_soft.copy()
                        kappa_B_soft[msB] = 0.0
                    mpB = np.abs(kappa_B_paw) < two_body_thr
                    if mpB.any():
                        kappa_B_paw = kappa_B_paw.copy()
                        kappa_B_paw[mpB] = 0.0

        meta = {
            "timestamp": datetime.datetime.now().isoformat(),
            "code": "GPAW+ASE",
            "gpaw_version": getattr(__import__("gpaw"), "__version__", "unknown"),
            "ase_version": getattr(__import__("ase"), "__version__", "unknown"),
            "numpy_version": np.__version__,
            "xc": getattr(self.calc, "xc", "unknown"),
            "mode": mode_name,
            "ecut_eV": ecut_ev,
            "nbands_requested": int(nbands) if nbands is not None else None,
            "include_compensation": bool(include_compensation),
            "thresholds": {"rho": self._thr_rho, "D": self._thr_D, "C": self._thr_C},
            "one_body": {"written": bool(write_one_body), "unit": one_body_unit, "k_diagonal": True},
            "two_body": {
                "written": bool(write_two_body),
                "unit": two_body_unit,
                "definition": "Reduced block-diagonal tensors: K_A=κ[k,i,k',l,k,j,k',l], K_B=κ[k,i,k,j,k',l,k',l]",
                "g0": two_body_g0,
                "k_diagonal": False,
                "split_written": bool(write_two_body_split),
                "two_body_threshold": float(two_body_thr),
            },
            "energy_scale_raw_to_eV": float(self._detect_energy_scale_to_ev(spin=spin)),
        }

        N_atoms = len(Ca_dict)

        write_hdf5(
            filepath,
            meta=meta,
            Npw=Npw,
            n_electrons=n_electrons,
            N_atoms=N_atoms,
            L_size=L_size,
            A=A,
            B2pi=B2pi,
            k_ibz_red=k_ibz_red,
            w_ibz=w_ibz,
            k_ibz_cart=k_ibz_cart,
            k_bz_red=k_bz_red,
            k_bz_cart=k_bz_cart,
            rho_tilde=rho_tilde,
            Ca_dict=Ca_dict,
            Da_dict=Da_dict,
            h_pq=h_pq,
            one_body_unit=one_body_unit,
            kappa_A=kappa_A,
            kappa_B=kappa_B,
            kappa_A_soft=kappa_A_soft,
            kappa_A_paw=kappa_A_paw,
            kappa_B_soft=kappa_B_soft,
            kappa_B_paw=kappa_B_paw,
            two_body_unit=two_body_unit,
            two_body_g0=two_body_g0,
            include_compensation=include_compensation,
            write_two_body_split=write_two_body_split,
            compression=compression,
            compression_opts=compression_opts,
            chunks=chunks,
        )
        return filepath


def write_hdf5(
    filepath: str,
    *,
    meta: dict,
    Npw: int,
    n_electrons: Optional[int] = None,
    N_atoms: Optional[int] = None,
    L_size: np.ndarray,
    A: np.ndarray,
    B2pi: np.ndarray,
    k_ibz_red: np.ndarray,
    w_ibz: np.ndarray,
    k_ibz_cart: np.ndarray,
    k_bz_red: np.ndarray,
    k_bz_cart: np.ndarray,
    rho_tilde: np.ndarray,
    Ca_dict: Dict[int, np.ndarray],
    Da_dict: Dict[int, np.ndarray],
    h_pq: Optional[np.ndarray] = None,
    one_body_unit: str = "eV",
    kappa_A: Optional[np.ndarray] = None,
    kappa_B: Optional[np.ndarray] = None,
    kappa_A_soft: Optional[np.ndarray] = None,
    kappa_A_paw: Optional[np.ndarray] = None,
    kappa_B_soft: Optional[np.ndarray] = None,
    kappa_B_paw: Optional[np.ndarray] = None,
    two_body_unit: str = "eV",
    two_body_g0: str = "zero",
    include_compensation: bool = True,
    write_two_body_split: bool = False,
    compression: str = "gzip",
    compression_opts: int = 4,
    chunks: bool = True,
) -> None:
    """Write PAW ingredients to an HDF5 file."""
    Nb = int(rho_tilde.shape[1])

    with h5py.File(filepath, "w") as h5:
        h5.attrs["meta_json"] = json.dumps(meta)

        h5.create_dataset("Npw", data=np.array(Npw, dtype=np.int64))
        if n_electrons is not None:
            h5.create_dataset("n_electrons", data=np.array(n_electrons, dtype=np.int64))
        if N_atoms is not None:
            h5.create_dataset("N_atoms", data=np.array(N_atoms, dtype=np.int64))
        h5.create_dataset("supercell_size", data=L_size)

        g_lat = h5.create_group("lattice")
        g_lat.create_dataset("A_direct_ang", data=A)
        g_lat.create_dataset("B_2pi_cart_invA", data=B2pi)

        g_k = h5.create_group("kmesh")
        g_k.create_dataset("ibz/reduced", data=k_ibz_red)
        g_k.create_dataset("ibz/cart_1_perA", data=k_ibz_cart)
        g_k.create_dataset("ibz/weights", data=w_ibz)
        g_k.create_dataset("bz/reduced", data=k_bz_red)
        g_k.create_dataset("bz/cart_1_perA", data=k_bz_cart)

        g_rho = h5.create_group("rho_tilde")
        Nk_rho = rho_tilde.shape[0]
        Nx, Ny, Nz = rho_tilde.shape[-3:]
        # --- SAFE chunking: ensure each chunk < 4 GiB (HDF5 limit) ---
        if chunks:
            bytes_per_val = np.dtype(rho_tilde.dtype).itemsize  # complex128 -> 16
            grid_size = int(Nx * Ny * Nz)
            target_bytes = 512 * 1024**2  # 512 MiB

            # want (1, b, 1, b, Nx, Ny, Nz)
            # bytes = bytes_per_val * (b*b) * grid_size
            max_b2 = max(1, target_bytes // max(1, bytes_per_val * grid_size))
            b = int(max(1, min(Nb, int(np.floor(np.sqrt(max_b2))))))

            rho_chunks = (1, b, 1, b, Nx, Ny, Nz)
        else:
            rho_chunks = None

        ds_rho = g_rho.create_dataset(
            "data",
            shape=rho_tilde.shape,
            dtype=rho_tilde.dtype,
            chunks=rho_chunks,
            compression=compression,
            compression_opts=compression_opts,
        )
        for ik in range(Nk_rho):
            for ikp in range(Nk_rho):
                ds_rho[ik, :, ikp, :, :, :, :] = rho_tilde[ik, :, ikp, :, :, :, :]

        g_C = h5.create_group("C_tensor")
        for a, C in Ca_dict.items():
            g_C.create_dataset(
                f"atom_{a:04d}",
                data=C,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks,
            )

        g_D = h5.create_group("D_tensor")
        for a, D in Da_dict.items():
            g_D.create_dataset(
                f"atom_{a:04d}",
                data=D,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks,
            )

        if h_pq is not None:
            g_H = h5.create_group("one_body")
            ds = g_H.create_dataset(
                "H_kikj",
                data=h_pq,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks,
            )
            ds.attrs["unit"] = one_body_unit
            ds.attrs["k_diagonal"] = True
            ds.attrs["shape"] = "(Nk, Nb, Nk, Nb)"

        if kappa_A is not None and kappa_B is not None:
            g2 = h5.create_group("two_body")

            # Shapes: (Nk, Nb, Nb, Nk, Nb)
            dsA = g2.create_dataset(
                "kappa_A_kikplkjkpl",
                data=kappa_A,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks,
            )
            dsA.attrs["unit"] = two_body_unit
            dsA.attrs["shape"] = "(Nk, Nb, Nb, Nk, Nb)"
            dsA.attrs["definition"] = "K_A[k,i,j,kp,l] = κ[k,i,kp,l,k,j,kp,l]"
            dsA.attrs["g0"] = two_body_g0
            dsA.attrs["include_compensation"] = bool(include_compensation)

            dsB = g2.create_dataset(
                "kappa_B_kikjkplkpl",
                data=kappa_B,
                compression=compression,
                compression_opts=compression_opts,
                chunks=chunks,
            )
            dsB.attrs["unit"] = two_body_unit
            dsB.attrs["shape"] = "(Nk, Nb, Nb, Nk, Nb)"
            dsB.attrs["definition"] = "K_B[k,i,j,kp,l] = κ[k,i,k,j,kp,l,kp,l]"
            dsB.attrs["g0"] = two_body_g0
            dsB.attrs["include_compensation"] = bool(include_compensation)

            if write_two_body_split:
                dsa_s = g2.create_dataset(
                    "kappa_A_soft",
                    data=kappa_A_soft,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                dsa_s.attrs["part"] = "soft"
                dsa_s.attrs["unit"] = two_body_unit

                dsa_p = g2.create_dataset(
                    "kappa_A_paw",
                    data=kappa_A_paw,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                dsa_p.attrs["part"] = "paw"
                dsa_p.attrs["unit"] = two_body_unit

                dsb_s = g2.create_dataset(
                    "kappa_B_soft",
                    data=kappa_B_soft,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                dsb_s.attrs["part"] = "soft"
                dsb_s.attrs["unit"] = two_body_unit

                dsb_p = g2.create_dataset(
                    "kappa_B_paw",
                    data=kappa_B_paw,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                dsb_p.attrs["part"] = "paw"
                dsb_p.attrs["unit"] = two_body_unit