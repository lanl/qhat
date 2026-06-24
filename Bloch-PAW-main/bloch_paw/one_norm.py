#!/usr/bin/env python
# coding: utf-8
"""
one_norm.py

Compute the Hamiltonian one-norm λ for the Bloch–UPAW LCU decomposition.

This module provides `OneNormCalculator`, which evaluates the one-norm defined
in Eq. (37) of the referenced paper. Conceptually, λ splits into:

  (1) One-body contribution:
        Σ_k Σ_i |ε_i(k)|
      where ε_i(k) are eigenvalues of the effective one-body matrix h'(k).

  (2) “Soft” (plane-wave / pseudo-density) two-body contribution:
        (1/4) Σ_J Σ_Q Σ_G ξ_G^(J)(Q)

  (3) PAW augmentation two-body contribution:
        (1/4) Σ_J Σ_Q Σ_a Σ_{r≤s} |ε^a_{rs}| χ_{rs}^{a,J}(Q)

Internally, the calculator follows a pipeline:

  Step 1: For each (Q, k) and plane-wave G, extract plane wave coefficients of ρ̃ and do SVD
          over band indices → soft spectra f^(J)(Q, G, k).
  Step 2: Diagonalize each on-site C^a tensor → eigenvalues ε^a and eigenvectors O.
  Step 3: Contract D^a with O, then SVD → PAW spectra f^{a,J}(Q, rs, k).
  Step 4: Accumulate ξ from soft spectra and Coulomb factors v'(G+Q).
  Step 5: Accumulate χ from PAW spectra.
  Step 6: Diagonalize effective one-body matrix (with optional two-body shift).
  Step 7: Assemble λ from the three terms above.

k-point folding (k ⊕ Q):
  k points are represented in coefficient space c such that:
      k_cart = c @ B
  where B has rows equal to reciprocal primitive vectors (with 2π convention).
  We fold c(k)+c(Q) back into the supplied mesh window and map to an existing
  mesh index using rounded hashes.

Notes on performance:
  - The “soft” Fourier transform is implemented as a direct sum over grid points
    (not an FFT). This is simple and robust, but can be expensive for large grids.
  - Several expensive intermediates are cached. Use `force=True` to recompute.

"""

from __future__ import annotations

import numpy as np
from numpy.linalg import eigh

# -----------------------------------------------------------------------------
# Numeric knobs (module-level defaults)
# -----------------------------------------------------------------------------

# When we fold k+Q in coefficient space and round to find a mesh point, we allow
# a small “window slack” so points lying on the boundary fold consistently.
_COEFF_WINDOW_TOL = 1e-8

# When rescaling matrices by Frobenius norm before SVD, we avoid division by 0
# by clamping the norm from below.
_SCALE_FLOOR = 1e-14

# In PAW mode weighting we take sqrt(|ε|); clamp |ε| from below to avoid 0.
_EPS_FLOOR = 1e-15


class OneNormCalculator:
    """
    Compute the Hamiltonian one-norm λ for the Bloch–UPAW LCU decomposition.

    The total one-norm is (paper Eq. 37):

        λ = Σ_k Σ_i |ε_i(k)|
          + (1/4) Σ_J Σ_Q [ Σ_G ξ^(J)_G(Q)
              + Σ_a Σ_{r≤s} |ε^a_{rs}| χ^{a,J}_{rs}(Q) ]

    Inputs (high level):
      - k_mesh: Cartesian k points (Nk, 3), closed under folding k ⊕ Q.
      - rho:   smooth pseudo pair-density in real space on a uniform grid:
               (Nk, Nb, Nk, Nb, Nx, Ny, Nz)
      - Ca_dict, Da_dict: PAW augmentation tensors.
      - h_pq and two-body integrals (either rank-8 kappa_pqrs or reduced blocks).

    """

    def __init__(
        self,
        k_mesh,
        L_size,
        a_vectors,
        N_pw,
        rho,
        Ca_dict,
        Da_dict,
        h_pq=None,
        kappa_pqrs=None,
        # Reduced block-diagonal two-body pieces (preferred, if available)
        kappa_A=None,
        kappa_B=None,
        sv_floor=1e-13,
        scale_floor=1e-14,
        thr_rank=1e-14,
        k_round_decimals=10,
    ):
        # -------------------------
        # Basic inputs / thresholds
        # -------------------------
        self.k_mesh_input = np.asarray(k_mesh, dtype=float)  # (Nk,3) Cartesian
        self.Nk = self.k_mesh_input.shape[0]
        self.thr_rank = float(thr_rank)
        self.sv_floor = float(sv_floor)
        self.scale_floor = float(scale_floor)
        self.k_round_decimals = int(k_round_decimals)

        # -------------------------
        # Real/reciprocal lattice
        # -------------------------
        self.Lx, self.Ly, self.Lz = L_size
        self.a = np.asarray(a_vectors, dtype=float)  # (3,3) primitive real-space lattice vectors

        # Supercell lattice vectors as rows: A0=Lx*a0, etc.
        self.A = np.vstack(
            [self.Lx * self.a[0], self.Ly * self.a[1], self.Lz * self.a[2]]
        )  # (3,3)
        self.V = abs(np.linalg.det(self.A))

        # Reciprocal basis B with rows G_i (2π convention). Coefficients c satisfy k = c @ B.
        self.B = self._reciprocal_basis(self.A)  # (3,3)
        self.Binv = np.linalg.inv(self.B)

        # -------------------------
        # k-mesh in coefficient space
        # -------------------------
        self.k_cart = self.k_mesh_input  # (Nk,3)
        self.c_mesh = self.k_cart @ self.Binv  # (Nk,3)

        # Per-axis coefficient window derived from the supplied mesh
        self.c_lo = np.min(self.c_mesh, axis=0) - _COEFF_WINDOW_TOL
        self.c_hi = np.max(self.c_mesh, axis=0) + _COEFF_WINDOW_TOL

        # Hash table: rounded coefficient triple -> mesh index (for fast k ⊕ Q mapping)
        self._coeff_to_index = {}
        for idx, c in enumerate(self.c_mesh):
            key = tuple(np.round(c, self.k_round_decimals))
            if key in self._coeff_to_index:
                raise ValueError(f"Duplicate k-point after rounding: {key}")
            self._coeff_to_index[key] = idx

        # Precompute the folded addition table: kp = k ⊕ Q for all (Q,k)
        self._build_k_add_table()

        # -------------------------
        # Core tensors / integrals
        # -------------------------
        self.N_pw = int(N_pw)  # number of retained plane waves (excluding G=0)
        if self.N_pw <= 0:
            raise ValueError("N_pw must be positive (number of plane-wave modes, excluding G=0).")

        self.rho = rho
        self.Ca_dict = Ca_dict
        self.Da_dict = Da_dict
        self.h_pq = h_pq

        # Two-body: prefer reduced representation if provided; keep rank-8 for compatibility
        self.kappa_pqrs = kappa_pqrs
        self.kappa_A = kappa_A
        self.kappa_B = kappa_B

        # -------------------------
        # Dimensions inferred from rho
        # -------------------------
        # rho shape: (Nk, Nb, Nk, Nb, Nx, Ny, Nz)
        self.Nb = rho.shape[1]
        self.grid_shape = rho.shape[-3:]  # (Nx,Ny,Nz)
        self.Nx, self.Ny, self.Nz = map(int, self.grid_shape)
        self.N_grid = int(self.Nx * self.Ny * self.Nz)

        # Build a fixed plane-wave set {G} (G != 0), then build real-space grid
        self.G_cart, self.n_pw = self._build_pw_basis()
        self.Npw = self.N_pw  # alias used throughout
        self._r_cart = self._build_real_space_grid()  # (N_grid,3)

        # -------------------------
        # Caches (expensive intermediates)
        # -------------------------
        self._soft_eigenvalues = None       # dict[J] -> (Nk_Q, Npw, Nk, 2*Nb)
        self._c_tensor_eigen = None         # dict[a] -> {"eps","O","pairs"}
        self._paw_eigenvalues = None        # dict[(a,J)] -> (Nk_Q, P, Nk, 2*Nb)
        self._one_body_eigenvalues = None   # (Nk, Nb)

        self._rho_fft_g = None              # reserved (not used currently)
        self._coulomb_kernel = None         # (Npw,)
        self._coulomb_kernel_shifted = None # (Nk, Npw)
        self._abs_c_eigenvalues = None      # dict[a] -> (na,na)
        self._xi = None                     # dict[J] -> (Nk_Q, Npw)
        self._chi = None                    # dict[(a,J)] -> (Nk_Q, na, na)
        self._c_eigen_signs = None          # dict[a] -> sign(eps_rs)

    @classmethod
    def from_reader(
        cls,
        reader,
        *,
        scale_floor=1e-14,
        sv_floor=1e-12,
        thr_rank=1e-5,
        kappa_pqrs=None,
        kappa_A=None,
        kappa_B=None,
        k_round_decimals=10,
    ):
        """
        Convenience constructor: build from a PawReader (or compatible object).

        The reader can provide attributes directly or a `.load()` interface.
        Two-body integrals are accepted either as:
          - reduced blocks (kappa_A, kappa_B) or
          - full rank-8 tensor (kappa_pqrs / kappa / kappa_pqrs-like).
        """

        def unpack(r):
            base = ("k_mesh", "L_size", "a_vectors", "N_pw", "rho", "Ca_dict", "Da_dict")

            if all(hasattr(r, a) for a in base):
                h = getattr(r, "h_pq", None)
                kA = getattr(r, "kappa_A", None)
                kB = getattr(r, "kappa_B", None)
                kappa = getattr(r, "kappa", None)
                if kappa is None:
                    kappa = getattr(r, "kappa_pqrs", None)

                return (r.k_mesh, r.L_size, r.a_vectors, r.N_pw, r.rho,
                        r.Ca_dict, r.Da_dict, h, kappa, kA, kB)

            if hasattr(r, "load"):
                k_mesh, L_size, a_vectors, N_pw, rho, Ca_dict, Da_dict, *_ = r.load()

                h = None
                if hasattr(r, "load_one_body"):
                    try:
                        h, *_ = r.load_one_body()
                    except Exception:
                        h = None

                kA = getattr(r, "kappa_A", None)
                kB = getattr(r, "kappa_B", None)

                kappa = None
                if hasattr(r, "load_two_body"):
                    try:
                        kappa, *_ = r.load_two_body()
                    except Exception:
                        kappa = None

                return (k_mesh, L_size, a_vectors, N_pw, rho, Ca_dict, Da_dict, h, kappa, kA, kB)

            raise TypeError("Reader must expose attributes or .load().")

        (
            k_mesh,
            L_size,
            a_vectors,
            N_pw,
            rho,
            Ca_dict,
            Da_dict,
            h_pq_found,
            kappa_found,
            kA_found,
            kB_found,
        ) = unpack(reader)

        # Prefer explicitly provided reduced blocks, else fall back to reader
        kA = kappa_A if (kappa_A is not None) else kA_found
        kB = kappa_B if (kappa_B is not None) else kB_found

        if kA is not None:
            kA = np.asarray(kA, dtype=np.complex128)
            if kA.ndim != 5:
                raise ValueError(f"kappa_A must have shape (Nk,Nb,Nb,Nk,Nb); got {kA.shape}")
        if kB is not None:
            kB = np.asarray(kB, dtype=np.complex128)
            if kB.ndim != 5:
                raise ValueError(f"kappa_B must have shape (Nk,Nb,Nb,Nk,Nb); got {kB.shape}")

        # If reduced blocks aren't available, allow legacy rank-8 tensor
        kappa = kappa_pqrs if (kappa_pqrs is not None) else kappa_found
        if (kA is None or kB is None) and (kappa is not None):
            kappa = np.asarray(kappa, dtype=np.complex128)
            if kappa.ndim != 8:
                raise ValueError(
                    f"kappa_pqrs must be rank-8 (Nk,Nb,Nk,Nb,Nk,Nb,Nk,Nb); got {kappa.shape}"
                )
        else:
            kappa = None

        L_size = tuple(int(x) for x in L_size)
        N_pw = int(N_pw)

        return cls(
            k_mesh,
            L_size,
            a_vectors,
            N_pw,
            rho,
            Ca_dict,
            Da_dict,
            h_pq_found,
            kappa,
            kappa_A=kA,
            kappa_B=kB,
            sv_floor=float(sv_floor),
            scale_floor=float(scale_floor),
            thr_rank=float(thr_rank),
            k_round_decimals=int(k_round_decimals),
        )

    # ----------------------------------------------------------------------
    # LATTICE / GRIDS
    # ----------------------------------------------------------------------
    @staticmethod
    def _reciprocal_basis(A):
        """
        Build reciprocal basis B with ROWS = G_i (2π convention).
        With this convention:
            k_cart = c @ B
        """
        vol = np.linalg.det(A)
        b1 = 2 * np.pi * np.cross(A[1], A[2]) / vol
        b2 = 2 * np.pi * np.cross(A[2], A[0]) / vol
        b3 = 2 * np.pi * np.cross(A[0], A[1]) / vol
        return np.vstack([b1, b2, b3])

    def _build_real_space_grid(self) -> np.ndarray:
        """
        Construct uniform real-space grid points r (Cartesian, Å) spanning the supercell.
        The grid matches the (Nx,Ny,Nz) sampling of ρ̃.
        """
        ix = np.arange(self.Nx, dtype=np.float64) / float(self.Nx)
        iy = np.arange(self.Ny, dtype=np.float64) / float(self.Ny)
        iz = np.arange(self.Nz, dtype=np.float64) / float(self.Nz)

        u1, u2, u3 = np.meshgrid(ix, iy, iz, indexing="ij")
        u = np.stack([u1, u2, u3], axis=-1).reshape(-1, 3)  # (N_grid,3)

        # A has rows A0,A1,A2, so r = u_x*A0 + u_y*A1 + u_z*A2
        r = u @ self.A
        return r.astype(np.float64)

    def _build_pw_basis(self):
        """
        Choose the N_pw smallest-|G| reciprocal lattice vectors (excluding G=0).

        Returns:
            G_cart: (Npw,3) Cartesian reciprocal vectors (Å^-1)
            n_pw:   (Npw,3) integer coefficients in the reciprocal basis (rows of B)
        """
        nmax = 1
        while True:
            n = np.arange(-nmax, nmax + 1, dtype=np.int32)
            n1, n2, n3 = np.meshgrid(n, n, n, indexing="ij")
            nn = np.stack([n1, n2, n3], axis=-1).reshape(-1, 3)

            G = nn @ self.B
            G2 = np.einsum("ij,ij->i", G, G, optimize=True)

            # exclude G=0
            mask = G2 > 0.0
            nn = nn[mask]
            G = G[mask]
            G2 = G2[mask]

            if G.shape[0] >= self.N_pw:
                # Deterministic ordering: primary |G|^2, then lexicographic n
                order = np.lexsort((nn[:, 2], nn[:, 1], nn[:, 0], G2))
                sel = order[: self.N_pw]
                n_pw = nn[sel]
                G_cart = G[sel]
                break

            nmax *= 2
            if nmax > 2**12:
                raise RuntimeError("Failed to build PW basis: search box grew too large.")

        return G_cart.astype(np.float64), n_pw.astype(int)

    # ----------------------------------------------------------------------
    # k ⊕ Q folding infrastructure
    # ----------------------------------------------------------------------
    def _build_k_add_table(self):
        """
        Precompute kp = (k ⊕ Q) for all (Q, k) pairs.

        Result:
            self.k_add has shape (Nk, Nk) and stores mesh indices.
            Access as kp = self.k_add[Q_index, k_index].
        """
        tbl = np.empty((self.Nk, self.Nk), dtype=np.int32)
        for k in range(self.Nk):
            for q in range(self.Nk):
                tbl[q, k] = self._k_add_mod(k, q)
        self.k_add = tbl

    def _fold_coeff(self, c_sum: np.ndarray):
        """
        Fold coefficient-space vector into the mesh window [c_lo, c_hi] componentwise.

        Returns:
            c_fold: folded coefficient vector
            ell:    integer shift applied (bookkeeping; not used downstream)
        """
        c_fold = np.array(c_sum, dtype=float, copy=True)
        ell = np.zeros(3, dtype=int)
        for i in range(3):
            if c_fold[i] < self.c_lo[i]:
                t = int(np.ceil(self.c_lo[i] - c_fold[i]))
                c_fold[i] += t
                ell[i] -= t
            elif c_fold[i] > self.c_hi[i]:
                t = int(np.ceil(c_fold[i] - self.c_hi[i]))
                c_fold[i] -= t
                ell[i] += t
        return c_fold, ell

    def _k_add_mod(self, k_idx, q_idx):
        """
        Compute (k ⊕ Q) by:
          1) adding coefficient vectors
          2) folding back into the mesh window
          3) mapping to an existing mesh point via rounding + hash lookup
        """
        c_sum = self.c_mesh[k_idx] + self.c_mesh[q_idx]
        c_fold, _ = self._fold_coeff(c_sum)
        key = tuple(np.round(c_fold, self.k_round_decimals))
        try:
            return self._coeff_to_index[key]
        except KeyError:
            raise KeyError(
                "Folded (k ⊕ Q) is not a member of the supplied k-mesh. "
                "Ensure the mesh is closed under lattice-vector wrapping, "
                f"or relax rounding. Missing coeff key={key}"
            )

    # ----------------------------------------------------------------------
    # Coulomb kernels and PAW weights helpers (cached)
    # ----------------------------------------------------------------------
    def _vprime_if_needed(self, force: bool = False):
        """
        Regularized Coulomb kernel on the retained PW grid:
            v'(G) = 4π / |G|^2  (with v'(0)=0, though G=0 is excluded here)
        """
        if (self._coulomb_kernel is not None) and not force:
            return self._coulomb_kernel

        G_norm = np.linalg.norm(self.G_cart, axis=1)
        vprime = np.zeros(self.Npw, dtype=np.float64)
        nz = G_norm > 0
        vprime[nz] = 4 * np.pi / (G_norm[nz] ** 2)
        self._coulomb_kernel = vprime
        return self._coulomb_kernel

    def _vprime_shifted_if_needed(self, force: bool = False):
        """
        Shifted Coulomb kernel:
            v'(G+Q) = 4π / |G+Q|^2
        computed for all Q in the mesh and retained plane waves.
        """
        if (self._coulomb_kernel_shifted is not None) and not force:
            return self._coulomb_kernel_shifted

        Q_cart = self.k_cart
        K = Q_cart[:, None, :] + self.G_cart[None, :, :]  # (Nk, Npw, 3)
        Knorm = np.linalg.norm(K, axis=-1)

        vQG = np.zeros_like(Knorm, dtype=np.float64)
        nz = Knorm > 0
        vQG[nz] = 4 * np.pi / (Knorm[nz] ** 2)

        self._coulomb_kernel_shifted = vQG
        return vQG

    def _abs_eps_a_if_needed(self, force: bool = False):
        """
        Build dense |ε^a_{rs}| matrices (na x na) from the eigenvalues of C^a.
        This is used when assembling the PAW term in λ.
        """
        if (self._abs_c_eigenvalues is not None) and not force:
            return self._abs_c_eigenvalues
        if self._c_tensor_eigen is None:
            self.diagonalize_Ca()

        AbsE = {}
        for a, eigs in self._c_tensor_eigen.items():
            na = self.Ca_dict[a].shape[0]
            mat = np.zeros((na, na), dtype=np.float64)
            for idx, (r, s) in enumerate(eigs["pairs"]):
                e = abs(eigs["eps"][idx])
                mat[r, s] = mat[s, r] = e
            AbsE[a] = mat

        self._abs_c_eigenvalues = AbsE
        return self._abs_c_eigenvalues

    @staticmethod
    def _assemble_Ca_weights_all(D, O):
        """
        Vectorized contraction for PAW:
          - Flatten D's last two indices into P=na*na
          - Multiply by eigenvector matrix O (P,P)

        Args:
            D: (Nk, Nb, Nk, Nb, na, na)
            O: (P, P)

        Returns:
            (Nk, Nb, Nk, Nb, P) complex array
        """
        Nk, Nb, _, _, na, _ = D.shape
        P = na * na
        Dflat = np.asarray(D, dtype=np.complex128, copy=False).reshape(Nk, Nb, Nk, Nb, P)
        return np.tensordot(Dflat, O, axes=([4], [0]))

    # ----------------------------------------------------------------------
    # Step 1: soft eigenvalues (direct Fourier sum + batched SVD)
    # ----------------------------------------------------------------------
    def compute_soft_modes(self, force: bool = False, g_batch: int = 256):
        """
        Compute soft (plane-wave) spectra f^(J)(Q, G, k).

        For each transferred momentum Q (indexed by q) and each k:
          - pick kp = k ⊕ Q
          - take the real-space pair-density block ρ̃(k, kp)
          - compute Fourier coefficients at (G+Q) for selected plane waves
          - do SVD over (Nb x Nb) band block to get singular values σ
          - build eigenvalues ±(1/2)σ as required by the decomposition

        Returns:
            dict {0: f0, 1: f1} where fJ has shape (Nk, Npw, Nk, 2*Nb)
        """
        if not force and getattr(self, "_soft_eigenvalues", None) is not None:
            return self._soft_eigenvalues

        Nk, Nb, thr_rank = self.Nk, self.Nb, self.thr_rank
        sv_floor = self.sv_floor
        r = self._r_cart
        invNg = 1.0 / float(self.N_grid)

        f = {
            0: np.zeros((Nk, self.Npw, Nk, 2 * Nb), dtype=np.float64),
            1: np.zeros((Nk, self.Npw, Nk, 2 * Nb), dtype=np.float64),
        }

        for q in range(Nk):
            Q_cart = self.k_cart[q]
            for k in range(Nk):
                kp = self.k_add[q, k]

                # Pull only the needed (k,kp) block from rho (works with ndarray or h5py-like)
                rho_block = np.asarray(self.rho[k, :, kp, :, :, :, :])  # (Nb, Nb, Nx, Ny, Nz)
                rho_flat = rho_block.reshape(Nb, Nb, self.N_grid)

                # Compute Fourier coefficients for batches of plane waves
                for g0 in range(0, self.Npw, int(g_batch)):
                    g1 = min(self.Npw, g0 + int(g_batch))
                    G_batch = self.G_cart[g0:g1]                  # (B,3)
                    K_batch = G_batch + Q_cart[None, :]           # (B,3)

                    # phase(r,G) = exp(-i (G+Q)·r) / N_grid
                    phase = np.exp(-1j * (r @ K_batch.T)) * invNg  # (N_grid,B)

                    # C_ij(G) = sum_r rho_ij(r) * phase(r,G)
                    C_batch = np.tensordot(rho_flat, phase, axes=([2], [0]))  # (Nb,Nb,B)

                    # SVD wants leading batch dimension: (B, Nb, Nb)
                    X = np.ascontiguousarray(np.transpose(C_batch, (2, 0, 1)))

                    # Rescale each (Nb x Nb) block by its Frobenius norm for stability
                    sF = np.linalg.norm(X, ord="fro", axis=(-1, -2))  # (B,)
                    scale = np.maximum(sF, self.scale_floor)
                    Xn = X / scale[:, None, None]

                    _, S, _ = np.linalg.svd(Xn, full_matrices=False)  # (B, Nb)
                    S = S * scale[:, None]
                    S[S < sv_floor] = 0.0

                    # eigenvalues are ±(1/2)σ, with negatives reversed to match ordering
                    evals = np.concatenate([-0.5 * S[..., ::-1], 0.5 * S], axis=-1)  # (B, 2Nb)

                    f[0][q, g0:g1, k, :] = evals
                    f[1][q, g0:g1, k, :] = evals

        # Optional small-value truncation
        if thr_rank > 0:
            for J in (0, 1):
                f[J][np.abs(f[J]) < thr_rank] = 0.0

        self._soft_eigenvalues = f
        return f

    # ----------------------------------------------------------------------
    # Step 2: diagonalize on-site C^a tensors
    # ----------------------------------------------------------------------
    def diagonalize_Ca(self, force: bool = False):
        """
        Diagonalize each C^a tensor (na,na,na,na) as a (P x P) matrix with P=na^2.

        Stores:
          - eps: eigenvalues ε^a (length P)
          - O:   eigenvectors (P x P)
          - pairs: mapping from flat index -> (r,s)
        """
        if (self._c_tensor_eigen is not None) and not force:
            return self._c_tensor_eigen

        out = {}
        for a, Ca in self.Ca_dict.items():
            na = Ca.shape[0]
            P = na * na
            M = np.asarray(Ca, dtype=float, copy=False).reshape(P, P)
            eps, O = eigh(M)
            pairs = [(p, q) for p in range(na) for q in range(na)]
            out[a] = {"eps": eps.astype(np.float64, copy=False), "O": O, "pairs": pairs}

        self._c_tensor_eigen = out

        # Invalidate caches that depend on C^a eigen-decomposition
        self._abs_c_eigenvalues = None
        self._paw_eigenvalues = None
        self._chi = None
        self._c_eigen_signs = None

        return out

    # ----------------------------------------------------------------------
    # Step 3: PAW eigenvalues (D contraction + batched SVD)
    # ----------------------------------------------------------------------
    def compute_paw_modes(self, eps_floor=_EPS_FLOOR, force: bool = False):
        """
        Compute PAW augmentation spectra f^{a,J}(Q, rs, k).

        For each atom a:
          - diagonalize C^a to get ε_rs and eigenvectors O
          - contract D^a with O to get C_all(...)
          - for each Q:
              gather paired (k, k⊕Q) blocks
              weight by sqrt(|ε|) and take SVD to get singular values σ
              map to eigenvalues ±(1/2)σ
        """
        sv_floor = self.sv_floor
        thr_rank = self.thr_rank

        if (self._paw_eigenvalues is not None) and not force:
            if self._c_eigen_signs is None:
                self._c_eigen_signs = {
                    a: np.sign(eigs["eps"]) for a, eigs in self._c_tensor_eigen.items()
                }
            return self._paw_eigenvalues, self._c_eigen_signs

        Ca_eigs = self.diagonalize_Ca(force=False)
        Nk, Nb = self.Nk, self.Nb
        f_paw = {}
        eps_sign = {}

        k_idx = np.arange(Nk)

        for a, D in self.Da_dict.items():
            eigs = Ca_eigs[a]
            O = eigs["O"]  # (P,P)
            eps = np.asarray(eigs["eps"], dtype=np.float64)
            P = len(eps)
            eps_sign[a] = np.sign(eps)

            # C_all: (Nk, Nb, Nk, Nb, P)
            C_all = self._assemble_Ca_weights_all(D, O)

            for J in (0, 1):
                fp = np.zeros((Nk, P, Nk, 2 * Nb), dtype=np.float64)

                for q in range(Nk):
                    kp_all = self.k_add[q]  # (Nk,) mapping k -> k⊕Q

                    # Pairwise gather: (Nk, Nb, Nb, P) for (k, kp[k]) blocks
                    C_qk = C_all[k_idx, :, kp_all, :, :]

                    # Reorder to (Nk, P, Nb, Nb) for batched SVD
                    C_qk = np.transpose(C_qk, (0, 3, 1, 2)).astype(np.complex128, copy=False)
                    C_qk = np.ascontiguousarray(C_qk)

                    # Weight by sqrt(|ε|) (with a small floor for stability)
                    amp = np.sqrt(np.maximum(np.abs(eps), eps_floor))[None, :, None, None]
                    X = amp * C_qk  # (Nk, P, Nb, Nb)

                    # Rescale each block before SVD
                    sF = np.linalg.norm(X, ord="fro", axis=(-1, -2))                 # (Nk,P)
                    scale = np.maximum(sF, self.scale_floor)[..., None, None]        # (Nk,P,1,1)
                    Xn = np.ascontiguousarray(X / scale)

                    _, S, _ = np.linalg.svd(Xn, full_matrices=False)                 # (Nk,P,Nb)
                    S = (S * scale[..., 0, 0][..., None]).astype(np.float64, copy=False)
                    S[S < sv_floor] = 0.0

                    pos = 0.5 * S
                    neg = -pos[..., ::-1]
                    evals = np.concatenate([neg, pos], axis=-1)                      # (Nk,P,2Nb)

                    # Store in fp with q leading: fp[q, rs, k, i]
                    fp[q, :, :, :] = np.transpose(evals, (1, 0, 2))                  # (P, Nk, 2Nb)

                if thr_rank > 0:
                    fp[np.abs(fp) < thr_rank] = 0.0
                f_paw[(a, J)] = fp

        self._paw_eigenvalues = f_paw
        self._c_eigen_signs = eps_sign
        self._chi = None  # dependent cache

        return f_paw, eps_sign

    # ----------------------------------------------------------------------
    # Step 4: ξ(Q,G) from soft spectra
    # ----------------------------------------------------------------------
    def compute_xi(self, force: bool = False):
        """
        Compute ξ^(J)_G(Q) used in the soft part of λ.

        Implementation follows:
          S(Q,G) = Σ_{k,i} | f^(J)_i(Q,G,k) |
          ξ(Q,G) = v'(G+Q) * [ S(Q,G) ]^2
        """
        if (self._xi is not None) and not force:
            return self._xi

        f_soft = self.compute_soft_modes(force=False)
        xi = {}

        # v'(G+Q) for all Q and retained G
        vprime_QG = self._vprime_shifted_if_needed(force=force)

        for J in (0, 1):
            fJ = f_soft[J]                          # (Nk_Q, Npw, Nk, 2Nb)
            S = np.sum(np.abs(fJ), axis=(2, 3))     # (Nk_Q, Npw)
            xi[J] = vprime_QG * (S ** 2)

        self._xi = xi
        return xi

    # ----------------------------------------------------------------------
    # Step 5: χ(Q,rs) from PAW spectra
    # ----------------------------------------------------------------------
    def compute_chi(self, force: bool = False):
        """
        Compute χ^{a,J}_{rs}(Q) used in the PAW part of λ.

        Implementation:
          S(Q,rs) = Σ_{k,i} | f^{a,J}_i(Q,rs,k) |
          χ(Q,rs) = [ S(Q,rs) ]^2
        Returned as dense (Nk_Q, na, na) matrices for each (a,J).
        """
        if (self._chi is not None) and not force:
            return self._chi

        f_paw, _ = self.compute_paw_modes(force=False)
        chi = {}

        for (a, J), fp in f_paw.items():
            P = fp.shape[1]
            S = np.sum(np.abs(fp), axis=(2, 3))  # (Nk_Q, P)

            na = self.Ca_dict[a].shape[0]
            M = np.zeros((self.Nk, na, na), dtype=np.float64)
            pairs = self._c_tensor_eigen[a]["pairs"]

            # Unflatten rs index back into (r,s) matrix
            for idx, (r, s) in enumerate(pairs):
                M[:, r, s] = S[:, idx]
                M[:, s, r] = S[:, idx]

            chi[(a, J)] = M ** 2

        self._chi = chi
        return chi

    # ----------------------------------------------------------------------
    # Step 6: one-body eigenvalues ε_i(k)
    # ----------------------------------------------------------------------
    def diagonalize_one_plus_two_body(self, force: bool = False):
        """
        Compute ε_i(k) for the one-body term in λ.

        If two-body data are available, we apply the mean-field-like shift used
        by the referenced formula (see code for exact contractions). If not, we
        diagonalize h_pq(k,k) directly.

        Returns:
            eps: (Nk, Nb) array of eigenvalues
        """
        if (self._one_body_eigenvalues is not None) and not force:
            return self._one_body_eigenvalues

        Nk, Nb = self.Nk, self.Nb
        eps = np.zeros((Nk, Nb))

        # No one-body input => zero contribution
        if self.h_pq is None:
            self._one_body_eigenvalues = eps
            return eps

        have_reduced = (self.kappa_A is not None) and (self.kappa_B is not None)

        # If we have no two-body info, just diagonalize h(k)
        if (not have_reduced) and (self.kappa_pqrs is None):
            for k in range(Nk):
                Hk = np.array(self.h_pq[k, :, k, :], dtype=float)
                Hk = 0.5 * (Hk + Hk.T)  # enforce Hermiticity in floating arithmetic
                w, _ = eigh(Hk)
                eps[k, :] = w
            self._one_body_eigenvalues = eps
            return eps

        # Otherwise include the (termA - termB) shift
        for k in range(Nk):
            Hk = np.array(self.h_pq[k, :, k, :], dtype=float)
            Hk = 0.5 * (Hk + Hk.T)

            if have_reduced:
                # Shapes assumed by reduced representation:
                # kappa_A[k, i, j, l, ?] etc. (kept as in original code)
                termA = np.sum(self.kappa_A[k, :, :, :, :], axis=(2, 3))
                termB = 2.0 * np.sum(self.kappa_B[k, :, :, :, :], axis=(2, 3))
            else:
                # Legacy rank-8 contractions
                termA = np.einsum("ikljkl->ij", self.kappa_pqrs[k, :, :, :, k, :, :, :])
                termB = 2.0 * np.einsum("ijklkl->ij", self.kappa_pqrs[k, :, k, :, :, :, :, :])

            shift = termA - termB
            w, _ = eigh(Hk - 0.5 * shift)
            eps[k, :] = w

        self._one_body_eigenvalues = eps
        return eps

    # ----------------------------------------------------------------------
    # Step 7: assemble λ
    # ----------------------------------------------------------------------
    def lambda_one_norm(self):
        """
        Compute the full Hamiltonian one-norm λ.

        Returns:
            float λ
        """
        xi = self.compute_xi()
        chi = self.compute_chi()
        eps_one = self.diagonalize_one_plus_two_body()

        term1 = np.sum(np.abs(eps_one))

        term2 = 0.0
        for J in (0, 1):
            term2 += 0.25 * np.sum(xi[J])

        term3 = 0.0
        AbsE = self._abs_eps_a_if_needed()
        for a in self._c_tensor_eigen:
            for J in (0, 1):
                # einsum('rs,Qrs->Q', AbsE, chi) then sum_Q
                term3 += 0.25 * np.sum(np.einsum("rs,Qrs->Q", AbsE[a], chi[(a, J)]))

        return float(term1 + term2 + term3)

    def compute_average_rank(self) -> float:
        """
        Compute average LCU rank indicators.

        Returns:
            (R_ell_not0_avg, R0_avg):
              - average rank across ℓ≠0 components
              - average one-body rank per k-point
        """
        thr_rank = self.thr_rank
        f_soft = self.compute_soft_modes()
        f_paw, _ = self.compute_paw_modes()
        eps = self.diagonalize_one_plus_two_body()

        Nk = self.Nk
        Npw = self.Npw

        # M counts “modes” per Q (soft plane-waves + PAW symmetric pairs)
        sum_pairs = 0.0
        for a in self.Ca_dict:
            na = self.Ca_dict[a].shape[0]
            sum_pairs += 0.5 * na * (na + 1)
        M = Npw + sum_pairs

        # L is the total number of ℓ≠0 terms (matches original code’s normalization)
        L = 2.0 * Nk * M

        rank_sum = 0
        for J in (0, 1):
            rank_sum += np.count_nonzero(np.abs(f_soft[J]) > thr_rank)
            for a in self.Ca_dict:
                rank_sum += np.count_nonzero(np.abs(f_paw[(a, J)]) > thr_rank)

        R0 = np.count_nonzero(np.abs(eps) > thr_rank)
        return float(rank_sum / L), float(R0 / Nk)