from __future__ import annotations

import numpy as np
import pytest

from bloch_paw.one_norm import OneNormCalculator


# ---------------------------------------------------------------------------
# Shared fixture: a small synthetic system (Nk=2, Nb=2, grid=3x3x3)
# ---------------------------------------------------------------------------

def _make_simple_cubic_system(
    Nk: int = 2,
    Nb: int = 2,
    grid: tuple[int, int, int] = (3, 3, 3),
    na: int = 2,
    seed: int = 42,
) -> dict:
    """Build a synthetic system on a simple-cubic lattice closed under k-point folding.

    The lattice constant is 3.0 Ang, with L_size = (Nk, 1, 1) so the
    first reciprocal direction is split into Nk k-points.
    """
    rng = np.random.default_rng(seed)
    a_lat = 3.0
    a_vectors = np.eye(3) * a_lat
    L_size = (Nk, 1, 1)

    A = np.diag([Nk * a_lat, a_lat, a_lat])
    vol = np.linalg.det(A)
    B = np.vstack([
        2 * np.pi * np.cross(A[1], A[2]) / vol,
        2 * np.pi * np.cross(A[2], A[0]) / vol,
        2 * np.pi * np.cross(A[0], A[1]) / vol,
    ])

    c_mesh = np.zeros((Nk, 3))
    for i in range(Nk):
        c_mesh[i, 0] = i / Nk
    k_mesh = c_mesh @ B

    Nx, Ny, Nz = grid
    N_pw = Nx * Ny * Nz

    rho = (rng.standard_normal((Nk, Nb, Nk, Nb, Nx, Ny, Nz))
           + 1j * rng.standard_normal((Nk, Nb, Nk, Nb, Nx, Ny, Nz)))

    P = na * na
    M = rng.standard_normal((P, P))
    M = 0.5 * (M + M.T)
    Ca = {0: M.reshape(na, na, na, na)}

    Da = {0: (rng.standard_normal((Nk, Nb, Nk, Nb, na, na))
              + 1j * rng.standard_normal((Nk, Nb, Nk, Nb, na, na)))}

    h_pq_mat = rng.standard_normal((Nk, Nb, Nk, Nb))
    h_pq = 0.5 * (h_pq_mat + h_pq_mat.transpose(2, 3, 0, 1))

    kappa_pqrs = (rng.standard_normal((Nk, Nb, Nk, Nb, Nk, Nb, Nk, Nb)) * 0.01
                  + 1j * rng.standard_normal((Nk, Nb, Nk, Nb, Nk, Nb, Nk, Nb)) * 0.01)

    return dict(
        k_mesh=k_mesh,
        L_size=L_size,
        a_vectors=a_vectors,
        N_pw=N_pw,
        rho=rho,
        Ca_dict=Ca,
        Da_dict=Da,
        h_pq=h_pq,
        kappa_pqrs=kappa_pqrs,
    )


@pytest.fixture
def small_system():
    return _make_simple_cubic_system()


@pytest.fixture
def calculator(small_system):
    return OneNormCalculator(**small_system)


# ---------------------------------------------------------------------------
# Construction and lattice setup
# ---------------------------------------------------------------------------

class TestConstruction:
    """Verify that the constructor correctly builds all derived lattice objects."""

    def test_basic_attributes(self, calculator, small_system):
        """Test that Nk, Nb, grid_shape, and N_pw are inferred from input arrays.

        These dimensions propagate to every downstream computation — wrong
        values silently produce wrong-shaped intermediates.
        """
        assert calculator.Nk == 2
        assert calculator.Nb == 2
        assert calculator.grid_shape == (3, 3, 3)
        assert calculator.N_pw == 27

    def test_lattice_and_reciprocal_vectors(self, calculator):
        """Test that real/reciprocal lattice vectors and k-point tables have correct shapes.

        The supercell A, reciprocal B, coefficient mesh c_mesh, and k-addition
        table k_add are all built in __init__ and must be self-consistent.
        """
        assert calculator.a.shape == (3, 3)
        assert calculator.A.shape == (3, 3)
        assert calculator.B.shape == (3, 3)
        assert calculator.V > 0
        assert calculator.c_mesh.shape == (2, 3)
        assert calculator.k_add.shape == (2, 2)

    def test_G_grid_and_masking(self, calculator):
        """Test that the FFT G-grid masks out exactly the G=0 origin.

        The one-norm formula (Eq. 37) sums over G != 0. If G=0 is not
        masked out, the Coulomb kernel 1/|G|^2 diverges.
        """
        N = 3 * 3 * 3
        assert calculator.G_cart.shape == (N, 3)
        assert calculator.n_ix.shape == (N, 3)
        assert calculator.G_mask.shape == (N,)
        assert calculator.G_norm.shape == (N,)
        # Exactly one G-vector (the origin) should be masked out
        assert np.sum(~calculator.G_mask) == 1
        zero_idx = np.where(~calculator.G_mask)[0][0]
        np.testing.assert_allclose(calculator.G_cart[zero_idx], 0.0, atol=1e-14)
        assert calculator.Npw == N - 1


# ---------------------------------------------------------------------------
# reciprocal_basis
# ---------------------------------------------------------------------------

class TestReciprocalBasis:
    """Verify the B = 2pi (A^T)^{-1} reciprocal basis computation."""

    def test_orthogonality_with_real_space(self, calculator):
        """Test that B @ A^T = 2pi I for the constructed supercell.

        This identity is the defining property of reciprocal lattice vectors
        and is assumed by all k-point folding operations.
        """
        product = calculator.B @ calculator.A.T
        np.testing.assert_allclose(product, 2 * np.pi * np.eye(3), atol=1e-12)

    def test_simple_cubic_values(self):
        """Test that a 5A cubic cell gives diagonal reciprocal vectors 2pi/5.

        A simple analytic case where the answer is known exactly.
        """
        a = 5.0
        A = np.eye(3) * a
        B = OneNormCalculator._reciprocal_basis(A)
        expected = np.eye(3) * (2 * np.pi / a)
        np.testing.assert_allclose(B, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# fold_coeff and k_add_mod
# ---------------------------------------------------------------------------

class TestKFolding:
    """Test the k-point folding algebra that maps (k+Q) back onto the mesh."""

    def test_fold_identity(self, calculator):
        """Test that folding a point already in the mesh window returns it unchanged.

        The lattice shift vector ell should be exactly zero.
        """
        c = calculator.c_mesh[0].copy()
        c_fold, ell = calculator._fold_coeff(c)
        np.testing.assert_allclose(c_fold, c, atol=1e-12)
        np.testing.assert_array_equal(ell, [0, 0, 0])

    def test_fold_shift(self, calculator):
        """Test that adding 1.0 to a coefficient component folds back by one lattice vector.

        Folding periodicity is essential for the LCU decomposition — if the
        mesh is not closed under addition modulo G, lambda_one_norm diverges.
        """
        c = calculator.c_mesh[0].copy()
        c[0] += 1.0
        c_fold, ell = calculator._fold_coeff(c)
        np.testing.assert_allclose(c_fold, calculator.c_mesh[0], atol=1e-10)
        assert ell[0] != 0

    def test_k_add_mod_closure(self, calculator):
        """Test that every (k + q) mod G maps to a valid mesh index.

        The k-mesh must form a group under folded addition, otherwise
        the Q-summation in the one-norm formula is undefined.
        """
        Nk = calculator.Nk
        for k in range(Nk):
            for q in range(Nk):
                kp = calculator._k_add_mod(k, q)
                assert 0 <= kp < Nk

    def test_k_add_table_matches_k_add_mod(self, calculator):
        """Test that the precomputed k_add table agrees with _k_add_mod pointwise.

        The table is built once in __init__ for speed; this checks it was
        built correctly.
        """
        Nk = calculator.Nk
        for q in range(Nk):
            for k in range(Nk):
                assert calculator.k_add[q, k] == calculator._k_add_mod(k, q)


# ---------------------------------------------------------------------------
# Compute methods: output shapes and properties
# ---------------------------------------------------------------------------

class TestComputeOutputShapes:
    """Verify that the multi-step pipeline produces arrays of correct shape.

    Each compute method implements a specific equation from the paper.
    Shape mismatches between steps cause hard-to-diagnose downstream errors.
    """

    def test_compute_soft_modes(self, calculator):
        """Test that soft mode SVD produces f[J] of shape (Nk, Npw, Nk, 2*Nb).

        compute_soft_modes (Step 1) FFTs rho and does batched SVD. The output
        has two J components (J=0,1 for real/imaginary split).
        """
        f = calculator.compute_soft_modes()
        Nk, Nb, Npw = calculator.Nk, calculator.Nb, calculator.Npw
        for J in (0, 1):
            assert f[J].shape == (Nk, Npw, Nk, 2 * Nb)

    def test_diagonalize_Ca(self, calculator):
        """Test that C^a eigendecomposition returns real eigenvalues and eigenvectors.

        Step 2 eigendecomposes the rank-4 C^a tensor. Eigenvalues must be real
        (C^a is symmetric), and the 'pairs' list drives the PAW mode computation.
        """
        out = calculator.diagonalize_Ca()
        assert 0 in out
        info = out[0]
        assert "eps" in info
        assert "O" in info
        assert "pairs" in info
        assert info["eps"].dtype == np.float64

    def test_compute_paw_modes(self, calculator):
        """Test that PAW mode SVD produces f_paw[(atom, J)] of correct shape.

        Step 3 contracts D^a with the C^a eigenvectors and does SVD per pair.
        """
        f_paw, eps_sign = calculator.compute_paw_modes()
        Nk, Nb, na = calculator.Nk, calculator.Nb, 2
        P = na * na
        for J in (0, 1):
            assert (0, J) in f_paw
            assert f_paw[(0, J)].shape == (Nk, P, Nk, 2 * Nb)
        assert 0 in eps_sign

    def test_compute_xi(self, calculator):
        """Test that xi (soft contribution) has shape (Nk, Npw) per J.

        Step 4 accumulates the soft one-norm contribution from the SVD
        eigenvalues across all Q-points.
        """
        xi = calculator.compute_xi()
        Nk, Npw = calculator.Nk, calculator.Npw
        for J in (0, 1):
            assert xi[J].shape == (Nk, Npw)

    def test_compute_chi(self, calculator):
        """Test that chi (PAW contribution) has shape (Nk, na, na) per (atom, J).

        Step 5 accumulates the PAW on-site one-norm contribution.
        """
        chi = calculator.compute_chi()
        Nk, na = calculator.Nk, 2
        for J in (0, 1):
            assert (0, J) in chi
            assert chi[(0, J)].shape == (Nk, na, na)

    def test_diagonalize_one_plus_two_body(self, calculator):
        """Test that one+two-body eigenvalues have shape (Nk, Nb).

        Step 6 diagonalises h'(k) = h(k) + mean-field shift. The eigenvalues
        enter term 1 of the one-norm formula.
        """
        eps = calculator.diagonalize_one_plus_two_body()
        assert eps.shape == (calculator.Nk, calculator.Nb)


# ---------------------------------------------------------------------------
# lambda_one_norm end-to-end
# ---------------------------------------------------------------------------

class TestLambdaOneNorm:
    """Test the top-level lambda_one_norm() that assembles all three terms."""

    def test_lambda_properties(self, calculator):
        """Test that lambda is a finite, non-negative float.

        lambda is a physical quantity (Hamiltonian one-norm) that must be
        non-negative and finite for the block encoding to be valid.
        """
        lam = calculator.lambda_one_norm()
        assert isinstance(lam, float)
        assert lam >= 0.0
        assert np.isfinite(lam)


# ---------------------------------------------------------------------------
# compute_average_rank
# ---------------------------------------------------------------------------

class TestComputeAverageRank:
    """Test the SVD rank averaging used for resource estimation."""

    def test_returns_nonnegative_floats(self, calculator):
        """Test that average rank and R0 are non-negative floats.

        These values feed into the Toffoli count formula — negative or
        non-float values would break the resource estimator.
        """
        rank_avg, r0 = calculator.compute_average_rank()
        assert isinstance(rank_avg, float)
        assert isinstance(r0, float)
        assert rank_avg >= 0.0
        assert r0 >= 0.0


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------

class TestCaching:
    """Test that expensive computations are cached and respect the force flag."""

    def test_soft_modes_cached(self, calculator):
        """Test that repeated calls return the same cached object.

        compute_soft_modes is the most expensive step (~O(Nk^2 * Npw * Nb^2)).
        Caching avoids redundant SVDs when lambda_one_norm calls it multiple times.
        """
        f1 = calculator.compute_soft_modes()
        f2 = calculator.compute_soft_modes()
        assert f1 is f2

    def test_soft_modes_force_recompute(self, calculator):
        """Test that force=True bypasses the cache and returns fresh arrays.

        Needed when users modify rho in-place (e.g. applying a new threshold)
        and want to recompute without creating a new calculator.
        """
        f1 = calculator.compute_soft_modes()
        f2 = calculator.compute_soft_modes(force=True)
        assert f1 is not f2


# ---------------------------------------------------------------------------
# from_reader (with a minimal mock reader)
# ---------------------------------------------------------------------------

class TestFromReader:
    """Test the two from_reader() code paths (attribute-style and load-style)."""

    def test_from_reader_attribute_style(self, small_system):
        """Test construction from a reader that exposes data as attributes.

        This is the PawReader path — all fields are set as attributes
        after load().
        """
        class AttrReader:
            pass

        reader = AttrReader()
        for key, val in small_system.items():
            setattr(reader, key, val)

        calc = OneNormCalculator.from_reader(reader)
        assert calc.Nk == 2
        assert calc.Nb == 2

    def test_from_reader_load_style(self, small_system):
        """Test construction from a reader with a .load() method returning a tuple.

        This is the legacy API path used by older notebook code.
        """
        class LoadReader:
            def __init__(self, data):
                self._data = data

            def load(self):
                d = self._data
                return (d["k_mesh"], d["L_size"], d["a_vectors"],
                        d["N_pw"], d["rho"], d["Ca_dict"], d["Da_dict"], {})

        reader = LoadReader(small_system)
        calc = OneNormCalculator.from_reader(
            reader,
            kappa_pqrs=small_system["kappa_pqrs"],
        )
        assert calc.Nk == 2

    def test_from_reader_rejects_bad_kappa_rank(self, small_system):
        """Test that passing a non-rank-8 kappa array raises ValueError.

        kappa_pqrs must be (Nk,Nb)^4 = rank 8. A rank-4 array indicates
        the user forgot to include k-point indices.
        """
        class AttrReader:
            pass

        reader = AttrReader()
        for key in ("k_mesh", "L_size", "a_vectors", "N_pw", "rho", "Ca_dict", "Da_dict"):
            setattr(reader, key, small_system[key])
        reader.h_pq = small_system["h_pq"]
        reader.kappa_pqrs = None

        with pytest.raises(ValueError, match="rank-8"):
            bad_kappa = np.zeros((2, 2, 2, 2))
            OneNormCalculator.from_reader(reader, kappa_pqrs=bad_kappa)


# ---------------------------------------------------------------------------
# One-body fallback paths (kappa=None, h_pq=None)
# ---------------------------------------------------------------------------

class TestOneBodyFallback:
    """Test the three code paths in diagonalize_one_plus_two_body().

    Path 1 (full): kappa present -> build mean-field h'(k) from kappa
    Path 2 (h_pq only): kappa=None, h_pq present -> diagonalise h(k) directly
    Path 3 (neither): both None -> return zeros
    """

    @pytest.fixture
    def system_no_kappa(self, small_system):
        """System with h_pq but no kappa_pqrs."""
        data = dict(small_system)
        data["kappa_pqrs"] = None
        return data

    @pytest.fixture
    def system_no_h_pq(self, small_system):
        """System with neither h_pq nor kappa_pqrs."""
        data = dict(small_system)
        data["h_pq"] = None
        data["kappa_pqrs"] = None
        return data

    def test_fallback_eigenvalues_match_eigh_of_h_pq(self, system_no_kappa):
        """Test that without kappa, eigenvalues equal numpy.linalg.eigh(h_pq[k,:,k,:]).

        This is the critical numerical test for the fallback path — it verifies
        that the fallback correctly diagonalises the one-body block h(k).
        """
        calc = OneNormCalculator(**system_no_kappa)
        eps = calc.diagonalize_one_plus_two_body()

        h_pq = system_no_kappa["h_pq"]
        Nk = calc.Nk
        for k in range(Nk):
            Hk = np.array(h_pq[k, :, k, :], dtype=float)
            expected, _ = np.linalg.eigh(Hk)
            np.testing.assert_allclose(eps[k, :], expected, atol=1e-12)

    def test_no_h_pq_returns_zeros(self, system_no_h_pq):
        """Test that eigenvalues are all zero when both h_pq and kappa are absent.

        This is the degenerate case where the one-body contribution to lambda
        is exactly zero (only soft + PAW terms contribute).
        """
        calc = OneNormCalculator(**system_no_h_pq)
        eps = calc.diagonalize_one_plus_two_body()
        np.testing.assert_array_equal(eps, 0.0)

    def test_fallback_lambda_nonzero_and_finite(self, system_no_kappa):
        """Test that lambda is positive and finite when h_pq is present but kappa absent.

        A zero lambda would indicate the one-body and soft+PAW terms
        cancelled exactly, which is not physical for random data.
        """
        calc = OneNormCalculator(**system_no_kappa)
        lam = calc.lambda_one_norm()
        assert lam > 0.0
        assert np.isfinite(lam)

    def test_no_h_pq_lambda_nonneg(self, system_no_h_pq):
        """Test that lambda is non-negative and finite with no one-body data.

        The soft+PAW terms should still produce a well-defined one-norm.
        """
        calc = OneNormCalculator(**system_no_h_pq)
        lam = calc.lambda_one_norm()
        assert lam >= 0.0
        assert np.isfinite(lam)

    def test_full_vs_fallback_different(self, small_system, system_no_kappa):
        """Test that the full path (with kappa) and fallback (without) give different eigenvalues.

        The mean-field correction from kappa shifts the effective one-body
        eigenvalues. If full and fallback agree, either kappa is being
        ignored or the correction is zero (both would be bugs).
        """
        calc_full = OneNormCalculator(**small_system)
        calc_fallback = OneNormCalculator(**system_no_kappa)

        eps_full = calc_full.diagonalize_one_plus_two_body()
        eps_fallback = calc_fallback.diagonalize_one_plus_two_body()

        assert not np.allclose(eps_full, eps_fallback)

    def test_compute_average_rank_no_kappa(self, system_no_kappa):
        """Test that compute_average_rank works without kappa.

        Resource estimation needs rank data even when the two-body term
        is omitted (fast path).
        """
        calc = OneNormCalculator(**system_no_kappa)
        rank_avg, r0 = calc.compute_average_rank()
        assert isinstance(rank_avg, float)
        assert isinstance(r0, float)
        assert rank_avg >= 0.0
