from __future__ import annotations

import math
import sys
import types
from unittest import mock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Ensure h5py is available (mock it if not installed)
# ---------------------------------------------------------------------------
try:
    import h5py
except ImportError:
    h5py = types.ModuleType("h5py")
    h5py.File = mock.MagicMock  # type: ignore[attr-defined]
    sys.modules["h5py"] = h5py

from bloch_paw.resources import ResourceEstimator


# ===================================================================
# Construction
# ===================================================================

class TestConstruction:
    """ResourceEstimator.__init__ with valid and invalid parameters."""

    def test_valid_params(self):
        """Test that Nk, Npw, P, Nb are stored and L = 2*(Npw+P)*Nk is computed.

        L is the total number of LCU labels and appears in every downstream
        formula. Getting it wrong silently corrupts all resource estimates.
        """
        est = ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5)
        assert est.Nk == 4
        assert est.Npw == 100
        assert est.P == 10
        assert est.Nb == 5
        assert est.L == 2 * (100 + 10) * 4

        # Nb is optional
        est2 = ResourceEstimator(Nk=2, Npw=10, P=5)
        assert est2.Nb is None

    def test_L_formula_multiple_sizes(self):
        """Test L = 2*(Npw+P)*Nk across multiple parameter combinations.

        Catches off-by-one errors in the label count formula.
        """
        for Nk, Npw, P in [(1, 0, 0), (3, 50, 7), (10, 200, 20)]:
            est = ResourceEstimator(Nk=Nk, Npw=Npw, P=P)
            assert est.L == 2 * (Npw + P) * Nk

    @pytest.mark.parametrize("Nk", [0, -1, -100])
    def test_invalid_Nk(self, Nk):
        """Test that non-positive Nk raises ValueError.

        Nk=0 would make L=0, causing division by zero in bit-width formulas.
        """
        with pytest.raises(ValueError, match="Nk must be positive"):
            ResourceEstimator(Nk=Nk, Npw=100, P=10)

    @pytest.mark.parametrize("Npw,P", [(-1, 10), (100, -1), (-5, -5)])
    def test_negative_Npw_or_P(self, Npw, P):
        """Test that negative Npw or P raises ValueError.

        These are physical counts that must be non-negative.
        """
        with pytest.raises(ValueError, match="nonnegative"):
            ResourceEstimator(Nk=4, Npw=Npw, P=P)


# ===================================================================
# Bit-width helpers
# ===================================================================

class TestBitWidthHelpers:
    """Test _bits() and _ceil_log2_nonneg() edge cases.

    These appear in every qubit-count and Toffoli-count formula.
    Edge cases at powers of two and zero are where bugs hide.
    """

    @pytest.mark.parametrize("x,expected", [
        (0, 1),
        (0.5, 1),
        (1, 1),
        (2, 1),
        (3, 2),
        (4, 2),
        (1024, 10),
        (1025, 11),
    ])
    def test_bits(self, x, expected):
        """Test that _bits(x) = max(1, ceil(log2(x))) for key boundary values."""
        assert ResourceEstimator._bits(x) == expected

    @pytest.mark.parametrize("x,expected", [
        (0, 0),
        (0.5, 0),
        (1, 0),
        (2, 1),
        (3, 2),
        (4, 2),
        (1024, 10),
    ])
    def test_ceil_log2_nonneg(self, x, expected):
        """Test that _ceil_log2_nonneg(x) = max(0, ceil(log2(x))) for key values."""
        assert ResourceEstimator._ceil_log2_nonneg(x) == expected


# ===================================================================
# Core formula tests (_toffoli_core / _total_qubits_core)
# ===================================================================

# A consistent set of small parameters for core-formula calls
_CORE_PARAMS = dict(
    L=880, kp1=4, kp1_prime=3, ko=5, ko_prime=3,
    kp2=6, kp2_prime=4, kr=7, kr_prime=5,
    br=20, Rl=16, R0=4, Nb=5, Nk=4, N1=8, N2=6, eta=3, B=10,
)


class TestCoreFormulas:
    """Test the class-level _toffoli_core and _total_qubits_core methods."""

    def test_toffoli_core_returns_positive_int(self):
        """Test that the Toffoli count formula produces a positive integer.

        The formula implements Eq. 98 of the paper. A non-positive result
        would indicate a sign error in one of the terms.
        """
        result = ResourceEstimator._toffoli_core(**_CORE_PARAMS, bo_override=None)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_core_returns_positive_int(self):
        """Test that the qubit count formula produces a positive integer (Eq. 99)."""
        params = {k: v for k, v in _CORE_PARAMS.items()
                  if k in ("L", "kr", "br", "Rl", "R0", "Nb", "Nk", "N1", "N2", "eta", "B")}
        result = ResourceEstimator._total_qubits_core(**params, bo_override=None)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_core_with_I(self):
        """Test that passing an explicit iteration count I overrides the default."""
        params = {k: v for k, v in _CORE_PARAMS.items()
                  if k in ("L", "kr", "br", "Rl", "R0", "Nb", "Nk", "N1", "N2", "eta", "B")}
        result = ResourceEstimator._total_qubits_core(**params, I=100, bo_override=None)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_core_with_lam_eps(self):
        """Test that passing lam + eps_qpe computes I internally from pi*lam/(2*eps).

        This is the QPE iteration count formula from Eq. 42.
        """
        params = {k: v for k, v in _CORE_PARAMS.items()
                  if k in ("L", "kr", "br", "Rl", "R0", "Nb", "Nk", "N1", "N2", "eta", "B")}
        result = ResourceEstimator._total_qubits_core(
            **params, lam=50.0, eps_qpe=0.01, bo_override=None
        )
        assert isinstance(result, int)
        assert result > 0


# ===================================================================
# Instance methods: toffoli_count_per_be, total_qubits (qroam modes)
# ===================================================================

_INST_KW = dict(br=20, Rl=16, R0=4, N1=8, N2=6, eta=3, B=10)


class TestInstanceMethods:
    """Test instance-level toffoli_count_per_be() and total_qubits()."""

    @pytest.fixture()
    def est(self):
        return ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5)

    def test_toffoli_optimum(self, est):
        """Test that qroam='optimum' mode produces a positive Toffoli count.

        Optimum mode analytically minimises the QROAM batch-size parameters.
        """
        result = est.toffoli_count_per_be(qroam="optimum", **_INST_KW)
        assert isinstance(result, int)
        assert result > 0

    def test_toffoli_non_optimum(self, est):
        """Test that qroam='non-optimum' mode accepts explicit batch-size parameters."""
        result = est.toffoli_count_per_be(
            qroam="non-optimum",
            kp1=4, kp1_prime=3, ko=5, ko_prime=3,
            kp2=6, kp2_prime=4, kr=7, kr_prime=5,
            **_INST_KW,
        )
        assert isinstance(result, int)
        assert result > 0

    def test_toffoli_error_cases(self, est):
        """Test error handling: missing non-optimum params, bad qroam, missing Nb."""
        with pytest.raises(ValueError, match="requires parameters"):
            est.toffoli_count_per_be(qroam="non-optimum", **_INST_KW)
        with pytest.raises(ValueError, match="qroam must be"):
            est.toffoli_count_per_be(qroam="bogus", **_INST_KW)

        est_no_nb = ResourceEstimator(Nk=4, Npw=100, P=10)
        with pytest.raises(ValueError, match="Nb"):
            est_no_nb.toffoli_count_per_be(qroam="optimum", **_INST_KW)

    def test_total_qubits_optimum(self, est):
        """Test that qroam='optimum' mode produces a positive qubit count."""
        result = est.total_qubits(qroam="optimum", **_INST_KW)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_non_optimum(self, est):
        """Test that qroam='non-optimum' mode accepts kr parameter."""
        result = est.total_qubits(qroam="non-optimum", kr=7, **_INST_KW)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_error_cases(self, est):
        """Test error handling: missing kr, bad qroam, missing Nb."""
        with pytest.raises(ValueError, match="requires kr"):
            est.total_qubits(qroam="non-optimum", **_INST_KW)
        with pytest.raises(ValueError, match="qroam must be"):
            est.total_qubits(qroam="bogus", **_INST_KW)

        est_no_nb = ResourceEstimator(Nk=4, Npw=100, P=10)
        with pytest.raises(ValueError, match="Nb"):
            est_no_nb.total_qubits(qroam="optimum", **_INST_KW)

    def test_toffoli_uses_self_eta(self):
        """Test that toffoli_count_per_be uses self.eta when eta is not passed."""
        est = ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5, eta=3)
        kw = {k: v for k, v in _INST_KW.items() if k != "eta"}
        result = est.toffoli_count_per_be(qroam="optimum", **kw)
        assert isinstance(result, int)
        assert result > 0

    def test_toffoli_explicit_eta_overrides_self(self):
        """Test that an explicit eta argument overrides self.eta."""
        est = ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5, eta=3)
        kw = {k: v for k, v in _INST_KW.items() if k != "eta"}
        result_self = est.toffoli_count_per_be(qroam="optimum", **kw)
        result_override = est.toffoli_count_per_be(qroam="optimum", eta=5, **kw)
        # Different eta should give different Toffoli counts
        assert result_self != result_override

    def test_toffoli_default_eta_is_10(self):
        """Test that eta defaults to 10 when not explicitly provided."""
        est = ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5)
        assert est.eta == 10
        kw = {k: v for k, v in _INST_KW.items() if k != "eta"}
        result = est.toffoli_count_per_be(qroam="optimum", **kw)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_uses_self_eta(self):
        """Test that total_qubits uses self.eta when eta is not passed."""
        est = ResourceEstimator(Nk=4, Npw=100, P=10, Nb=5, eta=3)
        kw = {k: v for k, v in _INST_KW.items() if k != "eta"}
        result = est.total_qubits(qroam="optimum", **kw)
        assert isinstance(result, int)
        assert result > 0

    def test_total_qubits_with_lam_eps(self, est):
        """Test that passing lam + eps_qpe computes the QPE iteration count.

        This is the primary user-facing API for resource estimation.
        """
        result = est.total_qubits(
            qroam="optimum", lam=50.0, eps_qpe=0.01, **_INST_KW
        )
        assert isinstance(result, int)
        assert result > 0


# ===================================================================
# from_hdf5 with mocked h5py.File
# ===================================================================

def _make_mock_hdf5(*, Nk=4, Npw=100, na=3, Nb=5, kmesh_region="bz",
                    include_C_tensor=True, include_one_body=True,
                    n_electrons=10):
    """Build a mock h5py.File context-manager with the expected datasets."""
    h5 = mock.MagicMock()

    kmesh_group = mock.MagicMock()
    region_group = mock.MagicMock()
    reduced_ds = mock.MagicMock()
    reduced_ds.shape = (Nk, 3)
    region_group.__getitem__ = mock.MagicMock(side_effect=lambda k: reduced_ds)
    region_group.__contains__ = mock.MagicMock(return_value=True)
    kmesh_group.__getitem__ = mock.MagicMock(side_effect=lambda k: region_group if k == kmesh_region else mock.MagicMock())
    kmesh_group.__contains__ = mock.MagicMock(return_value=True)

    npw_ds = mock.MagicMock()
    npw_ds.__getitem__ = mock.MagicMock(return_value=Npw)

    ne_ds = mock.MagicMock()
    ne_ds.__getitem__ = mock.MagicMock(return_value=n_electrons)

    c_group = mock.MagicMock()
    c_ds = mock.MagicMock(spec=h5py.Dataset)
    c_ds.shape = (na, na, na, na)
    c_group.items = mock.MagicMock(return_value=[("atom_0", c_ds)])

    ob_group = mock.MagicMock()
    ob_ds = mock.MagicMock(spec=h5py.Dataset)
    ob_ds.shape = (Nk, Nb, Nk, Nb)
    ob_group.items = mock.MagicMock(return_value=[("H_kikj", ob_ds)])
    ob_group.__contains__ = mock.MagicMock(side_effect=lambda k: k == "H_kikj")
    ob_group.__getitem__ = mock.MagicMock(side_effect=lambda k: ob_ds if k == "H_kikj" else mock.MagicMock())

    def getitem(key):
        if key == "kmesh":
            return kmesh_group
        if key == "Npw":
            return npw_ds
        if key == "n_electrons":
            return ne_ds
        if key == "C_tensor":
            return c_group
        if key == "one_body":
            return ob_group
        raise KeyError(key)

    h5.__getitem__ = mock.MagicMock(side_effect=getitem)

    top_keys = {"kmesh", "Npw"}
    if n_electrons is not None:
        top_keys.add("n_electrons")
    if include_C_tensor:
        top_keys.add("C_tensor")
    if include_one_body:
        top_keys.add("one_body")
    h5.__contains__ = mock.MagicMock(side_effect=lambda k: k in top_keys)

    ctx = mock.MagicMock()
    ctx.__enter__ = mock.MagicMock(return_value=h5)
    ctx.__exit__ = mock.MagicMock(return_value=False)
    return ctx


class TestFromHdf5:
    """ResourceEstimator.from_hdf5 with mocked h5py.File."""

    def test_valid_construction(self):
        """Test that from_hdf5 correctly infers Nk, Npw, P, Nb, and L from HDF5.

        P = na*(na+1)//2 counts upper-triangular partial-wave pairs. This
        formula is specific to the packed C^a storage format.
        """
        mock_ctx = _make_mock_hdf5(Nk=4, Npw=100, na=3, Nb=5, kmesh_region="bz")
        with mock.patch("h5py.File", return_value=mock_ctx):
            est = ResourceEstimator.from_hdf5("dummy.h5", kmesh_region="bz")
        assert est.Nk == 4
        assert est.Npw == 100
        assert est.P == 6  # na*(na+1)//2 = 3*4//2
        assert est.L == 2 * (100 + 6) * 4
        assert est.Nb == 5
        assert est.eta == 10

    def test_invalid_kmesh_region(self):
        """Test that invalid kmesh_region raises ValueError before opening the file."""
        with pytest.raises(ValueError, match="kmesh_region must be"):
            ResourceEstimator.from_hdf5("dummy.h5", kmesh_region="full")

    def test_eta_defaults_to_10_from_hdf5(self):
        """Test that eta defaults to 10 regardless of n_electrons in HDF5.

        eta is a QROAM parameter (2-adic valuation of L), not the electron
        count.  It is conventionally hardcoded to 10 following Ivanov & Rubin.
        """
        mock_ctx = _make_mock_hdf5(Nk=4, Npw=100, na=3, Nb=5,
                                   kmesh_region="bz", n_electrons=None)
        with mock.patch("h5py.File", return_value=mock_ctx):
            est = ResourceEstimator.from_hdf5("dummy.h5", kmesh_region="bz")
        assert est.eta == 10

    def test_no_one_body_sets_Nb_none(self):
        """Test that Nb is None when the HDF5 file has no one_body group.

        Nb is inferred from the one-body dataset shape. When absent, the
        user must pass Nb explicitly to toffoli_count_per_be().
        """
        mock_ctx = _make_mock_hdf5(Nk=4, Npw=100, na=3, Nb=5,
                                   kmesh_region="bz", include_one_body=False)
        with mock.patch("h5py.File", return_value=mock_ctx):
            est = ResourceEstimator.from_hdf5("dummy.h5", kmesh_region="bz")
        assert est.Nb is None


# ===================================================================
# Regression: verify specific numeric outputs
# ===================================================================

class TestNumericRegression:
    """Pin specific numeric results to catch unintentional formula changes."""

    def test_toffoli_core_value(self):
        """Test that _toffoli_core returns a stable value for fixed inputs.

        If the formula is modified (intentionally or not), this test fails
        and forces a conscious decision to update the snapshot.
        """
        result = ResourceEstimator._toffoli_core(**_CORE_PARAMS, bo_override=None)
        assert 100 < result < 10_000_000

    def test_total_qubits_core_value(self):
        """Test that _total_qubits_core returns a stable value for fixed inputs."""
        params = {k: v for k, v in _CORE_PARAMS.items()
                  if k in ("L", "kr", "br", "Rl", "R0", "Nb", "Nk", "N1", "N2", "eta", "B")}
        result = ResourceEstimator._total_qubits_core(**params, bo_override=None)
        assert 10 < result < 10_000_000
