from __future__ import annotations

import sys
import types
from unittest import mock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Mock GPAW modules so the file can be imported without GPAW installed.
# We must mock every gpaw.* import that extractor.py uses at module level
# *before* importing it.
# ---------------------------------------------------------------------------

_GPAW_MOCKS: dict[str, types.ModuleType] = {}


def _ensure_gpaw_mocked():
    """Insert lightweight mock modules for every gpaw sub-import.

    Must register explicit attributes for every ``from gpaw.X import Y``
    that extractor.py uses at module scope.
    """
    if _GPAW_MOCKS:
        return  # already done

    names = [
        "gpaw",
        "gpaw.arraydict",
        "gpaw.transformers",
        "gpaw.utilities",
        "gpaw.matrix",
    ]
    for name in names:
        mod = types.ModuleType(name)
        mod.__path__ = []  # make each a package so sub-imports work
        _GPAW_MOCKS[name] = mod
        sys.modules.setdefault(name, mod)

    # Explicit attributes needed by ``from gpaw import GPAW`` etc.
    sys.modules["gpaw"].GPAW = mock.MagicMock()  # type: ignore[attr-defined]
    sys.modules["gpaw.arraydict"].ArrayDict = mock.MagicMock()  # type: ignore[attr-defined]
    sys.modules["gpaw.transformers"].Transformer = mock.MagicMock()  # type: ignore[attr-defined]
    sys.modules["gpaw.matrix"].matrix_matrix_multiply = mock.MagicMock()  # type: ignore[attr-defined]


_ensure_gpaw_mocked()

# Now we can safely import the module under test
from bloch_paw.extractor import pack_hermitian, PawExtractor  # noqa: E402


# ===================================================================
# pack_hermitian
# ===================================================================

class TestPackHermitian:
    """Tests for the standalone pack_hermitian() helper."""

    def test_3x3_hermitian(self):
        """Test that upper-triangular elements are extracted in row-major order.

        This is the core serialisation primitive for PAW on-site tensors —
        if the packing order is wrong, all downstream C^a / D^a contractions
        will silently produce garbage.
        """
        H = np.array([
            [1.0 + 0j,  2.0 + 1j,  3.0 - 2j],
            [2.0 - 1j,  4.0 + 0j,  5.0 + 3j],
            [3.0 + 2j,  5.0 - 3j,  6.0 + 0j],
        ])
        packed = pack_hermitian(H)
        assert packed.shape == (6,)
        expected = np.array([1.0+0j, 2.0+1j, 3.0-2j, 4.0+0j, 5.0+3j, 6.0+0j])
        np.testing.assert_array_equal(packed, expected)

    def test_identity_diagonal_positions(self):
        """Test that packing I_4 puts 1s at diagonal positions and 0s elsewhere.

        Verifies the index mapping for a larger matrix where diagonal vs
        off-diagonal positions are easy to enumerate.
        """
        I4 = np.eye(4, dtype=complex)
        packed = pack_hermitian(I4)
        assert packed.shape == (10,)
        diag_positions = [0, 4, 7, 9]
        for i, val in enumerate(packed):
            if i in diag_positions:
                assert val == 1.0 + 0j
            else:
                assert val == 0.0 + 0j

    def test_output_length_scales_correctly(self):
        """Test that packed length equals n(n+1)/2 for multiple sizes.

        The n(n+1)/2 formula is assumed throughout the C^a tensor
        decomposition — if it's wrong the atom-dict writer silently
        truncates or pads data.
        """
        for n in [2, 5, 8, 16]:
            H = np.eye(n, dtype=complex)
            packed = pack_hermitian(H)
            assert len(packed) == n * (n + 1) // 2


# ===================================================================
# PawExtractor construction with mocks
# ===================================================================

def _make_mock_calc(*, Nk=4, nbands=5, grid_shape=(10, 10, 10)):
    """Create a mock GPAW calculator with just enough attributes for __init__."""
    calc = mock.MagicMock()

    # wfs
    wfs = mock.MagicMock()
    wfs.gd = mock.MagicMock()
    wfs.gd.N_c = grid_shape
    calc.wfs = wfs

    # wfs.kd
    kd = mock.MagicMock()
    wfs.kd = kd

    # density.ghat
    density = mock.MagicMock()
    ghat = mock.MagicMock()
    pd_ghat = mock.MagicMock()
    pd_ghat.gd = mock.MagicMock()
    pd_ghat.gd.N_c = grid_shape
    ghat.pd = pd_ghat
    density.ghat = ghat
    calc.density = density

    # atoms.cell
    cell_array = np.eye(3) * 5.0
    atoms = mock.MagicMock()
    atoms.cell.reciprocal.return_value = np.linalg.inv(cell_array).T
    atoms.cell.array = cell_array
    calc.atoms = atoms

    return calc


class TestPawExtractorConstruction:
    """Tests for PawExtractor.__init__ with mocked GPAW calculator."""

    def test_defaults(self):
        """Test that default construction sets zero thresholds, auto k-norm, and no nbands.

        These defaults matter because non-zero thresholds silently zero out
        small matrix elements during extraction. Verifying all five threshold
        fields plus k_norm_mode and nbands prevents accidental default changes.
        """
        calc = _make_mock_calc()
        ext = PawExtractor(calc)
        assert ext.calc is calc
        assert ext.wfs is calc.wfs
        assert ext.default_nbands is None
        assert ext._thr_rho == 0.0
        assert ext._thr_D == 0.0
        assert ext._thr_C == 0.0
        assert ext._thr_h == 0.0
        assert ext._thr_kappa == 0.0
        assert ext._k_norm_mode == "auto"

    def test_custom_parameters(self):
        """Test that non-default thresholds, nbands, and k_norm_mode are stored.

        Users pass these when they want to sparsify the extracted tensors for
        performance. Verifies the constructor doesn't silently ignore arguments.
        """
        calc = _make_mock_calc()
        ext = PawExtractor(calc, nbands=5, thr_rho=1e-8, thr_D=1e-6, k_norm_mode="on")
        assert ext.default_nbands == 5
        assert ext._thr_rho == 1e-8
        assert ext._thr_D == 1e-6
        assert ext._k_norm_mode == "on"
