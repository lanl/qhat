"""
Tests for RampedTrotterizedUnitary (the original Trotterization implementation).
"""

import numpy as np
import pytest

from qhat.common.trotter_original import build_ramped_trotterized_unitary


class TestPower:
    """Test that U**k is the unitary U^k (needed for fast-forwarding in phase estimation)."""

    @pytest.mark.parametrize("method", ["first order", "second order"])
    @pytest.mark.parametrize("k", [1, 2, 3])
    def test_power_matches_matrix_power(self, method, k):
        trotter = build_ramped_trotterized_unitary(
            [("XZ", 0.5), ("ZY", 0.3)], method, timestep=0.3, numsteps=2)
        expected = np.linalg.matrix_power(trotter.tensor_contract(), k)
        np.testing.assert_allclose((trotter ** k).tensor_contract(), expected, atol=1e-10)

    def test_power_preserves_time_step(self):
        trotter = build_ramped_trotterized_unitary(
            [("X", 1.0), ("Z", 1.0)], "second order", timestep=0.5, numsteps=4)
        powered = trotter ** 3
        assert powered.numsteps == 12
        assert powered.timestep == pytest.approx(1.5)


class TestPhaseErrorFromEnergyError:
    """The eigenphase (in turns) of U changes by phase_error_from_energy_error(dE) per dE."""

    def test_matches_eigenphase_shift(self):
        # A single ZI term is exact (two qubits; this implementation fails on one-qubit terms)
        c, t = 0.3, 0.7
        trotter = build_ramped_trotterized_unitary([("ZI", c)], "first order", timestep=t,
                                                   numsteps=1)
        phases = np.angle(np.linalg.eigvals(trotter.tensor_contract())) / (2 * np.pi)
        assert np.ptp(phases) == pytest.approx(trotter.phase_error_from_energy_error(2 * c))
