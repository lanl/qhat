"""
Tests for QPE construction in analysis/algorithm.py: the "QPE: QHAT flexible" and
"QPE: qualtran textbook" methods and phase-register sizing.
"""

import math

import pytest

from qhat.analysis.algorithm import (
    build_algorithm, compute_initial_phase_qubits, NewTextbookQPE)
from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.qpe_window_state import RectangularWindowState
from qhat.common.trotter_flattened import Trotterization


def make_config(method="QPE: QHAT flexible", **kwargs):
    config = AlgorithmConfiguration()
    config.method = method
    for key, value in kwargs.items():
        setattr(config, key, value)
    return config


@pytest.fixture
def unitary():
    return PauliStringEvolution("XY", coefficient=1.0, time=0.5)


@pytest.fixture
def trotter():
    return Trotterization.from_method(
        pauli_terms=[("XZ", 0.5), ("ZY", 0.3)], method="second order", time=0.4, num_steps=2)


def textbook_qubits(phase_error, p_fail):
    return math.ceil(math.log2(1 / phase_error)) + math.ceil(math.log2(2 + 1 / (2 * p_fail)))


class TestRouting:

    def test_defaults(self, unitary):
        algorithm = build_algorithm(make_config(num_phase_qubits=3), unitary)
        assert algorithm == FlexibleQPE(unitary, RectangularWindowState(3))

    def test_explicit_components_case_insensitive(self, unitary):
        config = make_config(num_phase_qubits=3, ctrl_state_prep="Rectangular", qft_inv="TextBook")
        assert build_algorithm(config, unitary) == FlexibleQPE(unitary, RectangularWindowState(3))

    def test_method_case_insensitive(self, unitary):
        config = make_config(method="qpe: qhat FLEXIBLE", num_phase_qubits=2)
        assert isinstance(build_algorithm(config, unitary), FlexibleQPE)

    @pytest.mark.parametrize("key", ["ctrl_state_prep", "qft_inv"])
    def test_unknown_component_raises(self, unitary, key):
        config = make_config(num_phase_qubits=3, **{key: "nonsense"})
        with pytest.raises(ValueError, match=key):
            build_algorithm(config, unitary)

    def test_qualtran_textbook(self, unitary):
        config = make_config(method="QPE: qualtran textbook", num_phase_qubits=3)
        algorithm = build_algorithm(config, unitary)
        assert isinstance(algorithm, NewTextbookQPE)
        assert algorithm.m_bits == 3


@pytest.mark.parametrize("method, qpe_class", [
    ("QPE: QHAT flexible", FlexibleQPE),
    ("QPE: qualtran textbook", NewTextbookQPE),
])
class TestPhaseQubits:
    """Both non-qubitized methods share the same sizing interface and results."""

    def test_explicit_num_phase_qubits_is_used_as_given(self, trotter, method, qpe_class):
        config = make_config(method=method, num_phase_qubits=5, energy_error=1e-6,
                             probability_of_failure=1e-6)
        algorithm = build_algorithm(config, trotter)
        assert isinstance(algorithm, qpe_class)
        assert algorithm.m_bits == 5

    def test_from_energy_error_and_probability_of_failure(self, trotter, method, qpe_class):
        dE, p_fail = 0.05, 0.1
        algorithm = build_algorithm(
            make_config(method=method, energy_error=dE, probability_of_failure=p_fail), trotter)
        phase_error = dE * trotter.time / (2 * math.pi)
        assert algorithm.m_bits == textbook_qubits(phase_error, p_fail)

    def test_requires_probability_of_failure(self, trotter, method, qpe_class):
        with pytest.raises(ValueError, match="probability_of_failure"):
            build_algorithm(make_config(method=method, energy_error=0.05), trotter)

    def test_requires_energy_error_or_num_phase_qubits(self, trotter, method, qpe_class):
        with pytest.raises(ValueError, match="num_phase_qubits"):
            build_algorithm(make_config(method=method, probability_of_failure=0.1), trotter)

    def test_unitary_without_conversion_requires_num_phase_qubits(self, unitary, method,
                                                                    qpe_class):
        config = make_config(method=method, energy_error=0.05, probability_of_failure=0.1)
        with pytest.raises(ValueError, match="phase_error_from_energy_error"):
            build_algorithm(config, unitary)


@pytest.mark.parametrize("method", ["QPE: QHAT flexible", "QPE: qualtran textbook"])
@pytest.mark.parametrize("W, dE", [(3.0, 0.1), (4.0, 0.5), (2.0, 0.25)])
def test_driver_evolution_time_gives_P0_precision_bits(method, W, dE):
    """With the driver's t = 2 pi / (2^P0 dE), the phase error is 2^-P0, so the register is
    P0 plus the textbook confidence bits."""
    p_fail = 0.1
    config = make_config(method=method, energy_error=dE, probability_of_failure=p_fail)
    P0, Elo, Ehi = compute_initial_phase_qubits(config, -W / 2, W / 2)
    t = 2 * math.pi / (Ehi - Elo)
    u = Trotterization.from_method(
        pauli_terms=[("XZ", 0.5), ("ZY", 0.3)], method="second order", time=t, num_steps=1)
    algorithm = build_algorithm(config, u)
    assert algorithm.m_bits == P0 + math.ceil(math.log2(2 + 1 / (2 * p_fail)))


def test_tensor_contract_shape(unitary):
    algorithm = build_algorithm(make_config(num_phase_qubits=2), unitary)
    assert algorithm.tensor_contract().shape == (2**4, 2**4)
