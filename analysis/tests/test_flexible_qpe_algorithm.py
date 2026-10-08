"""
Tests for QPE construction in analysis/algorithm.py: the "QPE: QHAT flexible" and
"QPE: qualtran textbook" methods and phase-register sizing.
"""

import math

import pytest
from qualtran.bloqs.phase_estimation import TextbookQPE

from qhat.analysis.algorithm import build_algorithm, compute_initial_phase_qubits
from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
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


def phase_qubits(algorithm):
    """Phase-register size: TextbookQPE calls it m_bits, FlexibleQPE num_ancilla_qubits."""
    if isinstance(algorithm, TextbookQPE):
        return algorithm.m_bits
    return algorithm.num_ancilla_qubits


class TestRouting:

    def test_defaults(self, unitary):
        algorithm = build_algorithm(make_config(num_phase_qubits=3), unitary)
        assert algorithm == FlexibleQPE(unitary, 3)

    def test_explicit_components_case_insensitive(self, unitary):
        config = make_config(num_phase_qubits=3, ancilla_prep="Rectangular", qft_inv="TextBook")
        assert build_algorithm(config, unitary) == FlexibleQPE(unitary, 3)

    def test_method_case_insensitive(self, unitary):
        config = make_config(method="qpe: qhat FLEXIBLE", num_phase_qubits=2)
        assert isinstance(build_algorithm(config, unitary), FlexibleQPE)

    @pytest.mark.parametrize("key", ["ancilla_prep", "qft_inv"])
    def test_unknown_component_raises(self, unitary, key):
        config = make_config(num_phase_qubits=3, **{key: "nonsense"})
        with pytest.raises(ValueError, match=key):
            build_algorithm(config, unitary)

    def test_qualtran_textbook(self, unitary):
        config = make_config(method="QPE: qualtran textbook", num_phase_qubits=3)
        assert build_algorithm(config, unitary) == TextbookQPE(unitary, 3)


@pytest.mark.parametrize("method, qpe_class", [
    ("QPE: QHAT flexible", FlexibleQPE),
    ("QPE: qualtran textbook", TextbookQPE),
])
class TestPhaseQubits:
    """Both non-qubitized methods take the same config and give the same sizes."""

    def test_explicit_num_phase_qubits_is_used_as_given(self, trotter, method, qpe_class):
        config = make_config(method=method, num_phase_qubits=5, energy_error=1e-6,
                             probability_of_failure=1e-6)
        algorithm = build_algorithm(config, trotter)
        assert isinstance(algorithm, qpe_class)
        assert phase_qubits(algorithm) == 5

    def test_from_energy_error_and_probability_of_failure(self, trotter, method, qpe_class):
        dE, p_fail = 0.05, 0.1
        algorithm = build_algorithm(
            make_config(method=method, energy_error=dE, probability_of_failure=p_fail), trotter)
        assert isinstance(algorithm, qpe_class)
        phase_error = dE * trotter.time / (2 * math.pi)
        assert phase_qubits(algorithm) == textbook_qubits(phase_error, p_fail)

    def test_requires_probability_of_failure(self, trotter, method, qpe_class):
        with pytest.raises(ValueError, match="probability_of_failure"):
            build_algorithm(make_config(method=method, energy_error=0.05), trotter)

    @pytest.mark.parametrize("p_fail", [0.0, 1.0, 1.5, -0.1])
    def test_probability_of_failure_out_of_range_raises(self, trotter, method, qpe_class, p_fail):
        config = make_config(method=method, energy_error=0.05, probability_of_failure=p_fail)
        with pytest.raises(ValueError, match="probability_of_failure"):
            build_algorithm(config, trotter)

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
    assert phase_qubits(algorithm) == P0 + math.ceil(math.log2(2 + 1 / (2 * p_fail)))


def test_tensor_contract_shape(unitary):
    algorithm = build_algorithm(make_config(num_phase_qubits=2), unitary)
    assert algorithm.tensor_contract().shape == (2**4, 2**4)
