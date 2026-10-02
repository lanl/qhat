"""
Tests for the "QPE: qhat textbook" algorithm method.
"""

import math

import pytest

from qhat.analysis.algorithm import build_algorithm
from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.qpe_window_state import RectangularWindowState


def make_config(**kwargs):
    config = AlgorithmConfiguration()
    config.method = "QPE: qhat textbook"
    for key, value in kwargs.items():
        setattr(config, key, value)
    return config


@pytest.fixture
def unitary():
    return PauliStringEvolution("XY", coefficient=1.0, time=0.5)


class TestRouting:

    def test_defaults(self, unitary):
        algorithm = build_algorithm(make_config(num_phase_qubits=3), unitary, P0=None)
        assert isinstance(algorithm, FlexibleQPE)
        assert algorithm.unitary is unitary
        assert algorithm.ctrl_state_prep == RectangularWindowState(3)
        assert algorithm.qft_inv == FlexibleQPE(unitary, RectangularWindowState(3)).qft_inv

    def test_explicit_components_case_insensitive(self, unitary):
        config = make_config(num_phase_qubits=3, ctrl_state_prep="Rectangular", qft_inv="TextBook")
        algorithm = build_algorithm(config, unitary, P0=None)
        assert algorithm == FlexibleQPE(unitary, RectangularWindowState(3))

    def test_method_case_insensitive(self, unitary):
        config = make_config(num_phase_qubits=2)
        config.method = "qpe: QHAT Textbook"
        assert isinstance(build_algorithm(config, unitary, P0=None), FlexibleQPE)

    @pytest.mark.parametrize("key", ["ctrl_state_prep", "qft_inv"])
    def test_unknown_component_raises(self, unitary, key):
        config = make_config(num_phase_qubits=3, **{key: "nonsense"})
        with pytest.raises(ValueError, match=key):
            build_algorithm(config, unitary, P0=None)


class TestPhaseQubits:

    def test_from_num_phase_qubits(self, unitary):
        algorithm = build_algorithm(make_config(num_phase_qubits=5), unitary, P0=None)
        assert algorithm.m_bits == 5

    def test_from_P0_and_probability_of_failure(self, unitary):
        p_fail = 0.1
        algorithm = build_algorithm(make_config(probability_of_failure=p_fail), unitary, P0=4)
        assert algorithm.m_bits == 4 + math.ceil(math.log2(2.0 + 0.5 / p_fail))


def test_tensor_contract_shape(unitary):
    algorithm = build_algorithm(make_config(num_phase_qubits=2), unitary, P0=None)
    assert algorithm.tensor_contract().shape == (2**4, 2**4)
