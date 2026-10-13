"""
Tests for QPE construction in analysis/algorithm.py: the "QPE: QHAT flexible" and
"QPE: qualtran textbook" methods and phase-register sizing.
"""

import functools
import itertools
import math

import numpy as np
import pytest
from qualtran.bloqs.phase_estimation import TextbookQPE

from qhat.analysis.algorithm import build_algorithm, qpe_energy_window
from qhat.analysis.config_types import (
    AlgorithmConfiguration, HamiltonianConfiguration, UnitaryConfiguration)
from qhat.analysis.file_io import load_state, save_state
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.matrix_eigendecomposition import compute_unitary_matrix
from qhat.analysis.numerical_simulation import numerical_simulation
from qhat.analysis.resource_estimation import estimate_resources
from qhat.analysis.unitary import encode_as_unitary
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.trotter_flattened import Trotterization

PAULIS = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]),
          "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1])}

# Commuting terms, so one Trotter step is exact. Eigenstates are the Bell states, with energies
# -1, 0, 0.4 and 0.6; the one-norm bound is [-1, 1].
BELL_TERMS = [("ZZ", 0.5), ("XX", 0.3), ("YY", 0.2)]


def hamiltonian_matrix(terms):
    return sum(c * functools.reduce(np.kron, [PAULIS[p] for p in s]) for s, c in terms)


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

    def test_default_method_is_valid(self, unitary):
        config = AlgorithmConfiguration()
        config.num_phase_qubits = 3
        assert build_algorithm(config, unitary) == TextbookQPE(unitary, 3)

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

    def test_explicit_num_phase_qubits_needs_no_conversion(self, unitary, method, qpe_class):
        """With num_phase_qubits set, energy_error is not converted (this unitary cannot)."""
        config = make_config(method=method, num_phase_qubits=3, energy_error=0.05)
        assert phase_qubits(build_algorithm(config, unitary)) == 3

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


class TestQPEEnergyWindow:

    def test_without_energy_error_is_unchanged(self):
        assert qpe_energy_window(make_config(), -1.5, 1.5) == (-1.5, 1.5)

    @pytest.mark.parametrize("Elo, Ehi, dE, lo, hi", [
        (-1.5, 1.5, 0.1, -1.6, 1.6),   # width 3 -> 2^5 * 0.1 = 3.2
        (0.0, 3.0, 0.1, -0.1, 3.1),    # off-center range: same width, same middle
        (-1.0, 1.0, 0.25, -1.0, 1.0),  # width 2 = 2^3 * 0.25 exactly: no widening
    ])
    def test_width_rounded_up_and_split_evenly(self, Elo, Ehi, dE, lo, hi):
        window = qpe_energy_window(make_config(energy_error=dE), Elo, Ehi)
        assert window == pytest.approx((lo, hi))


@pytest.mark.parametrize("method", ["QPE: QHAT flexible", "QPE: qualtran textbook"])
@pytest.mark.parametrize("W, dE", [
    (3.0, 0.1), (4.0, 0.5), (2.0, 0.25),
    (1.0, 0.165),  # round-off: phase error comes out as 0.12499999999999999, just below 2^-3
])
def test_driver_evolution_time_gives_P0_precision_bits(method, W, dE):
    """With t = 2 pi / (2^P0 dE) (the driver's formula with phase_scale_factor = 1), the phase
    error is 2^-P0, so the register is P0 plus the textbook confidence bits."""
    p_fail = 0.1
    config = make_config(method=method, energy_error=dE, probability_of_failure=p_fail)
    P0 = math.ceil(math.log2(W / dE))
    lo, hi = qpe_energy_window(config, -W / 2, W / 2)
    t = 2 * math.pi / (hi - lo)
    u = Trotterization.from_method(
        pauli_terms=[("XZ", 0.5), ("ZY", 0.3)], method="second order", time=t, num_steps=1)
    algorithm = build_algorithm(config, u)
    assert phase_qubits(algorithm) == P0 + math.ceil(math.log2(2 + 1 / (2 * p_fail)))


@pytest.mark.parametrize("Elo, Ehi", [(-1.0, 1.0), (-1.05, 0.65)], ids=["one-norm", "off-center"])
@pytest.mark.parametrize("dE, p_fail", [(0.1, 0.1), (0.15, 0.05)])
def test_energy_estimates_meet_requirements(Elo, Ehi, dE, p_fail):
    """End to end: window, t = 2 pi / (hi - lo), QPE sized from (dE, p_fail), decoded energies.

    For each eigenstate, Pr[|E_est - E| > dE] <= p_fail. U = exp(-i H t), so a phase-register
    value k decodes to the energy in [lo, hi) congruent to -(k / 2^m) (hi - lo).
    """
    config = make_config(energy_error=dE, probability_of_failure=p_fail)
    lo, hi = qpe_energy_window(config, Elo, Ehi)
    t = 2 * math.pi / (hi - lo)
    u = Trotterization.from_method(pauli_terms=BELL_TERMS, method="first order", time=t,
                                   num_steps=1)
    qpe = build_algorithm(config, u)
    m, n = qpe.num_ancilla_qubits, qpe.num_state_qubits
    energies, states = np.linalg.eigh(hamiltonian_matrix(BELL_TERMS))
    assert Elo <= energies.min() and energies.max() < Ehi  # guard: test Hamiltonian fits

    k = np.arange(2**m)
    estimates = lo + np.mod(-k / 2**m * (hi - lo) - lo, hi - lo)
    mat = qpe.tensor_contract()
    for E, psi in zip(energies, states.T):
        out = mat @ np.kron(np.eye(2**m)[0], psi)  # qpe_reg starts in |0>
        probs = np.sum(np.abs(out.reshape(2**m, 2**n))**2, axis=1)
        assert probs[np.abs(estimates - E) > dE].sum() <= p_fail


def diagonal_as_pauli_terms(diagonal):
    """I/Z Pauli strings and coefficients summing to diag(diagonal), with c_s = Tr(P_s H) / 2^n.

    The first character acts on the most significant bit. Example: diag(0, -1/2) = -0.25 I + 0.25 Z.
    """
    n = int(math.log2(len(diagonal)))
    terms = []
    for s in map("".join, itertools.product("IZ", repeat=n)):
        c = np.mean(np.diag(hamiltonian_matrix([(s, 1.0)])).real * diagonal)
        if abs(c) > 1e-12:
            terms.append((s, c))
    return terms


def qpe_from_pauli_file(path, time, num_trotter_steps, **algorithm):
    """Load a Pauli file, encode it with first-order Trotter, and build FlexibleQPE configured by
    the `algorithm` settings (e.g. num_phase_qubits=2).

    These are the driver's steps without its energy shift and choice of time, so the caller controls
    the eigenphases exactly.
    """
    config_hamiltonian = HamiltonianConfiguration()
    config_hamiltonian.load_pauli_strings(str(path))
    config_unitary = UnitaryConfiguration()
    config_unitary.encode_ramped_trotter(
        trotter_order="first order", trotter_steps=num_trotter_steps)
    unitary = encode_as_unitary(config_unitary, get_physical_hamiltonian(config_hamiltonian), time)
    return build_algorithm(make_config(**algorithm), unitary)


def write_diagonal_pauli_file(path, phases, energy_scale):
    """Write H = -energy_scale diag(phases) as a Pauli file. With t = 2 pi / energy_scale,
    U = exp(-i H t) has eigenphase phases[j] (in turns) on basis state j."""
    hamiltonian = -energy_scale * np.diag(phases)
    terms = diagonal_as_pauli_terms(np.diag(hamiltonian))
    np.testing.assert_allclose(hamiltonian_matrix(terms), hamiltonian)  # guard: conversion
    path.write_text("".join(f"{c:.17g} {s}\n" for s, c in terms))


@pytest.mark.parametrize("phases, target_state, expected_bits", [
    ([0, 1/2], [1, 0], "0"),
    ([0, 1/2], [0, 1], "1"),
    ([0, 1/4, 1/2, 3/4], [1, 0, 0, 0], "00"),
    ([0, 1/4, 1/2, 3/4], [0, 1, 0, 0], "01"),
    ([0, 1/4, 1/2, 3/4], [0, 0, 1, 0], "10"),
    ([0, 1/4, 1/2, 3/4], [0, 0, 0, 1], "11"),
])
@pytest.mark.parametrize("contraction", ["tensor_contract", "decomposition"])
def test_simulation_reads_exact_phases_from_pauli_file(
        tmp_path, phases, target_state, expected_bits, contraction):
    """U has eigenphases `phases` (in turns) on the computational basis states. With one phase
    qubit per target qubit these are exact register values, so QPE outputs |expected_bits>|target>.

    U = exp(-i H t) with t = 2 pi / W, so phase phi needs energy -W phi: H = -W diag(phases), e.g.
    phases [0, 1/2] -> H = diag(0, -0.85) = -0.425 I + 0.425 Z for W = 1.7.

    W is not 1 (t = 2 pi) because exp(-i 2 pi k P) = I for any integer k and Pauli string P, so a
    coefficient off by an integer would go unnoticed. Three Trotter steps (exact, as the terms
    commute) check that each step evolves for t / 3.
    """
    energy_scale = 1.7
    pauli_file = tmp_path / "hamiltonian.txt"
    write_diagonal_pauli_file(pauli_file, phases, energy_scale)

    m = len(expected_bits)
    qpe = qpe_from_pauli_file(pauli_file, time=2 * math.pi / energy_scale, num_trotter_steps=3,
                              num_phase_qubits=m)
    if contraction == "tensor_contract":
        mat = compute_unitary_matrix(qpe)  # what QHAT's numerical simulation uses
    else:
        mat = qpe.decompose_bloq().tensor_contract()  # the circuit, gate by gate

    phase_register_zero = np.eye(2**m)[0]  # qpe_reg comes first in the signature
    input_file = str(tmp_path / "input.npy")
    save_state(input_file, np.kron(phase_register_zero, target_state))
    result = numerical_simulation(input_file, qpe, mat)
    out = load_state(result["simulations"][0]["output_file"])

    probs = np.sum(np.abs(out.reshape(2**m, -1))**2, axis=1)
    assert probs[int(expected_bits, 2)] == pytest.approx(1.0, abs=1e-10)
    np.testing.assert_allclose(
        out, np.kron(np.eye(2**m)[int(expected_bits, 2)], target_state), atol=1e-10)


class TestConfigurationAndResources:

    def test_toml_saves_flexible_components(self):
        table = make_config(ancilla_prep="rectangular", qft_inv="textbook")._generate_TOML_table()
        assert table["ancilla_prep"] == "rectangular"
        assert table["qft_inv"] == "textbook"

    def test_toml_omits_unset_components(self):
        table = make_config()._generate_TOML_table()
        assert "ancilla_prep" not in table
        assert "qft_inv" not in table

    @pytest.mark.parametrize("estimator", ["pyLIQTR", "qualtran"])
    def test_resource_estimators_accept_flexible_qpe(self, trotter, estimator):
        algorithm = build_algorithm(make_config(num_phase_qubits=3), trotter)
        resources = estimate_resources(estimator, algorithm)[estimator.lower()]
        assert resources["qubit_count"] == 3 + 2
        assert resources["T_count"] > 0
