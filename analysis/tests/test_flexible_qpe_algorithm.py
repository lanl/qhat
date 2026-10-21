"""
Tests for QPE construction in analysis/algorithm.py: the "QPE: QHAT flexible" and
"QPE: qualtran textbook" methods and phase-register sizing.
"""

import functools
import itertools
import logging
import math

import numpy as np
import pytest
import tomlkit
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
from qhat.common.qpe_window_state import precision_bits
from qhat.common.tapered_qpe_sizing import fewest_runs, per_run_qpe_failure, size_tapered_qpe
from qhat.common.tapered_window_state import (
    BandwidthOptimizedWindowState, DPSSWindowState, TruncatedDPSSWindowState)
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


TAPERED = {
    "dpss": DPSSWindowState,
    "truncated dpss": TruncatedDPSSWindowState,
    "bandwidth-optimized": BandwidthOptimizedWindowState,
}


class TestTaperedWindows:

    @pytest.mark.parametrize("name, cls", [
        ("DPSS", DPSSWindowState),
        ("Truncated DPSS", TruncatedDPSSWindowState),
        ("Bandwidth-Optimized", BandwidthOptimizedWindowState),
    ])
    def test_routing_case_insensitive(self, unitary, name, cls):
        config = make_config(ancilla_prep=name, num_precision_qubits=2, num_confidence_qubits=3)
        assert build_algorithm(config, unitary) == FlexibleQPE(unitary, 5, cls(2, 3))

    @pytest.mark.parametrize("name, cls", TAPERED.items())
    @pytest.mark.parametrize("counts", [
        dict(num_phase_qubits=5, num_confidence_qubits=3),
        dict(num_phase_qubits=5, num_precision_qubits=2),
        dict(num_phase_qubits=5, num_precision_qubits=2, num_confidence_qubits=3),
    ], ids=["phase, confidence", "phase, precision", "all three"])
    def test_num_phase_qubits_may_replace_one_count(self, unitary, name, cls, counts):
        config = make_config(ancilla_prep=name, **counts)
        assert build_algorithm(config, unitary).ancilla_prep == cls(2, 3)

    def test_counts_take_precedence_over_requirements(self, trotter):
        config = make_config(ancilla_prep="dpss", num_precision_qubits=2, num_confidence_qubits=3,
                             energy_error=1e-6, probability_of_failure=1e-6)
        assert build_algorithm(config, trotter).ancilla_prep == DPSSWindowState(2, 3)

    @pytest.mark.parametrize("name", TAPERED)
    def test_num_phase_qubits_alone_raises(self, unitary, name):
        config = make_config(ancilla_prep=name, num_phase_qubits=5)
        with pytest.raises(ValueError, match=f'"{name}" with num_phase_qubits = 5: .*'
                           "cannot be sized from num_phase_qubits alone.*energy_error"):
            build_algorithm(config, unitary)

    @pytest.mark.parametrize("name", ["rectangular", *TAPERED])
    @pytest.mark.parametrize("counts, match", [
        (dict(num_precision_qubits=3), r"Cannot size .*\['num_precision_qubits'\]"),
        (dict(num_confidence_qubits=3), r"Cannot size .*\['num_confidence_qubits'\]"),
        (dict(num_phase_qubits=5, num_precision_qubits=3, num_confidence_qubits=3),
         "must equal"),
    ], ids=["precision only", "confidence only", "inconsistent"])
    def test_invalid_counts_raise(self, unitary, name, counts, match):
        with pytest.raises(ValueError, match=match):
            build_algorithm(make_config(ancilla_prep=name, **counts), unitary)

    def test_rectangular_accepts_a_split(self, unitary):
        config = make_config(num_precision_qubits=2, num_confidence_qubits=3)
        assert build_algorithm(config, unitary) == FlexibleQPE(unitary, 5)

    @pytest.mark.parametrize("name, cls", TAPERED.items())
    @pytest.mark.parametrize("overlap", [None, 0.3, 1.0])
    def test_from_energy_error_and_probability_of_failure(self, trotter, name, cls, overlap):
        dE, p_fail = 0.05, 1e-6
        config = make_config(ancilla_prep=name, energy_error=dE, probability_of_failure=p_fail,
                             overlap=overlap)
        phase_error = dE * trotter.time / (2 * math.pi)
        qpe = build_algorithm(config, trotter)
        plan = size_tapered_qpe(cls, phase_error, p_fail, 0.5 if overlap is None else overlap)
        assert qpe == FlexibleQPE(
            trotter, plan.num_phase_qubits,
            cls(plan.num_precision_qubits, plan.num_confidence_qubits))
        # One fewer precision bit than the rectangular window.
        assert qpe.ancilla_prep.num_precision_qubits == precision_bits(phase_error) - 1
        if overlap == 1.0:  # a single run: the single-run rule, with far fewer confidence bits
            assert qpe.ancilla_prep == cls.from_requirements(phase_error, p_fail)
            assert qpe.num_ancilla_qubits < textbook_qubits(phase_error, p_fail)


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


@pytest.mark.parametrize("ancilla_prep", ["dpss", "truncated dpss", "bandwidth-optimized"])
def test_tapered_qpe_meets_requirements_through_pauli_file(tmp_path, ancilla_prep):
    """End to end with a tapered window sized from (energy_error, probability_of_failure) for an
    initial state that is the target eigenstate (overlap 1), so a single run.

    U has off-grid eigenphases on the computational basis states (see
    `test_simulation_reads_exact_phases_from_pauli_file`). Here energy_error 0.2 is a phase error
    of 0.2 / 1.7 = 0.118 turns, so 3 precision qubits (delta = 1/16); p_fail = 1e-3 gives 4
    confidence qubits (DPSS) or 5.
    For each eigenstate, Pr[estimate more than 0.118 turns from the phase] <= 1e-3.
    """
    phases = [0.1234, 0.377, 0.6021, 0.871]
    energy_scale, dE, p_fail = 1.7, 0.2, 1e-3
    pauli_file = tmp_path / "hamiltonian.txt"
    write_diagonal_pauli_file(pauli_file, phases, energy_scale)
    settings = dict(time=2 * math.pi / energy_scale, num_trotter_steps=3)

    qpe = qpe_from_pauli_file(pauli_file, ancilla_prep=ancilla_prep, energy_error=dE,
                              probability_of_failure=p_fail, overlap=1.0, **settings)
    num_phase = qpe.num_ancilla_qubits
    assert qpe.ancilla_prep.num_precision_qubits == 3 and num_phase in (7, 8)
    # Guard: a rectangular window of the same size misses the target, so the taper matters.
    rectangular = qpe_from_pauli_file(pauli_file, num_phase_qubits=num_phase, **settings)

    phase_error = dE / energy_scale
    k = np.arange(2**num_phase)
    for qpe_, meets_target in [(qpe, True), (rectangular, False)]:
        mat = compute_unitary_matrix(qpe_)
        failures = []
        for j, phase in enumerate(phases):
            input_file = str(tmp_path / f"input_{j}.npy")
            save_state(input_file, np.kron(np.eye(2**num_phase)[0], np.eye(4)[j]))
            result = numerical_simulation(input_file, qpe_, mat)
            out = load_state(result["simulations"][0]["output_file"])
            probs = np.sum(np.abs(out.reshape(2**num_phase, -1))**2, axis=1)
            distance = np.abs(k / 2**num_phase - phase)
            distance = np.minimum(distance, 1 - distance)
            failures.append(probs[distance > phase_error].sum())
        assert (max(failures) <= p_fail) == meets_target


# The diagonal test problem: H = -1.7 diag(PHASES), t = 2 pi / 1.7, so energies are -1.7 x PHASES
# (in [-1.7, 0)) and the ground state is basis state 3, at -1.48. energy_error 0.2 is a phase error
# of 0.118 turns.
PHASES = [0.1234, 0.377, 0.6021, 0.871]
ENERGY_SCALE = 1.7
SETTINGS = dict(time=2 * math.pi / ENERGY_SCALE, num_trotter_steps=3)


@pytest.fixture
def pauli_file(tmp_path):
    path = tmp_path / "hamiltonian.txt"
    write_diagonal_pauli_file(path, PHASES, ENERGY_SCALE)
    return path


def direct_multi_run_failure(n, k, p_under, p_over):
    """P(the k-th smallest of n estimates is an under- or overestimate), as binomial sums: at
    least k underestimates, or at least n - k + 1 overestimates."""
    def at_least(j_min, q):
        return sum(math.comb(n, j) * q**j * (1 - q)**(n - j) for j in range(j_min, n + 1))
    return at_least(k, p_under) + at_least(n - k + 1, p_over)


@pytest.mark.parametrize("ancilla_prep, cls", TAPERED.items())
def test_tapered_qpe_multi_run_plan_through_pauli_file(tmp_path, pauli_file, ancilla_prep, cls):
    """End to end at the default overlap 0.5: the circuit is sized for N runs, reporting the K-th
    smallest energy estimate, with energy_error 0.2 and p_fail 1e-3 (3 precision and 4 confidence
    qubits; 14 runs for the DPSS, 15 otherwise; K = 2).

    The initial state has weight 0.5 on the ground state and 0.2, 0.3 on two excited states. From
    the simulated distribution of one run's energy estimate, the N-run failure is computed
    directly; it must be within the sizing's bound and p_fail.
    """
    dE, p_fail = 0.2, 1e-3
    qpe = qpe_from_pauli_file(pauli_file, ancilla_prep=ancilla_prep, energy_error=dE,
                              probability_of_failure=p_fail, **SETTINGS)
    plan = size_tapered_qpe(cls, dE / ENERGY_SCALE, p_fail, 0.5)
    assert qpe.ancilla_prep == cls(3, 4)
    assert (plan.num_precision_qubits, plan.num_confidence_qubits) == (3, 4)
    assert plan.num_runs == (14 if cls is DPSSWindowState else 15) and plan.quantile == 2

    m = qpe.num_ancilla_qubits
    input_file = str(tmp_path / "input.npy")
    save_state(input_file, np.kron(np.eye(2**m)[0], np.sqrt([0.2, 0.3, 0.0, 0.5])))
    result = numerical_simulation(input_file, qpe, compute_unitary_matrix(qpe))
    out = load_state(result["simulations"][0]["output_file"])
    probs = np.sum(np.abs(out.reshape(2**m, -1))**2, axis=1)

    # U = exp(-i H t): register value k decodes to the energy in [-1.7, 0) congruent to
    # -1.7 k / 2^m.
    estimates = -ENERGY_SCALE + np.mod(-ENERGY_SCALE * np.arange(2**m) / 2**m, ENERGY_SCALE)
    ground = -ENERGY_SCALE * PHASES[3]
    p_under = probs[estimates < ground - dE].sum()
    p_over = probs[estimates > ground + dE].sum()
    u = per_run_qpe_failure(cls, 3, 4)
    assert p_under <= u and 1 - p_under - p_over >= 0.5 * (1 - u)
    failure = direct_multi_run_failure(plan.num_runs, plan.quantile, p_under, p_over)
    assert failure <= plan.multi_run_failure <= p_fail
    # Guard: a single run fails about half the time.
    assert p_over > 0.45


@pytest.fixture
def qhat_log(caplog):
    """caplog for qhat's info-level records (the "qhat" logger does not propagate to root)."""
    qhat_logger = logging.getLogger("qhat")
    level = qhat_logger.level
    qhat_logger.addHandler(caplog.handler)
    qhat_logger.setLevel(logging.INFO)
    caplog.handler.setLevel(logging.INFO)
    yield caplog
    qhat_logger.removeHandler(caplog.handler)
    qhat_logger.setLevel(level)


def algorithm_messages(caplog, level=logging.INFO):
    return [r.getMessage() for r in caplog.records
            if r.name == "qhat.analysis.algorithm" and r.levelno == level]


TABLE_HEADER = "--     overlap   probability of failure"
NOTE = "-- NOTE: the analyses below report on a single run of the algorithm."


class TestMultiRunSizing:
    """The tapered path of "QPE: QHAT flexible" with the multi-run sizing, through a Pauli file."""

    def test_default_log_block(self, pauli_file, qhat_log):
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", energy_error=0.2,
                                  probability_of_failure=1e-3, **SETTINGS)
        assert qpe.ancilla_prep == DPSSWindowState(3, 4)
        info = algorithm_messages(qhat_log)
        start = info.index("-- overlap with target eigenstate = 0.5 (default)")
        table = [f"--     {p:<9.2f} <= {f}" for p, f in zip(
            [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
            ["9.92e-01", "5.85e-01", "1.98e-01", "4.75e-02", "8.11e-03", "9.19e-04", "6.56e-05",
             "8.40e-06", "6.80e-06", "6.79e-06", "6.79e-06"])]
        assert info[start:] == [
            "-- overlap with target eigenstate = 0.5 (default)",
            "-- report the K-th smallest of N runs: K chosen by sizing (default)",
            "-- tapered QPE sizing: 3 precision + 4 confidence qubits = 7 phase qubits",
            "-- if you run the circuit 1 time:",
            "--     probability of failure <= 5.00e-01, energy error <= 0.2",
            "-- if you run the circuit 14 times and take the 2nd-smallest estimate:",
            "--     probability of failure <= 9.19e-04, energy error <= 0.2",
            "-- if the true overlap differs from 0.5 (same 14 runs, 2nd-smallest estimate; "
            "energy error <= 0.2 in every row):",
            TABLE_HEADER,
            *table,
            NOTE,
        ]
        assert not algorithm_messages(qhat_log, logging.WARNING)

    def test_user_specified_options(self, pauli_file, qhat_log):
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="truncated dpss", energy_error=0.2,
                                  probability_of_failure=1e-3, overlap=0.9,
                                  force_lowest_estimate=True, **SETTINGS)
        plan = size_tapered_qpe(TruncatedDPSSWindowState, 0.2 / ENERGY_SCALE, 1e-3, 0.9, True)
        assert plan.quantile == 1
        assert qpe.ancilla_prep == TruncatedDPSSWindowState(plan.num_precision_qubits,
                                                            plan.num_confidence_qubits)
        info = algorithm_messages(qhat_log)
        assert "-- overlap with target eigenstate = 0.9 (user-specified)" in info
        assert ("-- report the smallest of N runs (K = 1): forced by "
                "algorithm.force_lowest_estimate (user-specified)") in info
        assert (f"-- if you run the circuit {plan.num_runs} times and take the smallest "
                "estimate:") in info

    def test_force_false_is_user_specified(self, pauli_file, qhat_log):
        qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", energy_error=0.2,
                            probability_of_failure=1e-3, force_lowest_estimate=False, **SETTINGS)
        assert ("-- report the K-th smallest of N runs: K chosen by sizing (user-specified)"
                in algorithm_messages(qhat_log))

    @pytest.mark.parametrize("overlap, error", [
        (0.0, ValueError), (-0.2, ValueError), (1.5, ValueError), (math.nan, ValueError),
        (0.5 + 0.1j, ValueError), ("0.5", TypeError), (True, TypeError),
    ])
    def test_invalid_overlap_raises(self, pauli_file, overlap, error):
        with pytest.raises(error, match="overlap"):
            qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", energy_error=0.2,
                                probability_of_failure=1e-3, overlap=overlap, **SETTINGS)

    @pytest.mark.parametrize("force", ["yes", 1, None])
    def test_invalid_force_lowest_estimate_raises(self, pauli_file, force):
        config = dict(ancilla_prep="dpss", energy_error=0.2, probability_of_failure=1e-3,
                      force_lowest_estimate=force)
        if force is None:  # unset is fine
            qpe_from_pauli_file(pauli_file, **config, **SETTINGS)
            return
        with pytest.raises(TypeError, match="force_lowest_estimate"):
            qpe_from_pauli_file(pauli_file, **config, **SETTINGS)

    def test_overlap_one_is_the_single_run_rule(self, pauli_file, qhat_log):
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="bandwidth-optimized",
                                  energy_error=0.2, probability_of_failure=1e-3, overlap=1.0,
                                  **SETTINGS)
        assert qpe.ancilla_prep == BandwidthOptimizedWindowState.from_requirements(
            0.2 / ENERGY_SCALE, 1e-3)
        assert "-- if you run the circuit 1 time and take the smallest estimate:" in \
            algorithm_messages(qhat_log)

    def test_fallback_when_no_plan(self, pauli_file, qhat_log):
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", energy_error=0.2,
                                  probability_of_failure=1e-3, overlap=1e-17, **SETTINGS)
        assert qpe.ancilla_prep == DPSSWindowState.from_requirements(0.2 / ENERGY_SCALE, 1e-3)
        errors = algorithm_messages(qhat_log, logging.ERROR)
        assert len(errors) == 1 and "FALLBACK" in errors[0] and "1e-17" in errors[0]
        info = algorithm_messages(qhat_log)
        assert ("-- FALLBACK tapered QPE sizing: 3 precision + 4 confidence qubits = 7 phase "
                "qubits") in info
        assert TABLE_HEADER not in info and info[-1] == NOTE

    @pytest.mark.parametrize("probability_of_failure", [1e-3, None])
    def test_given_counts(self, pauli_file, qhat_log, probability_of_failure):
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", num_precision_qubits=3,
                                  num_confidence_qubits=5, overlap=0.4,
                                  probability_of_failure=probability_of_failure, **SETTINGS)
        assert qpe == qpe_from_pauli_file(pauli_file, ancilla_prep="dpss",
                                          num_phase_qubits=8, num_confidence_qubits=5, **SETTINGS)
        info = algorithm_messages(qhat_log)
        assert ("-- tapered QPE: user-specified 3 precision + 5 confidence qubits = 8 phase "
                "qubits") in info
        assert info[-1] == NOTE
        warnings = algorithm_messages(qhat_log, logging.WARNING)
        if probability_of_failure is None:
            assert TABLE_HEADER not in info
            assert any("skipping the number of runs" in message for message in info)
            assert len(warnings) == 1 and "algorithm.overlap is ignored" in warnings[0]
            return
        plan = fewest_runs(DPSSWindowState, 3, 5, probability_of_failure, 0.4)
        # The bound comes from the precision qubits: 2^-4 turns is 1.7 / 16 = 0.106 in energy.
        assert (f"--     probability of failure <= {plan.multi_run_failure:.2e}, energy error "
                "<= 0.106") in info
        assert (f"-- if you run the circuit {plan.num_runs} times and take the "
                f"{'smallest' if plan.quantile == 1 else '2nd-smallest'} estimate:") in info
        assert plan.quantile in (1, 2)
        assert TABLE_HEADER in info and not warnings

    @pytest.mark.parametrize("energy_error, warns", [(0.2, False), (0.1, True)])
    def test_given_counts_warn_if_energy_error_unmet(self, pauli_file, qhat_log, energy_error,
                                                     warns):
        # 3 precision qubits bound the energy error by 0.106.
        qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", num_precision_qubits=3,
                            num_confidence_qubits=4, energy_error=energy_error,
                            probability_of_failure=1e-3, **SETTINGS)
        warnings = algorithm_messages(qhat_log, logging.WARNING)
        assert len(warnings) == warns
        if warns:
            assert "0.106" in warnings[0] and "energy_error = 0.1" in warnings[0]

    def test_given_counts_with_no_plan(self, pauli_file, qhat_log):
        # One confidence qubit: no bound on a run's failure, so no number of runs works.
        qpe = qpe_from_pauli_file(pauli_file, ancilla_prep="dpss", num_precision_qubits=3,
                                  num_confidence_qubits=1, probability_of_failure=1e-3,
                                  **SETTINGS)
        assert qpe.ancilla_prep == DPSSWindowState(3, 1)
        errors = algorithm_messages(qhat_log, logging.ERROR)
        assert len(errors) == 1 and "No number of runs" in errors[0]
        info = algorithm_messages(qhat_log)
        assert TABLE_HEADER not in info and info[-1] == NOTE

    @pytest.mark.parametrize("config", [
        dict(ancilla_prep="rectangular", energy_error=0.2, probability_of_failure=1e-3),
        dict(num_phase_qubits=4),
        dict(method="QPE: qualtran textbook", num_phase_qubits=4),
    ], ids=["rectangular from requirements", "rectangular from counts", "textbook"])
    def test_unused_options_warn(self, pauli_file, qhat_log, config):
        qpe_from_pauli_file(pauli_file, overlap=0.3, force_lowest_estimate=True, **config,
                            **SETTINGS)
        warnings = algorithm_messages(qhat_log, logging.WARNING)
        assert len(warnings) == 2
        assert "algorithm.overlap is ignored" in warnings[0]
        assert "algorithm.force_lowest_estimate is ignored" in warnings[1]
        assert NOTE not in algorithm_messages(qhat_log)


class TestConfigurationAndResources:

    def test_toml_saves_multi_run_options(self):
        config = make_config(ancilla_prep="dpss", overlap=0.25, force_lowest_estimate=True)
        table = tomlkit.parse(tomlkit.dumps({"algorithm": config._generate_TOML_table()}))
        assert table["algorithm"]["overlap"] == 0.25
        assert table["algorithm"]["force_lowest_estimate"] is True
        table = make_config(ancilla_prep="dpss")._generate_TOML_table()
        assert "overlap" not in table and "force_lowest_estimate" not in table

    def test_toml_saves_flexible_components(self):
        table = make_config(ancilla_prep="rectangular", qft_inv="textbook")._generate_TOML_table()
        assert table["ancilla_prep"] == "rectangular"
        assert table["qft_inv"] == "textbook"

    def test_toml_omits_unset_components(self):
        table = make_config()._generate_TOML_table()
        for key in ("ancilla_prep", "qft_inv", "num_precision_qubits", "num_confidence_qubits"):
            assert key not in table

    def test_toml_saves_qubit_split(self):
        config = make_config(ancilla_prep="dpss", num_precision_qubits=3, num_confidence_qubits=2)
        table = config._generate_TOML_table()
        assert (table["num_precision_qubits"], table["num_confidence_qubits"]) == (3, 2)

    @pytest.mark.parametrize("estimator", ["pyLIQTR", "qualtran"])
    def test_resource_estimators_accept_flexible_qpe(self, trotter, estimator):
        algorithm = build_algorithm(make_config(num_phase_qubits=3), trotter)
        resources = estimate_resources(estimator, algorithm)[estimator.lower()]
        assert resources["qubit_count"] == 3 + 2
        assert resources["T_count"] > 0

    @pytest.mark.parametrize("name", TAPERED)
    def test_pyliqtr_counts_the_window_state(self, trotter, name):
        """A tapered QPE costs the rectangular one plus its window circuit (the Hadamards of the
        rectangular window cost no T gates)."""
        def t_count(bloq):
            return estimate_resources("pyLIQTR", bloq)["pyliqtr"]["T_count"]
        config = make_config(ancilla_prep=name, num_precision_qubits=2, num_confidence_qubits=3)
        tapered = build_algorithm(config, trotter)
        rectangular = build_algorithm(make_config(num_phase_qubits=5), trotter)
        assert t_count(tapered.ancilla_prep) > 0
        assert t_count(tapered) == t_count(rectangular) + t_count(tapered.ancilla_prep)
        assert estimate_resources("pyLIQTR", tapered)["pyliqtr"]["qubit_count"] == 5 + 2
