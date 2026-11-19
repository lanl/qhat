"""Complex Hermitian Hamiltonians (e.g. relativistic chemistry) must produce real quantities.

Regression tests for https://github.com/lanl/qhat/issues/115, where complex-typed integrals made
the energy bounds, evolution time, and QPE phase error complex.
"""

import math
import warnings

import numpy as np
import pytest
from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.analysis.algorithm import qpe_energy_window, qpe_requirements
from qhat.analysis.config_types import (
    AlgorithmConfiguration,
    HamiltonianConfiguration,
    UnitaryConfiguration,
)
from qhat.analysis.hamiltonian import load_numpy
from qhat.analysis.unitary import encode_double_factorization, encode_ramped_trotter
from qhat.common.commuting_pauli_string_evolution import CommutingPauliStringEvolution
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.qpe_window_state import precision_bits
from .term_helpers import single_pauli_terms
from qhat.common.trotter_flattened import Trotterization

NORB = 3


def complex_hermitian_integrals(seed=0):
    rng = np.random.default_rng(seed)
    h1 = rng.normal(size=(NORB, NORB)) + 0.3j * rng.normal(size=(NORB, NORB))
    h1 = (h1 + h1.conj().T) / 2
    h2 = 0.1 * (rng.normal(size=(NORB,) * 4) + 1j * rng.normal(size=(NORB,) * 4))
    h2 = (h2 + h2.transpose(3, 2, 1, 0).conj()) / 2
    return h1, h2


def load(tmp_path, transform="JW", **arrays):
    filename = tmp_path / "hamiltonian.npz"
    np.savez(filename, **arrays)
    config = HamiltonianConfiguration()
    config.load_second_quantization(str(filename), fermion_to_qubit_transform=transform)
    return config, load_numpy(config)


@pytest.fixture
def complex_hamiltonian(tmp_path):
    h1, h2 = complex_hermitian_integrals()
    return load(tmp_path, ecore=np.array(2.5), h1e=h1, eri=h2)


class TestLoading:

    @pytest.mark.parametrize("transform", ["JW", "BK"])
    def test_pauli_coefficients_are_floats(self, tmp_path, transform):
        h1, h2 = complex_hermitian_integrals()
        _, H = load(tmp_path, transform, ecore=np.array(2.5), h1e=h1, eri=h2)
        terms = H.to_pauli_strings()
        assert len(terms) > 1
        assert all(type(c) is float for c in terms.values())
        assert all(type(term.coefficient) is float for term in H.terms())

    def test_energy_bounds_are_floats(self, complex_hamiltonian):
        config, H = complex_hamiltonian
        lo, hi = H.compute_initial_energy_bounds(config)
        assert type(lo) is float and type(hi) is float
        assert lo < hi

    def test_rejects_non_hermitian_one_body(self, tmp_path):
        h1, h2 = complex_hermitian_integrals()
        h1[0, 1] += 0.2j
        with pytest.raises(ValueError, match="one-body tensor is not Hermitian"):
            load(tmp_path, ecore=np.array(0.0), h1e=h1, eri=h2)

    def test_rejects_non_hermitian_two_body(self, tmp_path):
        h1, h2 = complex_hermitian_integrals()
        h2[0, 1, 2, 0] += 0.2j
        with pytest.raises(ValueError, match="two-body tensor is not Hermitian"):
            load(tmp_path, ecore=np.array(0.0), h1e=h1, eri=h2)

    def test_rejects_complex_constant(self, tmp_path):
        h1, h2 = complex_hermitian_integrals()
        with pytest.raises(ValueError, match="constant energy term must be real"):
            load(tmp_path, ecore=np.array(1.0 + 0.5j), h1e=h1, eri=h2)

    def test_accepts_complex_typed_constant_with_zero_imaginary_part(self, tmp_path):
        h1, h2 = complex_hermitian_integrals()
        _, H = load(tmp_path, ecore=np.array(1.0 + 0j), h1e=h1, eri=h2)
        assert H.core_operator().constant == 1.0

    def test_rejects_complex_energy_shift(self, complex_hamiltonian):
        _, H = complex_hamiltonian
        with pytest.raises(ValueError, match="energy shift must be real"):
            H.shift_energy(1.0 + 0.1j)


class TestNoComplexWarnings:
    """Run the ramped-Trotter + QPE sizing path with ComplexWarning turned into an error."""

    def test_trotter_qpe_path_is_real(self, complex_hamiltonian):
        config_ham, H = complex_hamiltonian
        energy_error = 0.05
        with warnings.catch_warnings():
            warnings.simplefilter("error", np.exceptions.ComplexWarning)

            lo, hi = H.compute_initial_energy_bounds(config_ham)
            E0 = (lo + hi) / 2
            H.shift_energy(-E0)
            config_alg = AlgorithmConfiguration()
            config_alg.method = "QPE: qhat flexible"
            config_alg.energy_error = energy_error
            config_alg.probability_of_failure = 0.1
            lo, hi = qpe_energy_window(config_alg, lo - E0, hi - E0)
            tevol_hbar = 2 * math.pi / (hi - lo)
            assert type(tevol_hbar) is float

            config_unitary = UnitaryConfiguration()
            config_unitary.encode_ramped_trotter(
                energy_error=energy_error, ordering_method="lexicographical",
                trotter_order="second order", trotter_steps=2)
            trotter = encode_ramped_trotter(config_unitary, H, tevol_hbar)
            assert type(trotter.time) is float
            assert all(type(term.coefficient) is float for term in trotter.hamiltonian_terms)
            trotter.decompose_bloq()

            phase_error, _ = qpe_requirements(config_alg, trotter)
            assert type(phase_error) is float
            assert precision_bits(phase_error) > 0


class TestSanityChecks:

    def test_precision_bits_rejects_complex(self):
        with pytest.raises(ValueError, match="phase_error must be real"):
            precision_bits(0.01 + 0.001j)

    def test_precision_bits_accepts_complex_typed_real(self):
        assert precision_bits(np.complex128(0.125 + 0j)) == 3

    def test_energy_window_rejects_complex_bounds(self):
        config = AlgorithmConfiguration()
        config.energy_error = 0.1
        with pytest.raises(ValueError, match="lower energy bound must be real"):
            qpe_energy_window(config, -1.0 + 0.5j, 1.0)

    def test_pauli_string_evolution_rejects_complex(self):
        with pytest.raises(ValueError, match="coefficient must be real"):
            PauliStringEvolution(pauli_string="XZ", coefficient=0.5 + 0.1j, time=1.0)
        with pytest.raises(ValueError, match="time must be real"):
            PauliStringEvolution(pauli_string="XZ", coefficient=0.5, time=1.0 + 0.1j)

    def test_trotterization_rejects_complex(self):
        terms = single_pauli_terms([("XI", 0.5), ("ZZ", 0.25)])
        with pytest.raises(ValueError, match="Trotterization time must be real"):
            Trotterization.from_method(terms, "second order", time=1.0 + 0.1j, num_steps=2)
        with pytest.raises(ValueError, match="coefficient of ZZ must be real"):
            single_pauli_terms([("XI", 0.5), ("ZZ", 0.25j)])

    def test_pauli_string_evolution_stores_floats(self):
        pse = PauliStringEvolution(pauli_string="XZ", coefficient=np.complex128(0.37 + 1e-19j),
                                   time=np.complex128(0.81 + 1e-20j), hbar=np.float64(1.0))
        assert type(pse.coefficient) is float and pse.coefficient == 0.37
        assert type(pse.time) is float and pse.time == 0.81
        assert type(pse.hbar) is float
        t_complexity(pse)  # builds the Cirq gate, which fails on complex exponents

    def test_commuting_pauli_string_evolution_stores_floats(self):
        cpse = CommutingPauliStringEvolution(
            pauli_terms=(("XZ", np.complex128(0.37 + 1e-19j)), ("ZX", -0.21)),
            time=np.complex128(0.81 + 1e-20j))
        assert all(type(c) is float for _, c in cpse.pauli_terms)
        assert cpse.pauli_terms[0] == ("XZ", 0.37)
        assert type(cpse.time) is float and type(cpse.hbar) is float
        t_complexity(cpse)

    def test_trotterization_stores_floats(self):
        trotter = Trotterization.from_method(
            single_pauli_terms([("XI", np.complex128(0.37 + 1e-19j)), ("ZZ", -0.21)]),
            "second order", time=np.complex128(0.81 + 1e-20j), num_steps=3)
        assert all(type(term.coefficient) is float for term in trotter.hamiltonian_terms)
        assert trotter.hamiltonian_terms[0].dense_string() == "XI"
        assert trotter.hamiltonian_terms[0].coefficient == 0.37
        assert all(type(c) is float for c in trotter.coefficients)
        assert type(trotter.time) is float and type(trotter.hbar) is float
        t_complexity(trotter)

    def test_ramped_trotter_accepts_complex_typed_energy_error(self, complex_hamiltonian):
        _, H = complex_hamiltonian
        config_unitary = UnitaryConfiguration()
        config_unitary.encode_ramped_trotter(
            energy_error=np.complex128(0.05 + 0j), trotter_order="second order")
        config_unitary.error_coeff_mode = "exact"
        with warnings.catch_warnings():
            warnings.simplefilter("error", np.exceptions.ComplexWarning)
            trotter = encode_ramped_trotter(config_unitary, H, 0.37)
        assert trotter.num_steps >= 1

    def test_ramped_trotter_rejects_complex_time(self, complex_hamiltonian):
        _, H = complex_hamiltonian
        config_unitary = UnitaryConfiguration()
        config_unitary.encode_ramped_trotter(
            energy_error=0.05, trotter_order="second order", trotter_steps=1)
        with pytest.raises(ValueError, match="evolution time t/hbar must be real"):
            encode_ramped_trotter(config_unitary, H, 0.5 + 0.1j)


class TestDoubleFactorizationGuard:

    def test_rejects_complex_integrals(self, complex_hamiltonian):
        _, H = complex_hamiltonian
        config_unitary = UnitaryConfiguration()
        config_unitary.encode_double_factorization(energy_error=1e-3)
        with pytest.raises(NotImplementedError, match="does not support complex integrals"):
            encode_double_factorization(config_unitary, H)
