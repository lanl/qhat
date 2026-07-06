"""
End-to-end integration tests for the prqpe estimator through the real QHAT seams.

Unlike test_prqpe_adapter.py, this module exercises the genuine pipeline objects
and therefore imports the full qhat stack (hamiltonian.py, analysis.py,
algorithm.py, unitary.py).  It is intended to be run separately, with the heavy
dependencies (pyLIQTR, qualtran, pyscf, ...) available.

A real Hamiltonian is built from data_pauli.json using the same loader pattern
as test_pauli_hamiltonian.py.  The test then walks the three integration seams:

- the unitary "none" no-op returns the Hamiltonian unchanged;
- the algorithm "QPE: partially randomized" method builds a PRQPEAlgorithm carrier;
- analyze_algorithm with resource_estimator="prqpe" produces a prqpe estimate.
"""

import os
import pytest

from qhat.analysis.config_types import (
    AlgorithmConfiguration,
    AnalysisConfiguration,
    GeneralConfiguration,
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    UnitaryConfiguration,
)
from qhat.analysis.hamiltonian import load_pauli
from qhat.analysis.unitary import encode_as_unitary
from qhat.analysis.algorithm import PRQPEAlgorithm, build_algorithm
from qhat.analysis.analysis import analyze_algorithm


# ==================================================================================
# Shared fixtures
# ==================================================================================

@pytest.fixture(scope="session")
def general_config():
    """Create a single GeneralConfiguration for the entire test session."""
    user_config = GeneralConfigurationUser()
    return GeneralConfiguration(user_config)


@pytest.fixture
def real_hamiltonian(general_config):
    """Load the real 4-qubit Hamiltonian from data_pauli.json."""
    test_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data_pauli.json",
    )
    config_ham = HamiltonianConfiguration()
    config_ham.load_pauli_strings(test_file)
    return load_pauli(config_ham)


# ==================================================================================
# Test: full prqpe integration through the real seams
# ==================================================================================

class TestPrqpeIntegration:
    """Exercise the unitary, algorithm, and analysis seams end to end."""

    def test_unitary_none_returns_hamiltonian_unchanged(self, real_hamiltonian):
        """The 'none' unitary method is a no-op: it returns the Hamiltonian."""
        config_unitary = UnitaryConfiguration()
        config_unitary.encode_none()
        assert config_unitary.method == "none"

        result = encode_as_unitary(config_unitary, real_hamiltonian, None)
        assert result is real_hamiltonian

    def test_build_algorithm_returns_prqpe_carrier(self, real_hamiltonian):
        """'QPE: partially randomized' builds a PRQPEAlgorithm carrier."""
        config_algorithm = AlgorithmConfiguration()
        config_algorithm.method = "QPE: partially randomized"
        config_algorithm.overlap = 1.0

        carrier = build_algorithm(config_algorithm, real_hamiltonian, None)

        assert isinstance(carrier, PRQPEAlgorithm)
        assert carrier.method == "partially_randomized"
        assert carrier.overlap == 1.0
        assert carrier.randomizer == "rte"

    def test_analyze_algorithm_produces_prqpe_estimate(self, real_hamiltonian):
        """analyze_algorithm dispatches to prqpe and yields a positive estimate."""
        config_algorithm = AlgorithmConfiguration()
        config_algorithm.method = "QPE: partially randomized"
        config_algorithm.overlap = 1.0
        carrier = build_algorithm(config_algorithm, real_hamiltonian, None)

        config_analysis = AnalysisConfiguration()
        config_analysis.resource_estimator = "prqpe"
        config_analysis.prqpe_target_precision = 0.0016
        # Auto C_gs: the 4-qubit system is well below prqpe_cgs_max_qubits.

        results = analyze_algorithm(
            config_analysis, carrier, hamiltonian=real_hamiltonian
        )

        assert "resource_estimates" in results
        estimate = results["resource_estimates"]
        assert estimate["estimator"] == "prqpe"
        assert estimate["toffoli_count"] > 0
        assert estimate["logical_qubits"] > 0
        assert estimate["C_gs"] > 0
