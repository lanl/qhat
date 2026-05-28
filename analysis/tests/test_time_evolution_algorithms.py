"""
Tests for time evolution algorithm methods.

Tests both "time evolution" and "controlled time evolution" algorithm methods
by verifying:
1. Basic routing and method selection
2. Correctness via matrix comparison against analytical results
3. Resource estimation produces reasonable results
"""

import pytest
import numpy as np
import sys
import os

from qhat.analysis.algorithm import build_algorithm
from qhat.analysis.analysis import analyze_algorithm
from qhat.analysis.config_types import AlgorithmConfiguration, AnalysisConfiguration
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.pauli_utils import analytical_evolution


# =============================================================================
# Test: Basic Functionality and Routing
# =============================================================================

class TestAlgorithmRouting:
    """Test that algorithm methods are correctly routed."""

    def test_time_evolution_returns_unitary(self):
        """Time evolution method should return unitary unchanged."""
        unitary = PauliStringEvolution("XYZ", coefficient=1.0, time=0.5)
        config = AlgorithmConfiguration()
        config.method = "time evolution"

        algorithm = build_algorithm(config, unitary, P0=None)

        assert algorithm is unitary

    def test_controlled_time_evolution_adds_control(self):
        """Controlled time evolution should add exactly one control qubit."""
        unitary = PauliStringEvolution("XY", coefficient=1.0, time=1.0)
        config = AlgorithmConfiguration()
        config.method = "controlled time evolution"

        original_qubits = unitary.num_qubits
        algorithm = build_algorithm(config, unitary, P0=None)

        assert algorithm is not unitary
        assert algorithm.signature.n_qubits() == original_qubits + 1

    def test_case_insensitive_routing(self):
        """Method names should be case-insensitive."""
        unitary = PauliStringEvolution("Z", coefficient=1.0, time=1.0)

        for method in ["time evolution", "TIME EVOLUTION", "Time Evolution"]:
            config = AlgorithmConfiguration()
            config.method = method
            algorithm = build_algorithm(config, unitary, P0=None)
            assert algorithm is unitary

    def test_invalid_method_raises_error(self):
        """Invalid method names should raise ValueError."""
        unitary = PauliStringEvolution("X", coefficient=1.0, time=1.0)
        config = AlgorithmConfiguration()
        config.method = "invalid method"

        with pytest.raises(ValueError, match="Invalid algorithm method"):
            build_algorithm(config, unitary, P0=None)


# =============================================================================
# Test: Matrix Correctness
# =============================================================================

class TestMatrixCorrectness:
    """Verify algorithm matrices match analytical results."""

    def test_time_evolution_single_pauli(self):
        """Test time evolution produces correct matrix for single Pauli."""
        coef = 0.5
        time = 1.0
        unitary = PauliStringEvolution("X", coefficient=coef, time=time)

        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        # Get actual matrix
        U_actual = algorithm.tensor_contract()

        # Get expected matrix
        U_expected = analytical_evolution("X", coef, time)

        assert np.allclose(U_actual, U_expected), \
            "Time evolution matrix doesn't match analytical result"

    def test_time_evolution_two_qubit(self):
        """Test time evolution for two-qubit Pauli string."""
        coef = 0.3
        time = 2.0
        unitary = PauliStringEvolution("XY", coefficient=coef, time=time)

        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U_actual = algorithm.tensor_contract()
        U_expected = analytical_evolution("XY", coef, time)

        assert np.allclose(U_actual, U_expected)

    def test_time_evolution_three_qubit(self):
        """Test time evolution for three-qubit Pauli string."""
        coef = 0.2
        time = 1.5
        unitary = PauliStringEvolution("XYZ", coefficient=coef, time=time)

        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U_actual = algorithm.tensor_contract()
        U_expected = analytical_evolution("XYZ", coef, time)

        assert np.allclose(U_actual, U_expected)

    def test_controlled_time_evolution_structure(self):
        """Test controlled evolution has correct block structure."""
        coef = 0.5
        time = 1.0
        unitary = PauliStringEvolution("X", coefficient=coef, time=time)

        config = AlgorithmConfiguration()
        config.method = "controlled time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        # Get matrices
        U_controlled = algorithm.tensor_contract()
        U_base = analytical_evolution("X", coef, time)

        # Controlled-U should have structure:
        # [[I, 0],
        #  [0, U]]
        # where I and U are 2x2 blocks for single qubit
        dim = U_base.shape[0]

        # Check upper-left block is identity
        I_block = U_controlled[:dim, :dim]
        assert np.allclose(I_block, np.eye(dim)), \
            "Upper-left block should be identity"

        # Check upper-right is zero
        zero_block = U_controlled[:dim, dim:]
        assert np.allclose(zero_block, 0), \
            "Upper-right block should be zero"

        # Check lower-left is zero
        zero_block = U_controlled[dim:, :dim]
        assert np.allclose(zero_block, 0), \
            "Lower-left block should be zero"

        # Check lower-right is U
        U_block = U_controlled[dim:, dim:]
        assert np.allclose(U_block, U_base), \
            "Lower-right block should be the base unitary"

    def test_unitarity_preserved(self):
        """Verify both methods produce unitary matrices."""
        unitary = PauliStringEvolution("XZ", coefficient=0.7, time=1.2)

        # Test time evolution
        config_te = AlgorithmConfiguration()
        config_te.method = "time evolution"
        algorithm_te = build_algorithm(config_te, unitary, P0=None)
        U_te = algorithm_te.tensor_contract()
        U_te_dag_U = U_te.conj().T @ U_te
        assert np.allclose(U_te_dag_U, np.eye(U_te.shape[0])), \
            "Time evolution doesn't preserve unitarity"

        # Test controlled time evolution
        config_cte = AlgorithmConfiguration()
        config_cte.method = "controlled time evolution"
        algorithm_cte = build_algorithm(config_cte, unitary, P0=None)
        U_cte = algorithm_cte.tensor_contract()
        U_cte_dag_U = U_cte.conj().T @ U_cte
        assert np.allclose(U_cte_dag_U, np.eye(U_cte.shape[0])), \
            "Controlled time evolution doesn't preserve unitarity"

    def test_time_evolution_identity(self):
        """Test evolution under identity Pauli string."""
        coef = 1.0
        time = 1.0
        unitary = PauliStringEvolution("I", coefficient=coef, time=time)

        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U_actual = algorithm.tensor_contract()
        # exp(-i c t I) = exp(-i c t) * I (global phase)
        U_expected = np.exp(-1j * coef * time) * np.eye(2)

        assert np.allclose(U_actual, U_expected)


# =============================================================================
# Test: Resource Estimation
# =============================================================================

class TestResourceEstimation:
    """Test resource estimation on time evolution algorithms."""

    def test_time_evolution_resource_estimation(self):
        """Resource estimation should work on time evolution algorithms."""
        unitary = PauliStringEvolution("XY", coefficient=1.0, time=1.0)

        algorithm_config = AlgorithmConfiguration()
        algorithm_config.method = "time evolution"
        algorithm = build_algorithm(algorithm_config, unitary, P0=None)

        analysis_config = AnalysisConfiguration()
        analysis_config.resource_estimator = "pyLIQTR"

        results = analyze_algorithm(analysis_config, algorithm)

        assert "resource_estimates" in results
        resources = results["resource_estimates"]
        assert "T_count" in resources
        assert "Clifford_count" in resources

        # Should have non-negative counts
        assert resources["T_count"] >= 0
        assert resources["Clifford_count"] >= 0

    def test_controlled_time_evolution_resource_estimation(self):
        """Resource estimation should work on controlled time evolution algorithms."""
        unitary = PauliStringEvolution("XY", coefficient=1.0, time=1.0)

        algorithm_config = AlgorithmConfiguration()
        algorithm_config.method = "controlled time evolution"
        algorithm = build_algorithm(algorithm_config, unitary, P0=None)

        analysis_config = AnalysisConfiguration()
        analysis_config.resource_estimator = "pyLIQTR"

        results = analyze_algorithm(analysis_config, algorithm)

        assert "resource_estimates" in results
        resources = results["resource_estimates"]
        assert "T_count" in resources
        assert "Clifford_count" in resources
        assert resources["T_count"] >= 0
        assert resources["Clifford_count"] >= 0

    def test_controlled_has_higher_resources(self):
        """Controlled version should have at least as many resources as base."""
        unitary = PauliStringEvolution("XY", coefficient=1.0, time=1.0)
        analysis_config = AnalysisConfiguration()
        analysis_config.resource_estimator = "pyLIQTR"

        # Time evolution resources
        config_te = AlgorithmConfiguration()
        config_te.method = "time evolution"
        algorithm_te = build_algorithm(config_te, unitary, P0=None)
        results_te = analyze_algorithm(analysis_config, algorithm_te)
        resources_te = results_te["resource_estimates"]

        # Controlled time evolution resources
        config_cte = AlgorithmConfiguration()
        config_cte.method = "controlled time evolution"
        algorithm_cte = build_algorithm(config_cte, unitary, P0=None)
        results_cte = analyze_algorithm(analysis_config, algorithm_cte)
        resources_cte = results_cte["resource_estimates"]

        # Controlled should have at least as many gates
        assert resources_cte["T_count"] >= resources_te["T_count"]
        assert resources_cte["Clifford_count"] >= resources_te["Clifford_count"]

    def test_resource_counts_are_reasonable(self):
        """Resource counts should be finite and non-negative."""
        analysis_config = AnalysisConfiguration()
        analysis_config.resource_estimator = "pyLIQTR"
        algorithm_config = AlgorithmConfiguration()
        algorithm_config.method = "time evolution"

        # Test with a few different algorithm sizes
        for pauli_string in ["XY", "XYZ", "XYZI"]:
            unitary = PauliStringEvolution(pauli_string, coefficient=1.0, time=1.0)
            algorithm = build_algorithm(algorithm_config, unitary, P0=None)
            results = analyze_algorithm(analysis_config, algorithm)

            t_count = results["resource_estimates"]["T_count"]
            clifford_count = results["resource_estimates"]["Clifford_count"]

            # Should be non-negative and finite
            assert 0 <= t_count < 1e10, f"T-count {t_count} unreasonable for {pauli_string}"
            assert 0 <= clifford_count < 1e10, f"Clifford count unreasonable for {pauli_string}"


# =============================================================================
# Test: Edge Cases
# =============================================================================

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_zero_coefficient(self):
        """Zero coefficient should give identity evolution."""
        unitary = PauliStringEvolution("XY", coefficient=0.0, time=1.0)
        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U = algorithm.tensor_contract()
        # exp(-i * 0 * XY * t) = I
        assert np.allclose(U, np.eye(U.shape[0]))

    def test_zero_time(self):
        """Zero time should give identity evolution."""
        unitary = PauliStringEvolution("XY", coefficient=1.0, time=0.0)
        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U = algorithm.tensor_contract()
        assert np.allclose(U, np.eye(U.shape[0]))

    def test_large_time(self):
        """Test evolution with large time values."""
        coef = 0.5
        time = 10.0  # Large time
        unitary = PauliStringEvolution("Z", coefficient=coef, time=time)
        config = AlgorithmConfiguration()
        config.method = "time evolution"
        algorithm = build_algorithm(config, unitary, P0=None)

        U_actual = algorithm.tensor_contract()
        U_expected = analytical_evolution("Z", coef, time)

        assert np.allclose(U_actual, U_expected)


# =============================================================================
# Test Runner
# =============================================================================

if __name__ == "__main__":
    print("="*70)
    print("TIME EVOLUTION ALGORITHM TESTS")
    print("="*70)
    print()

    # Run with pytest
    pytest.main([__file__, "-v"])
