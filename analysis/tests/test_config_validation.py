"""
Tests for configuration validation and autocomplete functionality.

Tests the validate_and_autocomplete_analysis_config function which:
1. Validates configuration consistency
2. Auto-enables dependent analyses where appropriate
"""

import pytest
from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.analysis import validate_and_autocomplete_analysis_config


def test_eigenvalue_error_requires_eigendecomposition():
    """Test that eigenvalue errors auto-enable eigendecomposition computation."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    # Should NOT raise an error - eigendecompositions will be computed automatically
    validate_and_autocomplete_analysis_config(config)

    # Verify that eigendecompositions will be computed (checked by requires_* functions)
    from qhat.analysis.analysis import requires_exact_eigendecomposition, requires_approximate_eigendecomposition
    assert requires_exact_eigendecomposition(config), "Exact eigendecomposition should be required"
    assert requires_approximate_eigendecomposition(config), "Approximate eigendecomposition should be required"


def test_eigenvalue_error_autocorrects_eigendecomposition_matrices():
    """Test that eigenvalue errors work even without eigendecomposition_matrices set."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True

    # Should work without setting eigendecomposition_matrices
    validate_and_autocomplete_analysis_config(config)

    # Eigendecompositions will be computed automatically via requires_* functions
    from qhat.analysis.analysis import requires_exact_eigendecomposition, requires_approximate_eigendecomposition
    assert requires_exact_eigendecomposition(config)
    assert requires_approximate_eigendecomposition(config)


def test_eigenvalue_error_with_exact_eigendecomposition():
    """Test that eigenvalue errors work with any eigendecomposition setting."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.eigendecomposition_matrices = 'exact'

    # Should work - eigendecompositions will be computed as needed
    validate_and_autocomplete_analysis_config(config)

    # Both eigendecompositions will be computed because enable_eigenvalue_errors is True
    from qhat.analysis.analysis import requires_exact_eigendecomposition, requires_approximate_eigendecomposition
    assert requires_exact_eigendecomposition(config)
    assert requires_approximate_eigendecomposition(config)


def test_eigenvalue_error_already_set_to_both():
    """Test that eigenvalue errors work when eigendecomposition_matrices is already 'both'."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 5
    config.eigendecomposition_matrices = 'both'

    # Should pass without modification
    validate_and_autocomplete_analysis_config(config)

    assert config.eigendecomposition_matrices == 'both'


def test_eigenvalue_error_with_all_eigenvalues():
    """Test that eigenvalue errors work regardless of eigendecomposition_matrices setting."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 'all'
    config.eigendecomposition_matrices = 'approximate'

    # Should work - eigendecompositions computed automatically via requires_* functions
    validate_and_autocomplete_analysis_config(config)

    # Both eigendecompositions will be computed because enable_eigenvalue_errors is True
    from qhat.analysis.analysis import requires_exact_eigendecomposition, requires_approximate_eigendecomposition
    assert requires_exact_eigendecomposition(config)
    assert requires_approximate_eigendecomposition(config)


def test_matrix_norm_errors_pass_validation():
    """Test that matrix norm errors pass validation (matrices computed automatically)."""
    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'

    # Should pass - matrices will be computed automatically in analyze_algorithm
    validate_and_autocomplete_analysis_config(config)


def test_state_errors_pass_validation():
    """Test that state errors pass validation (matrices computed automatically)."""
    config = AnalysisConfiguration()
    config.error_state_inputs = 'state.npy'

    # Should pass - matrices will be computed automatically in analyze_algorithm
    validate_and_autocomplete_analysis_config(config)


def test_multiple_error_types_with_eigenvalue():
    """Test combining eigenvalue errors with other error types."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 5
    config.eigendecomposition_matrices = 'approximate'
    config.error_matrix_norms = 'frobenius'
    config.error_state_inputs = 'state.npy'

    # Should work - eigendecompositions computed automatically
    validate_and_autocomplete_analysis_config(config)

    # Both eigendecompositions will be computed because enable_eigenvalue_errors is True
    from qhat.analysis.analysis import requires_exact_eigendecomposition, requires_approximate_eigendecomposition
    assert requires_exact_eigendecomposition(config)
    assert requires_approximate_eigendecomposition(config)


def test_no_analyses_configured():
    """Test that validation passes with no analyses configured (will fail later in analyze_algorithm)."""
    config = AnalysisConfiguration()

    # Should pass - the "no analyses requested" check happens in analyze_algorithm
    validate_and_autocomplete_analysis_config(config)


def test_opportunistic_matrix_output_for_numerical_simulation():
    """Test that matrix output is auto-enabled when matrices will be computed."""
    config = AnalysisConfiguration()
    config.numerical_simulation_inputs = 'state.npy'

    # Before validation
    assert config.matrix_output_file is None

    # Validate - should auto-enable matrix output
    validate_and_autocomplete_analysis_config(config)

    # After validation
    assert config.matrix_output_file == 'unitary_matrix.npz'


def test_opportunistic_exact_matrix_output_for_error_analysis():
    """Test that exact matrix output is auto-enabled for error analyses."""
    config = AnalysisConfiguration()
    config.error_matrix_norms = 'frobenius'

    # Before validation
    assert config.exact_matrix_output_file is None
    assert config.matrix_output_file is None

    # Validate - should auto-enable both matrix outputs
    validate_and_autocomplete_analysis_config(config)

    # After validation
    assert config.exact_matrix_output_file == 'exact_hamiltonian.npz'
    assert config.matrix_output_file == 'unitary_matrix.npz'


def test_no_opportunistic_enabling_when_already_set():
    """Test that opportunistic enabling respects existing settings."""
    config = AnalysisConfiguration()
    config.numerical_simulation_inputs = 'state.npy'
    config.matrix_output_file = 'my_custom_name.npz'

    # Validate - should NOT change the custom filename
    validate_and_autocomplete_analysis_config(config)

    # Should keep custom name
    assert config.matrix_output_file == 'my_custom_name.npz'
