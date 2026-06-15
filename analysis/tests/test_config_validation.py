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
    """Test that eigenvalue errors require eigendecomposition to be configured."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    # num_eigenvalues is 0 by default

    with pytest.raises(ValueError, match="enable_eigenvalue_errors requires eigendecomposition"):
        validate_and_autocomplete_analysis_config(config)


def test_eigenvalue_error_autocorrects_eigendecomposition_matrices():
    """Test that eigenvalue errors auto-set eigendecomposition_matrices to 'both'."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 5
    config.eigendecomposition_matrices = 'approximate'

    # Should auto-correct to 'both'
    validate_and_autocomplete_analysis_config(config)

    assert config.eigendecomposition_matrices == 'both'


def test_eigenvalue_error_with_exact_eigendecomposition():
    """Test that eigenvalue errors auto-set eigendecomposition_matrices from 'exact' to 'both'."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 5
    config.eigendecomposition_matrices = 'exact'

    # Should auto-correct to 'both'
    validate_and_autocomplete_analysis_config(config)

    assert config.eigendecomposition_matrices == 'both'


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
    """Test that eigenvalue errors work with num_eigenvalues='all'."""
    config = AnalysisConfiguration()
    config.enable_eigenvalue_errors = True
    config.num_eigenvalues = 'all'
    config.eigendecomposition_matrices = 'approximate'

    # Should auto-correct to 'both'
    validate_and_autocomplete_analysis_config(config)

    assert config.eigendecomposition_matrices == 'both'


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

    # Should auto-correct eigendecomposition_matrices to 'both'
    validate_and_autocomplete_analysis_config(config)

    assert config.eigendecomposition_matrices == 'both'


def test_no_analyses_configured():
    """Test that validation passes with no analyses configured (will fail later in analyze_algorithm)."""
    config = AnalysisConfiguration()

    # Should pass - the "no analyses requested" check happens in analyze_algorithm
    validate_and_autocomplete_analysis_config(config)
