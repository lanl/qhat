"""
Tests for resource estimation with multiple methods.

Tests include:
- Single resource estimation method
- Multiple resource estimation methods
- Verification that all requested methods return results
"""

import pytest
from unittest.mock import Mock, patch
from qhat.analysis.resource_estimation import estimate_resources


# =============================================================================
# Mock Resource Estimation Functions
# =============================================================================

def mock_pyliqtr_resources(algorithm):
    """Mock pyLIQTR resource estimation."""
    return {
        "Clifford_count": 100,
        "T_count": 50,
        "qubit_count": 10
    }


def mock_qualtran_resources(algorithm):
    """Mock Qualtran resource estimation."""
    return {
        "Clifford_count": 105,
        "T_count": 52,
        "qubit_count": 10
    }


def mock_cirq_resources(algorithm):
    """Mock Cirq resource estimation (not implemented yet)."""
    raise NotImplementedError


# =============================================================================
# Resource Estimation Tests
# =============================================================================

def test_single_resource_estimator():
    """Test resource estimation with a single method."""
    mock_algorithm = Mock()

    with patch('qhat.analysis.resource_estimation.resource_estimation_pyliqtr',
               side_effect=mock_pyliqtr_resources):
        results = estimate_resources('pyliqtr', mock_algorithm)

    # Should return dict with one key
    assert isinstance(results, dict)
    assert 'pyliqtr' in results
    assert len(results) == 1

    # Check the nested structure
    assert results['pyliqtr']['Clifford_count'] == 100
    assert results['pyliqtr']['T_count'] == 50
    assert results['pyliqtr']['qubit_count'] == 10
    print("✓ Test passed: single_resource_estimator")


def test_multiple_resource_estimators():
    """Test resource estimation with multiple methods."""
    mock_algorithm = Mock()

    with patch('qhat.analysis.resource_estimation.resource_estimation_pyliqtr',
               side_effect=mock_pyliqtr_resources), \
         patch('qhat.analysis.resource_estimation.resource_estimation_qualtran',
               side_effect=mock_qualtran_resources):

        results = estimate_resources(['pyliqtr', 'qualtran'], mock_algorithm)

    # Should return dict with two keys
    assert isinstance(results, dict)
    assert 'pyliqtr' in results
    assert 'qualtran' in results
    assert len(results) == 2

    # Check pyliqtr results
    assert results['pyliqtr']['Clifford_count'] == 100
    assert results['pyliqtr']['T_count'] == 50
    assert results['pyliqtr']['qubit_count'] == 10

    # Check qualtran results
    assert results['qualtran']['Clifford_count'] == 105
    assert results['qualtran']['T_count'] == 52
    assert results['qualtran']['qubit_count'] == 10
    print("✓ Test passed: multiple_resource_estimators")


def test_case_insensitive_estimator_names():
    """Test that estimator names are case-insensitive."""
    mock_algorithm = Mock()

    with patch('qhat.analysis.resource_estimation.resource_estimation_pyliqtr',
               side_effect=mock_pyliqtr_resources):

        results1 = estimate_resources('PyLIQTR', mock_algorithm)
        results2 = estimate_resources('PYLIQTR', mock_algorithm)
        results3 = estimate_resources('pyliqtr', mock_algorithm)

    # All should have the same structure (normalized to lowercase key)
    assert 'pyliqtr' in results1
    assert 'pyliqtr' in results2
    assert 'pyliqtr' in results3
    print("✓ Test passed: case_insensitive_estimator_names")


def test_invalid_estimator_raises_error():
    """Test that invalid estimator name raises ValueError."""
    mock_algorithm = Mock()

    with pytest.raises(ValueError, match="Invalid resource estimator method"):
        estimate_resources('invalid_method', mock_algorithm)

    with pytest.raises(ValueError, match="Invalid resource estimator method"):
        estimate_resources(['invalid_method'], mock_algorithm)

    print("✓ Test passed: invalid_estimator_raises_error")


def test_string_input_converted_to_list():
    """Test that string input is treated the same as single-element list."""
    mock_algorithm = Mock()

    with patch('qhat.analysis.resource_estimation.resource_estimation_pyliqtr',
               side_effect=mock_pyliqtr_resources):

        results_string = estimate_resources('pyliqtr', mock_algorithm)
        results_list = estimate_resources(['pyliqtr'], mock_algorithm)

    # Both should produce the same structure
    assert results_string.keys() == results_list.keys()
    assert results_string['pyliqtr'] == results_list['pyliqtr']
    print("✓ Test passed: string_input_converted_to_list")


if __name__ == "__main__":
    test_single_resource_estimator()
    test_multiple_resource_estimators()
    test_case_insensitive_estimator_names()
    test_invalid_estimator_raises_error()
    test_string_input_converted_to_list()
    print("\n✓ All resource estimation tests passed!")
