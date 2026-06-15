"""
Test that configuration validation runs early in the driver, before Hamiltonian loading.

This ensures fail-fast behavior when configuration is invalid.
"""

import pytest
import tempfile
import os
from unittest.mock import patch, MagicMock


def test_validation_called_before_hamiltonian_load():
    """Test that validation is called before get_physical_hamiltonian."""
    from qhat.analysis import driver

    # Track the order of calls
    call_order = []

    # Mock the functions we care about
    with patch('qhat.analysis.driver.load_configuration') as mock_load_config, \
         patch('qhat.analysis.driver.configure_logging') as mock_configure_logging, \
         patch('qhat.analysis.driver.validate_and_autocomplete_analysis_config') as mock_validate, \
         patch('qhat.analysis.driver.get_physical_hamiltonian') as mock_get_hamiltonian, \
         patch('qhat.analysis.driver.encode_as_unitary') as mock_encode_unitary, \
         patch('qhat.analysis.driver.build_algorithm') as mock_build_algorithm, \
         patch('qhat.analysis.driver.analyze_algorithm') as mock_analyze_algorithm:

        # Create a mock state object
        mock_state = MagicMock()
        mock_state.config_general.loglevel = 'info'
        mock_state.config_general.logfile = 'test.log'
        mock_state.config_general.git_hash = 'test123'
        mock_state.config_unitary.method = 'not ramped trotter'  # Skip Trotter computation
        mock_load_config.return_value = mock_state

        # Track call order
        def track_validate(*args, **kwargs):
            call_order.append('validate')

        def track_hamiltonian(*args, **kwargs):
            call_order.append('hamiltonian')
            return MagicMock()

        mock_validate.side_effect = track_validate
        mock_get_hamiltonian.side_effect = track_hamiltonian

        # Mock other functions to return reasonable values
        mock_encode_unitary.return_value = MagicMock()
        mock_build_algorithm.return_value = MagicMock()
        mock_analyze_algorithm.return_value = {}

        # Run the driver
        driver.run()

        # Verify validation was called before Hamiltonian loading
        assert call_order == ['validate', 'hamiltonian'], \
            f"Expected ['validate', 'hamiltonian'], got {call_order}"


def test_validation_error_prevents_hamiltonian_load():
    """Test that validation errors prevent expensive Hamiltonian loading."""
    from qhat.analysis import driver
    from qhat.analysis.config_types import AnalysisConfiguration

    # Track whether Hamiltonian was loaded
    hamiltonian_loaded = []

    with patch('qhat.analysis.driver.load_configuration') as mock_load_config, \
         patch('qhat.analysis.driver.configure_logging') as mock_configure_logging, \
         patch('qhat.analysis.driver.get_physical_hamiltonian') as mock_get_hamiltonian:

        # Create a mock state with invalid configuration
        mock_state = MagicMock()
        mock_state.config_general.loglevel = 'info'
        mock_state.config_general.logfile = 'test.log'
        mock_state.config_general.git_hash = 'test123'
        mock_state.config_unitary.method = 'not ramped trotter'  # Avoid Trotter code path

        # Create invalid config: eigenvalue errors without eigendecomposition
        invalid_config = AnalysisConfiguration()
        invalid_config.enable_eigenvalue_errors = True
        invalid_config.num_eigenvalues = 0  # Invalid!
        mock_state.config_analysis = invalid_config

        mock_load_config.return_value = mock_state

        # Track if Hamiltonian loading is attempted
        def track_hamiltonian(*args, **kwargs):
            hamiltonian_loaded.append(True)
            return MagicMock()

        mock_get_hamiltonian.side_effect = track_hamiltonian

        # Run should fail with ValueError during validation
        with pytest.raises(ValueError, match="enable_eigenvalue_errors requires eigendecomposition"):
            driver.run()

        # Hamiltonian should NOT have been loaded
        assert len(hamiltonian_loaded) == 0, \
            "Hamiltonian was loaded despite validation error (should fail fast)"


def test_validation_autocorrects_config():
    """Test that validation auto-corrects configuration as expected."""
    from qhat.analysis import driver
    from qhat.analysis.config_types import AnalysisConfiguration

    with patch('qhat.analysis.driver.load_configuration') as mock_load_config, \
         patch('qhat.analysis.driver.configure_logging') as mock_configure_logging, \
         patch('qhat.analysis.driver.get_physical_hamiltonian') as mock_get_hamiltonian, \
         patch('qhat.analysis.driver.encode_as_unitary') as mock_encode_unitary, \
         patch('qhat.analysis.driver.build_algorithm') as mock_build_algorithm, \
         patch('qhat.analysis.driver.analyze_algorithm') as mock_analyze_algorithm:

        # Create config that needs auto-correction
        config = AnalysisConfiguration()
        config.enable_eigenvalue_errors = True
        config.num_eigenvalues = 5
        config.eigendecomposition_matrices = 'approximate'  # Should be auto-corrected to 'both'

        mock_state = MagicMock()
        mock_state.config_general.loglevel = 'info'
        mock_state.config_general.logfile = 'test.log'
        mock_state.config_general.git_hash = 'test123'
        mock_state.config_analysis = config
        mock_state.config_unitary.method = 'not ramped trotter'
        mock_load_config.return_value = mock_state

        # Mock other functions
        mock_get_hamiltonian.return_value = MagicMock()
        mock_encode_unitary.return_value = MagicMock()
        mock_build_algorithm.return_value = MagicMock()
        mock_analyze_algorithm.return_value = {}

        # Run the driver
        driver.run()

        # Verify config was auto-corrected
        assert config.eigendecomposition_matrices == 'both', \
            f"Expected 'both', got '{config.eigendecomposition_matrices}'"
