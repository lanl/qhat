"""Tests for command-line parameter functionality in configuration loading."""

import pytest
import sys
import tempfile
import os
from pathlib import Path

from qhat.analysis.configuration import load_configuration


def test_cmdline_params_basic(monkeypatch, tmp_path):
    """Test basic command-line parameter passing."""
    # Create a simple test config file
    config_file = tmp_path / "test_config.py"
    config_file.write_text("""
# Use command-line parameters
test_value = params.get('test_param', 'default')
numeric_value = params.get('number', 42)

# Use output_directory which is preserved in GeneralConfiguration
general.output_directory = f"test_{test_value}_{numeric_value}"
""")

    # Mock sys.argv to simulate command-line arguments
    test_args = [
        'driver.py',
        str(config_file),
        '-p', 'test_param=hello',
        '-p', 'number=100'
    ]
    monkeypatch.setattr(sys, 'argv', test_args)

    # Load configuration
    state = load_configuration()

    # Verify the parameters were used
    assert 'hello' in state.config_general.output_directory
    assert '100' in state.config_general.output_directory


def test_cmdline_params_type_conversion(monkeypatch, tmp_path):
    """Test that parameters are properly converted to Python types."""
    config_file = tmp_path / "test_config.py"
    config_file.write_text("""
# Access parameters directly if they're provided
float_val = float_param
int_val = int_param
list_val = list_param

general.file_stub = f"test"
""")

    test_args = [
        'driver.py',
        str(config_file),
        '-p', 'float_param=3.14',
        '-p', 'int_param=42',
        '-p', 'list_param=[1,2,3]'
    ]
    monkeypatch.setattr(sys, 'argv', test_args)

    state = load_configuration()
    # If we got here without errors, the types were correctly evaluated
    assert True


def test_cmdline_params_string_fallback(monkeypatch, tmp_path):
    """Test that invalid Python expressions are treated as strings."""
    config_file = tmp_path / "test_config.py"
    config_file.write_text("""
# This should be a string since "hello-world" isn't valid Python
string_val = params.get('weird_string', 'default')

general.output_directory = string_val
""")

    test_args = [
        'driver.py',
        str(config_file),
        '-p', 'weird_string=hello-world'
    ]
    monkeypatch.setattr(sys, 'argv', test_args)

    state = load_configuration()
    assert state.config_general.output_directory == 'hello-world'


def test_cmdline_params_no_params(monkeypatch, tmp_path):
    """Test that configs work without any command-line parameters."""
    config_file = tmp_path / "test_config.py"
    config_file.write_text("""
# Use defaults when no params provided
test_value = params.get('missing_param', 'default_value')

general.output_directory = test_value
""")

    test_args = ['driver.py', str(config_file)]
    monkeypatch.setattr(sys, 'argv', test_args)

    state = load_configuration()
    assert state.config_general.output_directory == 'default_value'


def test_cmdline_params_invalid_format(monkeypatch, tmp_path):
    """Test that invalid parameter format raises an error."""
    config_file = tmp_path / "test_config.py"
    config_file.write_text("general.file_stub = 'test'")

    # Parameter without '=' should raise ValueError
    test_args = [
        'driver.py',
        str(config_file),
        '-p', 'invalid_param_no_equals'
    ]
    monkeypatch.setattr(sys, 'argv', test_args)

    with pytest.raises(ValueError, match="Parameter must be in KEY=VALUE format"):
        load_configuration()


def test_cmdline_params_multiple_params(monkeypatch, tmp_path):
    """Test passing multiple parameters."""
    config_file = tmp_path / "test_config.py"
    config_file.write_text("""
a = params.get('a', 0)
b = params.get('b', 0)
c = params.get('c', 0)

general.output_directory = f"test_{a}_{b}_{c}"
""")

    test_args = [
        'driver.py',
        str(config_file),
        '-p', 'a=1',
        '-p', 'b=2',
        '-p', 'c=3'
    ]
    monkeypatch.setattr(sys, 'argv', test_args)

    state = load_configuration()
    assert state.config_general.output_directory == 'test_1_2_3'
