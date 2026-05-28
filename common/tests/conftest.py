"""
Pytest configuration for common tests.

Configures logging before tests run so that custom VERBOSE level is available.
"""

import pytest
from qhat.common.logging_utils import configure_logging


@pytest.fixture(scope="session", autouse=True)
def configure_logging_for_tests():
    """Configure logging once for all tests in the session."""
    configure_logging(level="warning")  # Set to warning to reduce test output
    yield
