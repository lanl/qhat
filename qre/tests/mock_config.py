"""
Mock configuration for testing.

Provides a simple mock GeneralConfiguration object that can be used
in tests without requiring the full QRE infrastructure.
"""


class MockConfig:
    """Mock configuration object for testing."""

    def log_verbose(self, message):
        """Mock log_verbose - does nothing (suppresses output during tests)."""
        pass

    def log(self, message):
        """Mock log - does nothing (suppresses output during tests)."""
        pass


# Singleton instance for convenience
mock_config = MockConfig()
