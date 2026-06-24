"""Shared test configuration and fixtures."""
from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# Markers
# ---------------------------------------------------------------------------

def pytest_configure(config):
    config.addinivalue_line(
        "markers", "integration: marks tests that require GPAW data files"
    )


def pytest_collection_modifyitems(config, items):
    """Skip integration tests unless the user explicitly selects them via -m."""
    markexpr = getattr(config.option, "markexpr", "")
    if "integration" in markexpr:
        return  # user asked for integration tests; respect the filter
    skip = pytest.mark.skip(reason="integration test — run with: pytest -m integration")
    for item in items:
        if item.get_closest_marker("integration"):
            item.add_marker(skip)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return __import__("numpy").random.default_rng(42)
