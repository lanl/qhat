"""Integration-test conftest: restore the real GPAW package before each test.

The top-level conftest stubs ``gpaw`` as a bare namespace module for unit-test
isolation.  Integration tests need the real GPAW installation, so we swap the
stub back out just before each integration test runs.
"""
from __future__ import annotations

import importlib
import sys
import types

import pytest


def pytest_runtest_setup(item):
    """Before each integration test, ensure the real GPAW package is loaded."""
    if not item.get_closest_marker("integration"):
        return

    # Remove any bare-stub entries for gpaw and its sub-namespaces.
    # A stub is a plain types.ModuleType with no __version__ or __file__.
    for key in list(sys.modules):
        if key == "gpaw" or key.startswith("gpaw."):
            mod = sys.modules[key]
            if isinstance(mod, types.ModuleType) and not getattr(mod, "__version__", None) and not getattr(mod, "__file__", None):
                del sys.modules[key]

    # Now import the real package; skip the test if GPAW is not installed.
    try:
        importlib.import_module("gpaw")
    except ImportError:
        pytest.skip("GPAW not installed")
