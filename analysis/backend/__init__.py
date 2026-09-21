"""Backend system for QHAT.

This module provides the infrastructure for pluggable quantum computing backends.
Backends wrap various frameworks (Qualtran, PennyLane, Qiskit, etc.) with a
unified interface.

Main components:
- Backend: Protocol defining required backend operations
- Unitary: Abstract base class for framework-agnostic quantum operators
- BackendRegistry: System for discovering and loading backends
- get_backend(): Convenience function for obtaining backend instances

Example:
    >>> from qhat.analysis.backend import get_backend
    >>> backend = get_backend("pennylane")
    >>> unitary = backend.encode_pauli_trotter(...)
    >>> resources = unitary.estimate_resources()
"""

from typing import Dict, List, Optional, Any, Type
import logging

from qhat.analysis.backend.protocol import Backend
from qhat.analysis.backend.base import Unitary
from qhat.analysis.backend.types import ResourceEstimate, UnsupportedOperationError

logger = logging.getLogger(__name__)


class BackendRegistry:
    """Registry for discovering and loading backends.

    This singleton class manages backend registration and instantiation.
    Backends can be registered explicitly or loaded lazily on first use.

    Lazy loading allows QHAT to work without all backend dependencies installed.
    A backend is only imported when requested, and import errors are gracefully
    handled.
    """

    _backends: Dict[str, Type] = {}
    _loaded_backends: Dict[str, Any] = {}

    @classmethod
    def register(cls, name: str, backend_class: Type):
        """Register a backend class.

        Args:
            name: Backend identifier (e.g., "qualtran")
            backend_class: Class implementing Backend protocol
        """
        cls._backends[name] = backend_class
        logger.debug(f"Registered backend: {name}")

    @classmethod
    def get_backend(cls, name: str, **config) -> Backend:
        """Get or create a backend instance.

        Backends are cached based on name and configuration. Requesting the
        same backend with the same config returns the cached instance.

        Args:
            name: Backend identifier
            **config: Backend-specific configuration options

        Returns:
            Backend instance

        Raises:
            ValueError: If backend not found or failed to load
            ImportError: If backend dependencies not installed

        Example:
            >>> backend = BackendRegistry.get_backend("qualtran")
            >>> backend = BackendRegistry.get_backend(
            ...     "pennylane",
            ...     device="default.qubit"
            ... )
        """
        # Create cache key from name and config
        # Use frozenset for hashable config representation
        try:
            cache_key = f"{name}:{hash(frozenset(config.items()))}"
        except TypeError:
            # Config contains unhashable values, don't cache
            cache_key = None

        # Check cache
        if cache_key and cache_key in cls._loaded_backends:
            logger.debug(f"Using cached backend: {name}")
            return cls._loaded_backends[cache_key]

        # Try lazy import if not registered
        if name not in cls._backends:
            cls._try_lazy_load(name)

        if name not in cls._backends:
            available = cls.list_available()
            raise ValueError(
                f"Backend '{name}' not found. "
                f"Available backends: {available}"
            )

        # Instantiate backend
        backend_class = cls._backends[name]
        try:
            backend = backend_class(**config)
            logger.info(f"Loaded backend: {name}")

            # Cache if possible
            if cache_key:
                cls._loaded_backends[cache_key] = backend

            return backend

        except Exception as e:
            logger.error(f"Failed to load backend '{name}': {e}")
            raise ValueError(f"Failed to load backend '{name}': {e}") from e

    @classmethod
    def _try_lazy_load(cls, name: str):
        """Attempt to lazy-load a backend.

        This tries to import the backend module and register it. If the import
        fails (dependencies not installed), we log a warning but don't error.

        Args:
            name: Backend identifier to try loading
        """
        try:
            if name == "qualtran":
                from qhat.analysis.backend.qualtran_backend import QualtranBackend
                cls.register("qualtran", QualtranBackend)
            elif name == "pennylane":
                from qhat.analysis.backend.pennylane_backend import PennyLaneBackend
                cls.register("pennylane", PennyLaneBackend)
            elif name == "qiskit":
                from qhat.analysis.backend.qiskit_backend import QiskitBackend
                cls.register("qiskit", QiskitBackend)
            else:
                logger.debug(f"No lazy loader defined for backend: {name}")

        except ImportError as e:
            logger.warning(
                f"Backend '{name}' cannot be loaded (dependencies not installed): {e}"
            )

    @classmethod
    def list_available(cls) -> List[str]:
        """List all registered backends.

        Attempts to lazy-load known backends first, then returns the list
        of successfully registered backends.

        Returns:
            List of backend names

        Example:
            >>> BackendRegistry.list_available()
            ['qualtran', 'pennylane', 'qiskit']
        """
        # Try to lazy-load all known backends
        for name in ["qualtran", "pennylane", "qiskit"]:
            if name not in cls._backends:
                cls._try_lazy_load(name)

        return sorted(cls._backends.keys())

    @classmethod
    def check_backend_available(cls, name: str) -> bool:
        """Check if a backend is available (dependencies installed).

        Args:
            name: Backend identifier

        Returns:
            True if backend can be loaded, False otherwise

        Example:
            >>> if BackendRegistry.check_backend_available("pennylane"):
            ...     backend = BackendRegistry.get_backend("pennylane")
        """
        try:
            cls.get_backend(name)
            return True
        except (ValueError, ImportError):
            return False

    @classmethod
    def clear_cache(cls):
        """Clear cached backend instances.

        Useful for testing or when backend configuration changes.
        """
        cls._loaded_backends.clear()
        logger.debug("Cleared backend cache")


# Convenience function
def get_backend(name: str, **config) -> Backend:
    """Get a backend instance.

    This is the primary entry point for obtaining backends in QHAT.

    Args:
        name: Backend identifier ('qualtran', 'pennylane', 'qiskit')
        **config: Backend-specific configuration

    Returns:
        Backend instance

    Raises:
        ValueError: If backend not found
        ImportError: If backend dependencies not installed

    Example:
        >>> from qhat.analysis.backend import get_backend
        >>> backend = get_backend("qualtran")
        >>> print(f"Using {backend.name} with capabilities: {backend.capabilities}")
    """
    return BackendRegistry.get_backend(name, **config)


# Export public API
__all__ = [
    "Backend",
    "Unitary",
    "ResourceEstimate",
    "UnsupportedOperationError",
    "BackendRegistry",
    "get_backend",
]
