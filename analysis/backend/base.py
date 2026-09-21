"""Base classes for framework-agnostic quantum operators."""

from abc import ABC, abstractmethod
from typing import Any, Optional
import numpy as np

from qhat.analysis.backend.types import ResourceEstimate


class Unitary(ABC):
    """Abstract base class for framework-agnostic unitary operators.

    This represents a quantum operator that can be composed, controlled,
    and analyzed. Concrete implementations wrap framework-specific objects
    (Qualtran Bloq, PennyLane QNode, Qiskit QuantumCircuit, etc.).

    The Unitary abstraction allows QHAT's front-end code to work with quantum
    operators without depending on any specific framework. Backend implementations
    provide concrete Unitary subclasses that delegate to their native representations.

    Example:
        >>> backend = get_backend("qualtran")
        >>> unitary = backend.encode_pauli_trotter(...)
        >>> controlled_unitary = unitary.controlled(num_controls=2)
        >>> resources = unitary.estimate_resources()
        >>> print(f"T gates: {resources.t_gates}")
    """

    def __init__(self, backend_name: str, num_qubits: int):
        """Initialize unitary.

        Args:
            backend_name: Identifier of the backend that created this
            num_qubits: Number of qubits this operator acts on
        """
        self._backend_name = backend_name
        self._num_qubits = num_qubits

    @property
    def backend_name(self) -> str:
        """Name of backend that created this unitary."""
        return self._backend_name

    @property
    def num_qubits(self) -> int:
        """Number of qubits this operator acts on."""
        return self._num_qubits

    @abstractmethod
    def controlled(self, num_controls: int = 1) -> 'Unitary':
        """Generate controlled version of this unitary.

        Returns a new unitary U' that applies U controlled on num_controls qubits:
        - If all controls are |1>, apply U to target qubits
        - Otherwise, apply identity

        The returned unitary acts on (num_controls + self.num_qubits) qubits.

        Args:
            num_controls: Number of control qubits to add

        Returns:
            New Unitary with control qubits added

        Raises:
            NotImplementedError: If backend doesn't support controlled operations
        """
        ...

    @abstractmethod
    def adjoint(self) -> 'Unitary':
        """Generate adjoint (Hermitian conjugate) of this unitary.

        Returns U† such that U† U = U U† = I.

        Returns:
            Adjoint unitary acting on same qubits

        Raises:
            NotImplementedError: If backend doesn't support adjoint
        """
        ...

    @abstractmethod
    def power(self, exponent: float) -> 'Unitary':
        """Raise unitary to a power: U^exponent.

        Critical for "fast-forwardable" phase estimation where U^(2^k) can
        be computed more efficiently than applying U repeatedly.

        For time evolution operators U = exp(-iHt), this often means:
        U^k = exp(-iHkt) which may have similar cost to U itself.

        Args:
            exponent: Power to raise operator to

        Returns:
            U^exponent

        Raises:
            NotImplementedError: If backend doesn't support power operation
        """
        ...

    @abstractmethod
    def to_matrix(self,
                  max_qubits: int = 20,
                  sparse: bool = False) -> np.ndarray:
        """Convert to explicit matrix representation.

        Constructs the 2^N x 2^N unitary matrix where N = self.num_qubits.
        For large systems this is memory-intensive (2^N)^2 * 16 bytes.

        Args:
            max_qubits: Safety limit to prevent memory overflow
            sparse: Return scipy sparse matrix if True, dense numpy array if False

        Returns:
            Matrix representation of this unitary

        Raises:
            ValueError: If num_qubits > max_qubits
            MemoryError: If allocation fails
            NotImplementedError: If backend doesn't support matrix conversion
        """
        ...

    @abstractmethod
    def estimate_resources(self) -> ResourceEstimate:
        """Estimate quantum resources required to implement this operator.

        Analyzes the circuit/algorithm and counts gates, qubits, depth, etc.
        The accuracy and detail of estimates vary by backend.

        Returns:
            Resource estimates (qubits, gates, depth, etc.)

        Raises:
            NotImplementedError: If backend doesn't support resource estimation
        """
        ...

    @abstractmethod
    def get_native_object(self) -> Any:
        """Get the underlying framework-specific object.

        Use this escape hatch when you need to perform framework-specific
        operations not covered by the Unitary abstract interface.

        Returns:
            Native object (Bloq, QNode, QuantumCircuit, etc.)

        Example:
            >>> unitary = backend.encode_pauli_trotter(...)
            >>> if unitary.backend_name == "qualtran":
            ...     bloq = unitary.get_native_object()  # Returns Bloq
            ...     # Use Qualtran-specific features
            ...     t_complexity = bloq.t_complexity()
        """
        ...

    def __repr__(self) -> str:
        return (f"<{self.__class__.__name__} "
                f"backend={self.backend_name} "
                f"qubits={self.num_qubits}>")

    # Optional: Operators for composition
    # These provide convenience but aren't required by all backends

    def __matmul__(self, other: 'Unitary') -> 'Unitary':
        """Compose unitaries: (self @ other) means apply other then self.

        Not all backends may support this. Default implementation raises NotImplementedError.

        Args:
            other: Unitary to compose with

        Returns:
            Composed unitary

        Raises:
            NotImplementedError: If backend doesn't support composition
            TypeError: If other is from incompatible backend
        """
        raise NotImplementedError(
            f"{self.backend_name} backend doesn't support unitary composition via @"
        )
