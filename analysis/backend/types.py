"""Shared types for backend system."""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional


@dataclass
class ResourceEstimate:
    """Unified resource estimation across backends.

    This dataclass provides a common format for quantum resource estimates
    regardless of which backend computed them.

    Attributes:
        num_qubits: Number of qubits required
        t_gates: Number of T gates (primary cost metric for fault-tolerant QC)
        clifford_gates: Number of Clifford gates (H, S, CNOT, etc.)
        rotation_gates: Number of arbitrary rotation gates
        measurements: Number of measurement operations
        depth: Circuit depth (critical path length), if computable
        two_qubit_gates: Number of two-qubit gates, if tracked separately
        backend_specific: Dict for backend-specific metrics not in common schema
    """
    num_qubits: int
    t_gates: int
    clifford_gates: int
    rotation_gates: int = 0
    measurements: int = 0
    depth: Optional[int] = None
    two_qubit_gates: Optional[int] = None
    backend_specific: Dict[str, Any] = field(default_factory=dict)

    def total_gates(self) -> int:
        """Total gate count (excluding measurements)."""
        return self.t_gates + self.clifford_gates + self.rotation_gates

    def __str__(self) -> str:
        lines = [
            "Resource Estimate:",
            f"  Qubits: {self.num_qubits}",
            f"  T gates: {self.t_gates}",
            f"  Clifford gates: {self.clifford_gates}",
        ]
        if self.rotation_gates > 0:
            lines.append(f"  Rotation gates: {self.rotation_gates}")
        if self.measurements > 0:
            lines.append(f"  Measurements: {self.measurements}")
        if self.depth is not None:
            lines.append(f"  Depth: {self.depth}")
        if self.two_qubit_gates is not None:
            lines.append(f"  Two-qubit gates: {self.two_qubit_gates}")
        if self.backend_specific:
            lines.append(f"  Backend-specific: {self.backend_specific}")
        return "\n".join(lines)


class UnsupportedOperationError(Exception):
    """Raised when backend doesn't support requested operation.

    Backends are not required to implement all operations. When a user
    requests an unsupported operation, the backend should raise this
    exception with a clear message about what's not supported and
    optionally which backends do support it.
    """
    pass
