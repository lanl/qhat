"""PennyLane backend implementation.

This backend provides QHAT functionality using PennyLane, a popular framework
for quantum machine learning and differentiable quantum computing.
"""

from typing import Dict, Set, Any, List, Tuple
import logging
import numpy as np

from qhat.analysis.backend.protocol import Backend
from qhat.analysis.backend.base import Unitary
from qhat.analysis.backend.types import ResourceEstimate, UnsupportedOperationError

logger = logging.getLogger(__name__)


class PennyLaneUnitary(Unitary):
    """Wrapper around PennyLane quantum operations.

    Unlike Qualtran's Bloq which is a high-level circuit object, PennyLane
    represents circuits as sequences of operations applied to a QNode. We
    store the operation list and can reconstruct QNodes as needed.
    """

    def __init__(self, operations: List[Any], num_qubits: int, device_name: str = "default.qubit"):
        """Initialize from PennyLane operations.

        Args:
            operations: List of PennyLane operations
            num_qubits: Number of qubits
            device_name: PennyLane device to use
        """
        super().__init__("pennylane", num_qubits)
        self._operations = operations
        self._device_name = device_name

    def controlled(self, num_controls: int = 1) -> 'PennyLaneUnitary':
        """Generate controlled operations using qml.ctrl."""
        import pennylane as qml

        try:
            # Wrap each operation in qml.ctrl
            # Control qubits are prepended (indices 0 to num_controls-1)
            # Target qubits are shifted by num_controls
            controlled_ops = []

            for op in self._operations:
                # Create controlled version
                # qml.ctrl wraps the operation with control wires
                ctrl_op = qml.ctrl(op, control=list(range(num_controls)))
                controlled_ops.append(ctrl_op)

            return PennyLaneUnitary(
                controlled_ops,
                self.num_qubits + num_controls,
                self._device_name
            )

        except Exception as e:
            logger.error(f"Failed to create controlled PennyLane operations: {e}")
            raise NotImplementedError(f"PennyLane backend failed to create controlled unitary: {e}")

    def adjoint(self) -> 'PennyLaneUnitary':
        """Generate adjoint using qml.adjoint."""
        import pennylane as qml

        try:
            # Apply adjoint to each operation and reverse order
            adjoint_ops = [qml.adjoint(op) for op in reversed(self._operations)]

            return PennyLaneUnitary(adjoint_ops, self.num_qubits, self._device_name)

        except Exception as e:
            logger.error(f"Failed to create adjoint PennyLane operations: {e}")
            raise NotImplementedError(f"PennyLane backend failed to create adjoint: {e}")

    def power(self, exponent: float) -> 'PennyLaneUnitary':
        """Raise to power by repeating operations.

        For integer powers, we repeat the operations. For fractional powers,
        this is approximate and may not be accurate.
        """
        import pennylane as qml

        if not isinstance(exponent, int) or exponent < 0:
            logger.warning(
                f"PennyLane power operation with exponent={exponent} may be inaccurate. "
                "Consider using a backend that supports symbolic powers."
            )

        try:
            # For integer exponents, repeat operations
            if isinstance(exponent, int) and exponent >= 0:
                powered_ops = self._operations * exponent
                return PennyLaneUnitary(powered_ops, self.num_qubits, self._device_name)
            else:
                # Fractional or negative powers not well-supported
                raise NotImplementedError(
                    f"PennyLane backend doesn't support non-integer powers (exponent={exponent})"
                )

        except Exception as e:
            logger.error(f"Failed to raise PennyLane unitary to power: {e}")
            raise NotImplementedError(f"PennyLane power operation failed: {e}")

    def to_matrix(self, max_qubits: int = 20, sparse: bool = False) -> np.ndarray:
        """Convert to matrix using qml.matrix."""
        if self.num_qubits > max_qubits:
            raise ValueError(
                f"Cannot convert {self.num_qubits}-qubit operator to matrix "
                f"(exceeds max_qubits={max_qubits})"
            )

        try:
            import pennylane as qml

            # Create a QNode that applies the operations
            dev = qml.device(self._device_name, wires=self.num_qubits)

            @qml.qnode(dev)
            def circuit():
                for op in self._operations:
                    qml.apply(op)
                return qml.state()

            # Get the unitary matrix
            matrix = qml.matrix(circuit)()

            if sparse:
                from scipy.sparse import csr_matrix
                return csr_matrix(matrix)

            return matrix

        except Exception as e:
            logger.error(f"Failed to convert PennyLane operations to matrix: {e}")
            raise NotImplementedError(f"PennyLane matrix conversion failed: {e}")

    def estimate_resources(self) -> ResourceEstimate:
        """Estimate resources by counting operations.

        PennyLane doesn't have built-in T-complexity analysis, so we count
        operations by type and make rough estimates.
        """
        import pennylane as qml

        t_gates = 0
        clifford_gates = 0
        rotation_gates = 0
        two_qubit_gates = 0

        for op in self._operations:
            op_name = op.name if hasattr(op, 'name') else str(type(op).__name__)

            # Categorize gates
            if op_name in ['T', 'Adjoint(T)']:
                t_gates += 1
            elif op_name in ['Hadamard', 'PauliX', 'PauliY', 'PauliZ', 'S', 'CNOT', 'CZ', 'SWAP']:
                clifford_gates += 1
            elif op_name in ['RX', 'RY', 'RZ', 'Rot', 'PhaseShift', 'PauliRot']:
                rotation_gates += 1
                # Estimate T-gate cost of rotations (rough approximation)
                # Arbitrary rotations typically require ~50 T gates for ε=10^-3 precision
                t_gates += 50

            # Count two-qubit gates
            if hasattr(op, 'wires') and len(op.wires) == 2:
                two_qubit_gates += 1

        return ResourceEstimate(
            num_qubits=self.num_qubits,
            t_gates=t_gates,
            clifford_gates=clifford_gates,
            rotation_gates=rotation_gates,
            two_qubit_gates=two_qubit_gates,
            backend_specific={
                'pennylane_num_operations': len(self._operations),
                'pennylane_device': self._device_name,
                'note': 'T-gate counts for rotations are rough estimates'
            }
        )

    def get_native_object(self) -> List[Any]:
        """Return the list of PennyLane operations."""
        return self._operations


class PennyLaneBackend:
    """Backend implementation using PennyLane.

    PennyLane is a framework for differentiable quantum computing with strong
    support for quantum machine learning, automatic differentiation, and
    execution on various hardware backends.
    """

    def __init__(self, device: str = "default.qubit", **config):
        """Initialize PennyLane backend.

        Args:
            device: PennyLane device name (default: "default.qubit")
            **config: Additional device configuration
        """
        self.device_name = device
        self.config = config
        self._validate_imports()

    def _validate_imports(self):
        """Ensure PennyLane is available."""
        try:
            import pennylane as qml
        except ImportError as e:
            raise ImportError(
                f"PennyLane backend requires pennylane package: {e}"
            ) from e

    @property
    def name(self) -> str:
        return "pennylane"

    @property
    def capabilities(self) -> Set[str]:
        return {
            "pauli_trotter",
            "pauli_lcu",  # Limited support
            "textbook_qpe",  # Possible but not optimized
            "resource_estimation",  # Basic counting
            "matrix_conversion",
            "controlled_operations",
            "adjoint_operations",
        }

    def _pauli_tuple_to_pennylane(
        self,
        pauli_tuple: Tuple[Tuple[int, str], ...],
        num_qubits: int
    ) -> str:
        """Convert sparse Pauli tuple to PennyLane Pauli word.

        Args:
            pauli_tuple: Sparse format ((qubit_idx, 'X'/'Y'/'Z'), ...)
            num_qubits: Total number of qubits

        Returns:
            Dense Pauli string for PennyLane (e.g., "XIZI")
        """
        dense = ['I'] * num_qubits
        for idx, op in pauli_tuple:
            dense[idx] = op
        return ''.join(dense)

    def encode_pauli_trotter(
        self,
        pauli_strings: Dict[tuple, float],
        trotter_order: str,
        evolution_time: float,
        num_steps: int,
        num_qubits: int,
        **kwargs
    ) -> PennyLaneUnitary:
        """Encode using Trotterization with PennyLane's PauliRot.

        PennyLane provides qml.PauliRot which efficiently implements
        exp(-i * coef * t * PauliString) for any Pauli string.
        """
        import pennylane as qml

        operations = []
        dt = evolution_time / num_steps

        # Map trotter order to formula
        if trotter_order.lower() == "first order":
            # First order: (e^(-iH_1 dt) ... e^(-iH_N dt))^n
            for step in range(num_steps):
                for pauli_tuple, coef in pauli_strings.items():
                    pauli_word = self._pauli_tuple_to_pennylane(pauli_tuple, num_qubits)
                    # PauliRot implements exp(-i * angle/2 * PauliString)
                    # We want exp(-i * coef * dt * PauliString)
                    # So angle = 2 * coef * dt
                    operations.append(
                        qml.PauliRot(2 * coef * dt, pauli_word, wires=range(num_qubits))
                    )

        elif trotter_order.lower() == "second order":
            # Second order Suzuki formula: forward then backward
            for step in range(num_steps):
                # Forward sweep
                for pauli_tuple, coef in pauli_strings.items():
                    pauli_word = self._pauli_tuple_to_pennylane(pauli_tuple, num_qubits)
                    operations.append(
                        qml.PauliRot(coef * dt, pauli_word, wires=range(num_qubits))
                    )
                # Backward sweep
                for pauli_tuple, coef in reversed(list(pauli_strings.items())):
                    pauli_word = self._pauli_tuple_to_pennylane(pauli_tuple, num_qubits)
                    operations.append(
                        qml.PauliRot(coef * dt, pauli_word, wires=range(num_qubits))
                    )

        elif trotter_order.lower() == "fourth order":
            # Fourth order Suzuki formula (more complex)
            # Using standard 4th order formula with 5 sub-steps
            # S4 = S2(p)^2 S2(1-4p) S2(p)^2 where p = 1/(4-4^(1/3))
            p = 1.0 / (4.0 - 4.0 ** (1.0 / 3.0))

            def second_order_step(time_fraction):
                """Apply second-order step with given time scaling."""
                scaled_dt = dt * time_fraction
                # Forward
                for pauli_tuple, coef in pauli_strings.items():
                    pauli_word = self._pauli_tuple_to_pennylane(pauli_tuple, num_qubits)
                    operations.append(
                        qml.PauliRot(coef * scaled_dt, pauli_word, wires=range(num_qubits))
                    )
                # Backward
                for pauli_tuple, coef in reversed(list(pauli_strings.items())):
                    pauli_word = self._pauli_tuple_to_pennylane(pauli_tuple, num_qubits)
                    operations.append(
                        qml.PauliRot(coef * scaled_dt, pauli_word, wires=range(num_qubits))
                    )

            for step in range(num_steps):
                second_order_step(p)
                second_order_step(p)
                second_order_step(1.0 - 4.0 * p)
                second_order_step(p)
                second_order_step(p)

        else:
            raise ValueError(f"Unsupported Trotter order for PennyLane: {trotter_order}")

        return PennyLaneUnitary(operations, num_qubits, self.device_name)

    def encode_pauli_lcu(
        self,
        pauli_strings: Dict[tuple, float],
        num_qubits: int,
        prepare_type: str = 'AS',
        probability_eps: float = 0.002,
        **kwargs
    ) -> PennyLaneUnitary:
        """LCU encoding with PennyLane.

        LCU requires SELECT and PREPARE oracles which are more naturally
        expressed in frameworks like Qualtran. PennyLane can implement them
        but it's not the most efficient approach.

        For now, we raise UnsupportedOperationError and suggest using Qualtran.
        """
        raise UnsupportedOperationError(
            "PennyLane backend doesn't yet support LCU block encoding. "
            "Use 'qualtran' backend for LCU-based methods."
        )

    def build_textbook_qpe(
        self,
        unitary: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> PennyLaneUnitary:
        """Build textbook QPE with PennyLane.

        This requires:
        1. Hadamards on phase register
        2. Controlled-U^(2^k) operations
        3. Inverse QFT

        PennyLane can do this but it's not optimized. We provide basic implementation.
        """
        import pennylane as qml

        if not isinstance(unitary, PennyLaneUnitary):
            raise TypeError(
                f"PennyLane backend requires PennyLaneUnitary, got {type(unitary).__name__}"
            )

        operations = []
        total_qubits = num_phase_qubits + unitary.num_qubits

        # Phase register: qubits 0 to num_phase_qubits-1
        # Target register: qubits num_phase_qubits to total_qubits-1

        # Step 1: Hadamards on phase register
        for i in range(num_phase_qubits):
            operations.append(qml.Hadamard(wires=i))

        # Step 2: Controlled-U^(2^k) operations
        # This is tricky because we need to raise U to powers efficiently
        for k in range(num_phase_qubits):
            # We need controlled-U^(2^k) controlled on qubit k
            # For simplicity, we repeat U operations 2^k times
            # This is NOT fast-forwardable and will be expensive
            power = 2 ** k

            # Wrap unitary operations with control on qubit k
            for _ in range(power):
                for op in unitary.get_native_object():
                    ctrl_op = qml.ctrl(op, control=k)
                    operations.append(ctrl_op)

        # Step 3: Inverse QFT on phase register
        # Implement basic QFT inverse
        for i in range(num_phase_qubits // 2):
            operations.append(qml.SWAP(wires=[i, num_phase_qubits - 1 - i]))

        for i in range(num_phase_qubits):
            operations.append(qml.Hadamard(wires=i))
            for j in range(i):
                angle = -np.pi / (2 ** (i - j))
                # Controlled phase rotation
                operations.append(qml.ctrl(qml.PhaseShift(angle, wires=i), control=j))

        logger.warning(
            "PennyLane textbook QPE is not optimized and may be very expensive. "
            "Consider using 'qualtran' backend for efficient QPE."
        )

        return PennyLaneUnitary(operations, total_qubits, self.device_name)

    def build_qubitized_qpe(
        self,
        block_encoding: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> PennyLaneUnitary:
        """Qubitized QPE not supported in PennyLane backend."""
        raise UnsupportedOperationError(
            "PennyLane backend doesn't support qubitized QPE. "
            "Use 'qualtran' backend for qubitization-based methods."
        )
