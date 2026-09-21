"""Qualtran backend implementation.

This backend wraps QHAT's existing Qualtran and pyLIQTR functionality,
preserving all existing behavior while conforming to the Backend protocol.
"""

from typing import Dict, Set, Any, Optional
import logging
import math
import numpy as np

from qhat.analysis.backend.protocol import Backend
from qhat.analysis.backend.base import Unitary
from qhat.analysis.backend.types import ResourceEstimate, UnsupportedOperationError

logger = logging.getLogger(__name__)


class QualtranUnitary(Unitary):
    """Wrapper around Qualtran Bloq objects.

    This class adapts Qualtran's Bloq to QHAT's Unitary interface, allowing
    the front-end code to work with Bloqs without direct Qualtran dependencies.
    """

    def __init__(self, bloq: 'Bloq', num_qubits: int):
        """Initialize from a Qualtran Bloq.

        Args:
            bloq: Qualtran Bloq object
            num_qubits: Number of qubits (must match bloq's signature)
        """
        super().__init__("qualtran", num_qubits)
        self._bloq = bloq

    def controlled(self, num_controls: int = 1) -> 'QualtranUnitary':
        """Generate controlled Bloq using Qualtran's controlled() method."""
        try:
            controlled_bloq = self._bloq.controlled(n=num_controls)
            # Controlled bloq adds num_controls qubits
            return QualtranUnitary(controlled_bloq, self.num_qubits + num_controls)
        except Exception as e:
            logger.error(f"Failed to create controlled bloq: {e}")
            raise NotImplementedError(f"Qualtran backend failed to create controlled bloq: {e}")

    def adjoint(self) -> 'QualtranUnitary':
        """Generate adjoint Bloq using Qualtran's adjoint() method."""
        try:
            adjoint_bloq = self._bloq.adjoint()
            return QualtranUnitary(adjoint_bloq, self.num_qubits)
        except Exception as e:
            logger.error(f"Failed to create adjoint bloq: {e}")
            raise NotImplementedError(f"Qualtran backend failed to create adjoint bloq: {e}")

    def power(self, exponent: float) -> 'QualtranUnitary':
        """Raise to power via cirq.pow.

        Qualtran Bloqs can be raised to powers using Cirq's pow protocol.
        This enables fast-forwardable phase estimation.
        """
        try:
            import cirq
            powered_bloq = cirq.pow(self._bloq, exponent)
            return QualtranUnitary(powered_bloq, self.num_qubits)
        except Exception as e:
            logger.error(f"Failed to raise bloq to power {exponent}: {e}")
            raise NotImplementedError(f"Qualtran backend failed to compute power: {e}")

    def to_matrix(self, max_qubits: int = 20, sparse: bool = False) -> np.ndarray:
        """Convert Bloq to matrix via Cirq interop."""
        if self.num_qubits > max_qubits:
            raise ValueError(
                f"Cannot convert {self.num_qubits}-qubit operator to matrix "
                f"(exceeds max_qubits={max_qubits}). "
                f"Matrix would require {2**(2*self.num_qubits) * 16 / 1e9:.1f} GB."
            )

        try:
            from qualtran.cirq_interop import BloqAsCirqGate
            import cirq

            # Convert Bloq to Cirq gate
            gate = BloqAsCirqGate(self._bloq)

            # Get unitary matrix
            unitary_matrix = cirq.unitary(gate)

            if sparse:
                from scipy.sparse import csr_matrix
                return csr_matrix(unitary_matrix)

            return unitary_matrix

        except Exception as e:
            logger.error(f"Failed to convert bloq to matrix: {e}")
            raise NotImplementedError(f"Qualtran backend failed to generate matrix: {e}")

    def estimate_resources(self) -> ResourceEstimate:
        """Estimate resources using Qualtran's t_complexity.

        Qualtran provides detailed T-complexity analysis which we convert
        to the unified ResourceEstimate format.
        """
        try:
            t_complexity = self._bloq.t_complexity()

            return ResourceEstimate(
                num_qubits=self.num_qubits,
                t_gates=t_complexity.t,
                clifford_gates=t_complexity.clifford,
                rotation_gates=getattr(t_complexity, 'rotations', 0),
                backend_specific={
                    'qualtran_t_complexity': str(t_complexity),
                    'qualtran_repr': repr(t_complexity)
                }
            )

        except Exception as e:
            logger.warning(f"Failed to estimate resources via t_complexity: {e}")
            # Return minimal estimate
            return ResourceEstimate(
                num_qubits=self.num_qubits,
                t_gates=0,
                clifford_gates=0,
                backend_specific={'error': str(e)}
            )

    def get_native_object(self) -> 'Bloq':
        """Return the underlying Qualtran Bloq."""
        return self._bloq


class QualtranBackend:
    """Backend implementation using Qualtran and pyLIQTR.

    This backend preserves all existing QHAT functionality while conforming
    to the Backend protocol. It delegates to existing code in qhat.common
    and qhat.analysis modules where possible.
    """

    def __init__(self, **config):
        """Initialize Qualtran backend.

        Args:
            **config: Configuration options (currently unused, reserved for future)
        """
        self.config = config
        self._validate_imports()

    def _validate_imports(self):
        """Ensure required packages are available."""
        try:
            import qualtran
            import pyLIQTR
            import cirq
        except ImportError as e:
            raise ImportError(
                f"Qualtran backend requires qualtran, pyLIQTR, and cirq packages: {e}"
            ) from e

    @property
    def name(self) -> str:
        return "qualtran"

    @property
    def capabilities(self) -> Set[str]:
        return {
            "pauli_trotter",
            "pauli_lcu",
            "textbook_qpe",
            "qubitized_qpe",
            "double_factorization",
            "resource_estimation",
            "matrix_conversion",
            "controlled_operations",
            "adjoint_operations",
            "power_operations",
        }

    def encode_pauli_trotter(
        self,
        pauli_strings: Dict[tuple, float],
        trotter_order: str,
        evolution_time: float,
        num_steps: int,
        num_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Encode using QHAT's existing Trotter implementation.

        Delegates to qhat.common.trotter_flattened which contains the optimized
        Trotterization code.
        """
        # Import existing implementation
        trotter_impl = kwargs.get('trotter_implementation', 'flattened')

        if trotter_impl == 'flattened':
            from qhat.common.trotter_flattened import build_ramped_trotterized_unitary
        elif trotter_impl == 'original':
            from qhat.common.trotter_original import build_ramped_trotterized_unitary
        else:
            raise ValueError(f"Unknown trotter_implementation: {trotter_impl}")

        # Convert to format expected by existing code (list of (pauli, coef) tuples)
        pauli_items = list(pauli_strings.items())

        # Build the Bloq using existing code
        try:
            if trotter_impl == 'flattened':
                bloq = build_ramped_trotterized_unitary(
                    pauli_items,
                    trotter_order,
                    evolution_time,
                    num_steps,
                    combine_terms=kwargs.get('combine_terms', True),
                    tensor_contraction_method=kwargs.get('tensor_contraction_method', None)
                )
            else:  # original
                bloq = build_ramped_trotterized_unitary(
                    pauli_items,
                    trotter_order,
                    evolution_time,
                    num_steps
                )

            return QualtranUnitary(bloq, num_qubits)

        except Exception as e:
            logger.error(f"Failed to build Trotter circuit: {e}")
            raise ValueError(f"Qualtran Trotter encoding failed: {e}") from e

    def encode_pauli_lcu(
        self,
        pauli_strings: Dict[tuple, float],
        num_qubits: int,
        prepare_type: str = 'AS',
        probability_eps: float = 0.002,
        **kwargs
    ) -> QualtranUnitary:
        """Encode using LCU block encoding.

        Uses the PauliStringLCU implementation from qhat.analysis.unitary.
        """
        # Import from existing module
        from qhat.analysis.unitary import PauliStringLCU

        # Create a minimal hamiltonian-like object that PauliStringLCU expects
        class _HamiltonianAdapter:
            """Adapter to make pauli_strings dict look like a Hamiltonian."""

            def __init__(self, pauli_dict, nq):
                self._pauli_dict = pauli_dict
                self._nq = nq

            def get_all_pauli_strings(self, return_as='strings'):
                """Convert sparse pauli tuples to dense strings."""
                if return_as == 'strings':
                    import cirq
                    result = {}
                    for pauli_tuple, coef in self._pauli_dict.items():
                        # Convert tuple format to cirq.DensePauliString
                        dense_pauli = ['I'] * self._nq
                        for idx, op in pauli_tuple:
                            dense_pauli[idx] = op
                        dense_str = ''.join(dense_pauli)
                        result[dense_str] = coef
                    return result
                return self._pauli_dict

            def num_qubits(self):
                return self._nq

        hamiltonian_adapter = _HamiltonianAdapter(pauli_strings, num_qubits)

        try:
            bloq = PauliStringLCU(
                hamiltonian_adapter,
                prepare_type=prepare_type,
                probability_eps=probability_eps
            )

            return QualtranUnitary(bloq, num_qubits)

        except Exception as e:
            logger.error(f"Failed to build LCU encoding: {e}")
            raise ValueError(f"Qualtran LCU encoding failed: {e}") from e

    def build_textbook_qpe(
        self,
        unitary: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Build QPE using Qualtran's TextbookQPE.

        Uses the implementation from qualtran.bloqs.phase_estimation.
        """
        from qualtran.bloqs.phase_estimation import TextbookQPE

        if not isinstance(unitary, QualtranUnitary):
            raise TypeError(
                f"Qualtran backend requires QualtranUnitary, got {type(unitary).__name__}"
            )

        try:
            # Get the underlying Bloq
            inner_bloq = unitary.get_native_object()

            # Build QPE circuit
            qpe_bloq = TextbookQPE(inner_bloq, num_phase_qubits)

            # QPE adds num_phase_qubits ancillas to the system
            total_qubits = unitary.num_qubits + num_phase_qubits

            return QualtranUnitary(qpe_bloq, total_qubits)

        except Exception as e:
            logger.error(f"Failed to build textbook QPE: {e}")
            raise ValueError(f"Qualtran textbook QPE failed: {e}") from e

    def build_qubitized_qpe(
        self,
        block_encoding: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Build qubitized QPE.

        Uses qubitization walk operator and the NewQubitizationQPE from
        qhat.analysis.algorithm (which includes a bugfix backported from Qualtran 0.5.0).
        """
        from qualtran.bloqs.qubitization.qubitization_walk_operator import QubitizationWalkOperator
        from qhat.analysis.algorithm import NewQubitizationQPE

        if not isinstance(block_encoding, QualtranUnitary):
            raise TypeError(
                f"Qualtran backend requires QualtranUnitary, got {type(block_encoding).__name__}"
            )

        try:
            # Get the underlying LCU Bloq
            lcu_bloq = block_encoding.get_native_object()

            # LCU bloqs should have _select_gate and _prepare_gate attributes
            if not (hasattr(lcu_bloq, '_select_gate') and hasattr(lcu_bloq, '_prepare_gate')):
                raise ValueError(
                    "Block encoding must be an LCU with select and prepare gates"
                )

            # Build quantum walk operator
            walk_operator = QubitizationWalkOperator(
                lcu_bloq._select_gate,
                lcu_bloq._prepare_gate
            )

            # Build qubitized QPE
            qpe_bloq = NewQubitizationQPE(walk_operator, num_phase_qubits)

            total_qubits = block_encoding.num_qubits + num_phase_qubits

            return QualtranUnitary(qpe_bloq, total_qubits)

        except Exception as e:
            logger.error(f"Failed to build qubitized QPE: {e}")
            raise ValueError(f"Qualtran qubitized QPE failed: {e}") from e
