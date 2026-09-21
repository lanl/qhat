"""Backend protocol definition.

This module defines the Backend protocol that all backend implementations
must satisfy. Using Protocol (PEP 544) allows structural subtyping, meaning
we can wrap existing framework objects without inheritance.
"""

from typing import Protocol, Set, Dict, Any, runtime_checkable

from qhat.analysis.backend.types import UnsupportedOperationError


@runtime_checkable
class Backend(Protocol):
    """Protocol defining required backend operations.

    Backends are not required to implement all methods. For unsupported
    operations, they should raise UnsupportedOperationError with a helpful
    message.

    All backend implementations should:
    1. Define 'name' property returning a string identifier
    2. Define 'capabilities' property returning supported operation names
    3. Implement the operations listed in 'capabilities'
    4. Raise UnsupportedOperationError for operations not in 'capabilities'
    """

    @property
    def name(self) -> str:
        """Backend identifier (e.g., 'qualtran', 'pennylane', 'qiskit').

        This name is used in configuration files and logging.
        """
        ...

    @property
    def capabilities(self) -> Set[str]:
        """Set of supported operation names.

        Common capability strings:
        - 'pauli_trotter': encode_pauli_trotter method
        - 'pauli_lcu': encode_pauli_lcu method
        - 'textbook_qpe': build_textbook_qpe method
        - 'qubitized_qpe': build_qubitized_qpe method
        - 'double_factorization': encode_double_factorization method
        - 'resource_estimation': accurate resource estimation
        - 'matrix_conversion': to_matrix on Unitary objects
        - 'controlled_operations': controlled() on Unitary objects

        Returns:
            Set of capability strings
        """
        ...

    def encode_pauli_trotter(
        self,
        pauli_strings: Dict[tuple, float],
        trotter_order: str,
        evolution_time: float,
        num_steps: int,
        num_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Encode Hamiltonian as Trotterized time evolution.

        Constructs U = (e^(-iH_1 dt) ... e^(-iH_N dt))^n where:
        - H = sum_i H_i is the Hamiltonian decomposition
        - dt = evolution_time / num_steps
        - n = num_steps

        Args:
            pauli_strings: Dict mapping sparse Pauli tuples to coefficients.
                Format: {((qubit_idx, 'X'/'Y'/'Z'), ...): coefficient}
                Example: {((0, 'X'), (1, 'Z')): 0.5} represents 0.5*X_0*Z_1
            trotter_order: Product formula order:
                "first order", "second order", "fourth order", etc.
            evolution_time: Total evolution time t in exp(-iHt)
            num_steps: Number of Trotter steps
            num_qubits: System size
            **kwargs: Backend-specific options, may include:
                - combine_terms (bool): Group commuting terms
                - ordering_method (str): Term ordering strategy
                - tensor_contraction_method (str): Internal optimization

        Returns:
            Unitary operator implementing Trotterized evolution

        Raises:
            UnsupportedOperationError: If backend doesn't support Trotterization
            ValueError: If parameters are invalid (e.g., unknown trotter_order)
        """
        ...

    def encode_pauli_lcu(
        self,
        pauli_strings: Dict[tuple, float],
        num_qubits: int,
        prepare_type: str = 'AS',
        probability_eps: float = 0.002,
        **kwargs
    ) -> 'Unitary':
        """Encode Hamiltonian as LCU (Linear Combination of Unitaries) block encoding.

        Constructs block encoding of H = sum_i alpha_i U_i where:
        - U_i are Pauli string unitaries
        - alpha_i are coefficients
        - Uses SELECT and PREPARE oracles

        Args:
            pauli_strings: Dict mapping sparse Pauli tuples to coefficients
            num_qubits: System size
            prepare_type: State preparation method:
                'AS' - Alias sampling (common choice)
                'rotation' - Rotation-based preparation
                Backend-specific types may be supported
            probability_eps: Probability error tolerance for state preparation
            **kwargs: Backend-specific options

        Returns:
            Unitary implementing LCU block encoding

        Raises:
            UnsupportedOperationError: If backend doesn't support LCU
            ValueError: If parameters are invalid
        """
        ...

    def build_textbook_qpe(
        self,
        unitary: 'Unitary',
        num_phase_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Build textbook phase estimation circuit.

        Standard QPE algorithm:
        1. Prepare phase register in |+>^m state
        2. Apply controlled-U^(2^k) operations
        3. Apply inverse QFT

        Args:
            unitary: Time evolution or quantum walk operator to analyze
            num_phase_qubits: Precision bits (m in the literature)
            **kwargs: Backend-specific options, may include:
                - qft_implementation: Choice of QFT algorithm
                - ancilla_prep: Alternative to Hadamards

        Returns:
            Complete QPE circuit as Unitary

        Raises:
            UnsupportedOperationError: If backend doesn't support textbook QPE
            TypeError: If unitary is from incompatible backend
        """
        ...

    def build_qubitized_qpe(
        self,
        block_encoding: 'Unitary',
        num_phase_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Build qubitized phase estimation circuit.

        Qubitization-based QPE for block-encoded Hamiltonians.
        Uses quantum walk operator constructed from SELECT and PREPARE.

        Args:
            block_encoding: LCU or other block encoding from encode_pauli_lcu
            num_phase_qubits: Precision bits
            **kwargs: Backend-specific options

        Returns:
            Qubitized QPE circuit

        Raises:
            UnsupportedOperationError: If backend doesn't support qubitized QPE
            TypeError: If block_encoding is not an LCU-type encoding
        """
        ...
