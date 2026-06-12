"""
Matrix operations for exact Hamiltonian computation.

This module provides utilities for converting Pauli string representations
of Hamiltonians into matrix form, supporting both dense (small systems)
and matrix-free (large systems) approaches.

The choice between dense and sparse representation is based on a memory threshold:
    - Dense matrices require (2^num_qubits)^2 * 16 bytes (complex128)
    - Default threshold: 16 GB (allows dense up to 15 qubits, dimension 32768)
    - Calculations performed in log space to safely handle large qubit counts

Dense representation (memory ≤ threshold):
    - Use dense matrix representation via pauli_dict_to_matrix()
    - Fast, simple, allows direct matrix operations
    - Full matrix stored in memory

Sparse/matrix-free representation (memory > threshold):
    - Use PauliStringOperator class for matrix-free operations
    - Only materializes matrix-vector products, not full matrix
    - Compatible with scipy.sparse.linalg iterative methods
    - Scales to 25-30+ qubits
"""

import numpy as np
from scipy.sparse.linalg import LinearOperator
import logging

logger = logging.getLogger(__name__)


# =================================================================================================
# Pauli Matrix Definitions
# =================================================================================================

# 2×2 Pauli matrices
PAULI_I = np.array([[1, 0], [0, 1]], dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)

PAULI_MATRICES = {
    'I': PAULI_I,
    'X': PAULI_X,
    'Y': PAULI_Y,
    'Z': PAULI_Z
}


# =================================================================================================
# Dense Matrix Operations (Small Systems)
# =================================================================================================

def pauli_char_to_matrix(char):
    """
    Convert a single Pauli character to its 2×2 matrix representation.

    Parameters:
        char: str, one of 'I', 'X', 'Y', 'Z'

    Returns:
        2×2 complex numpy array

    Raises:
        ValueError: If char is not a valid Pauli character
    """
    char = char.upper()
    if char not in PAULI_MATRICES:
        raise ValueError(f"Invalid Pauli character: '{char}'. Must be one of I, X, Y, Z.")
    return PAULI_MATRICES[char]


def pauli_string_to_matrix(pauli_string):
    """
    Convert a Pauli string to its full matrix representation.

    A Pauli string like "XYZ" represents X ⊗ Y ⊗ Z (tensor product).

    Parameters:
        pauli_string: str, string of Pauli characters (e.g., "XYZ", "IXZI")

    Returns:
        2^n × 2^n complex numpy array, where n = len(pauli_string)

    Example:
        >>> mat = pauli_string_to_matrix("XY")
        >>> mat.shape
        (4, 4)
    """
    if not pauli_string:
        return np.array([[1]], dtype=complex)

    # Start with the first Pauli matrix
    result = pauli_char_to_matrix(pauli_string[0])

    # Tensor product with remaining matrices
    for char in pauli_string[1:]:
        result = np.kron(result, pauli_char_to_matrix(char))

    return result


def pauli_dict_to_matrix(pauli_dict, num_qubits):
    """
    Convert a dictionary of Pauli strings to a dense Hamiltonian matrix.

    Computes H = Σ_i coefficient_i * P_i, where P_i are Pauli string matrices.

    Parameters:
        pauli_dict: dict, mapping Pauli strings to coefficients
                   e.g., {"XYZ": 0.5, "ZZI": -0.3, "III": 1.0}
        num_qubits: int, number of qubits (determines matrix dimension 2^n)

    Returns:
        2^num_qubits × 2^num_qubits complex numpy array

    Example:
        >>> pauli_dict = {"ZZ": 0.5, "XX": -0.5, "II": 1.0}
        >>> H = pauli_dict_to_matrix(pauli_dict, 2)
        >>> H.shape
        (4, 4)
    """
    dimension = 2 ** num_qubits
    H = np.zeros((dimension, dimension), dtype=complex)

    for pauli_string, coefficient in pauli_dict.items():
        # Validate Pauli string length
        if len(pauli_string) != num_qubits:
            raise ValueError(
                f"Pauli string '{pauli_string}' has length {len(pauli_string)}, "
                f"but num_qubits={num_qubits}"
            )

        # Convert to matrix and add to Hamiltonian
        pauli_matrix = pauli_string_to_matrix(pauli_string)
        H += coefficient * pauli_matrix

    return H


# =================================================================================================
# Matrix-Free Linear Operator (Large Systems)
# =================================================================================================

class PauliStringOperator(LinearOperator):
    """
    Matrix-free linear operator for applying Pauli string Hamiltonians to state vectors.

    This class stores the Hamiltonian as a list of (coefficient, Pauli_string) pairs
    and implements matrix-vector multiplication without ever materializing the full matrix.

    Compatible with scipy.sparse.linalg iterative methods like eigsh(), which only
    need matvec() operations.

    Memory: O(num_terms + 2^n) instead of O(2^(2n)) for dense matrix.

    Parameters:
        pauli_dict: dict, mapping Pauli strings to coefficients
        num_qubits: int, number of qubits

    Attributes:
        shape: tuple, (dimension, dimension) where dimension = 2^num_qubits
        dtype: numpy dtype, always complex128

    Example:
        >>> pauli_dict = {"ZZ": 0.5, "XX": -0.5}
        >>> H_op = PauliStringOperator(pauli_dict, 2)
        >>> state = np.array([1, 0, 0, 0], dtype=complex)
        >>> result = H_op @ state  # Calls matvec internally
    """

    def __init__(self, pauli_dict, num_qubits):
        """Initialize the Pauli string operator."""
        self.pauli_dict = pauli_dict
        self.num_qubits = num_qubits
        self.dimension = 2 ** num_qubits

        # Validate all Pauli strings
        for pauli_string in pauli_dict.keys():
            if len(pauli_string) != num_qubits:
                raise ValueError(
                    f"Pauli string '{pauli_string}' has length {len(pauli_string)}, "
                    f"but num_qubits={num_qubits}"
                )

        # Initialize parent LinearOperator
        super().__init__(
            shape=(self.dimension, self.dimension),
            dtype=np.complex128
        )

    def _matvec(self, v):
        """
        Apply the Hamiltonian to a state vector: H|v⟩.

        Parameters:
            v: 1D numpy array of length 2^num_qubits

        Returns:
            1D numpy array of length 2^num_qubits
        """
        if len(v) != self.dimension:
            raise ValueError(
                f"Input vector has length {len(v)}, expected {self.dimension}"
            )

        result = np.zeros(self.dimension, dtype=complex)

        for pauli_string, coefficient in self.pauli_dict.items():
            # Apply each Pauli string term to the state
            result += coefficient * self._apply_pauli_string(pauli_string, v)

        return result

    def _apply_pauli_string(self, pauli_string, v):
        """
        Apply a single Pauli string operator to a state vector.

        This implements the action of a tensor product of Pauli operators
        without materializing the full matrix.

        Parameters:
            pauli_string: str, e.g., "XYZ"
            v: 1D numpy array, state vector

        Returns:
            1D numpy array, result of applying Pauli string
        """
        result = v.copy()

        for qubit_idx, pauli_char in enumerate(pauli_string):
            result = self._apply_single_qubit_pauli(result, qubit_idx, pauli_char)

        return result

    def _apply_single_qubit_pauli(self, v, qubit_idx, pauli_char):
        """
        Apply a single-qubit Pauli operator to the specified qubit.

        Parameters:
            v: 1D numpy array, state vector
            qubit_idx: int, which qubit to apply to (0 = first qubit)
            pauli_char: str, one of 'I', 'X', 'Y', 'Z'

        Returns:
            1D numpy array, result after applying Pauli operator
        """
        pauli_char = pauli_char.upper()

        if pauli_char == 'I':
            # Identity: no change
            return v

        # Create result array
        result = np.zeros_like(v)

        # Bit mask for the target qubit
        bit_mask = 1 << qubit_idx

        if pauli_char == 'X':
            # X flips the qubit: |0⟩ ↔ |1⟩
            for i in range(len(v)):
                j = i ^ bit_mask  # Flip the qubit bit
                result[i] = v[j]

        elif pauli_char == 'Y':
            # Y = -iXZ: flips and adds phase
            # Y|0⟩ = -i|1⟩, Y|1⟩ = i|0⟩
            for i in range(len(v)):
                j = i ^ bit_mask  # Flip the qubit bit
                # Phase depends on whether qubit was |0⟩ or |1⟩ BEFORE flip
                if j & bit_mask:  # Original (before flip) was |1⟩, result is |0⟩
                    result[i] = 1j * v[j]
                else:  # Original (before flip) was |0⟩, result is |1⟩
                    result[i] = -1j * v[j]

        elif pauli_char == 'Z':
            # Z adds phase: |0⟩ → |0⟩, |1⟩ → -|1⟩
            for i in range(len(v)):
                if i & bit_mask:  # Qubit is |1⟩
                    result[i] = -v[i]
                else:  # Qubit is |0⟩
                    result[i] = v[i]

        else:
            raise ValueError(f"Invalid Pauli character: '{pauli_char}'")

        return result

    def _rmatvec(self, v):
        """
        Apply the adjoint (conjugate transpose) of the Hamiltonian.

        For Hermitian Hamiltonians, H† = H, so this is the same as _matvec.
        """
        # For Hermitian operators, adjoint is the same as forward
        # But we need to conjugate the coefficients
        result = np.zeros(self.dimension, dtype=complex)

        for pauli_string, coefficient in self.pauli_dict.items():
            result += np.conj(coefficient) * self._apply_pauli_string(pauli_string, v)

        return result


# =================================================================================================
# Utility Functions
# =================================================================================================

def get_operator_type(num_qubits, memory_threshold_gb, force_dense=None, force_sparse=None):
    """
    Determine whether to use dense matrix or matrix-free operator.

    Parameters:
        num_qubits: int, number of qubits
        memory_threshold_gb: float, memory threshold in GB for dense representation
                           Dense matrix memory = (2^num_qubits)^2 * 16 bytes
        force_dense: bool or None, force dense representation
        force_sparse: bool or None, force sparse/matrix-free representation

    Returns:
        str, either "dense" or "sparse"

    Raises:
        ValueError: If both force_dense and force_sparse are True
    """
    if force_dense and force_sparse:
        raise ValueError("Cannot force both dense and sparse representations")

    if force_dense:
        return "dense"
    if force_sparse:
        return "sparse"

    # Automatic selection based on memory requirements
    # Dense matrix requires: dimension^2 * 16 bytes (complex128)
    # where dimension = 2^num_qubits
    #
    # Work in log space to avoid overflow:
    # log2(threshold_bytes) = log2(threshold_gb * 2^30) = log2(threshold_gb) + 30
    # log2(required_bytes) = log2(2^(2*num_qubits) * 16) = 2*num_qubits + 4
    log2_threshold_bytes = np.log2(memory_threshold_gb) + 30
    log2_required_bytes = 2 * num_qubits + 4
    memory_gb = 2 ** (log2_required_bytes - 30)
    if log2_required_bytes <= log2_threshold_bytes:
        # Compute actual memory for logging
        logger.verbose(
            f"Using dense matrix representation for {num_qubits} qubits "
            f"(dim={2**num_qubits}, memory={memory_gb:.3f} GB)"
        )
        return "dense"
    else:
        # Compute what memory would be required
        logger.verbose(
            f"Using sparse/matrix-free representation for {num_qubits} qubits "
            f"(would require {memory_gb:.3f} GB)"
        )
        return "sparse"


def create_hamiltonian_operator(pauli_dict, num_qubits, memory_threshold_gb, force_dense=None, force_sparse=None):
    """
    Create either a dense matrix or matrix-free operator for the Hamiltonian.

    Parameters:
        pauli_dict: dict, mapping Pauli strings to coefficients
        num_qubits: int, number of qubits
        memory_threshold_gb: float, memory threshold in GB for dense representation
        force_dense: bool or None, force dense representation
        force_sparse: bool or None, force sparse/matrix-free representation

    Returns:
        numpy array (dense) or PauliStringOperator (sparse)
    """
    operator_type = get_operator_type(num_qubits, memory_threshold_gb, force_dense, force_sparse)

    if operator_type == "dense":
        return pauli_dict_to_matrix(pauli_dict, num_qubits)
    else:
        return PauliStringOperator(pauli_dict, num_qubits)
