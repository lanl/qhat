"""
Matrix and eigendecomposition computation, conversion, and output operations.

This module handles:
- Matrix computation (exact Hamiltonian, approximate unitary)
- Eigendecomposition computation and processing
- Conversion between unitary eigenvalues and eigenenergies
- Matrix and eigendecomposition output to files
- Flexible operator output system using OperatorRepresentation
"""

import numpy as np
import logging
from pathlib import Path
import datetime

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import save_matrix, save_eigendecomposition
from qhat.analysis.matrix_operations import compute_unitarity_error
from qhat.analysis.operators import convert_unitary_eigenvalues_to_eigenenergies

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------
# Matrix computation
# -------------------------------------------------------------------------------------------------

# TODO: Do I need this function?
def _compute_unitary_matrix(algorithm):
    """
    Compute the unitary matrix representation of the algorithm.

    Parameters:
        algorithm: The algorithm bloq to analyze

    Returns:
        The unitary matrix as a numpy array

    Raises:
        AttributeError: If algorithm doesn't have tensor_contract() method
        Exception: If tensor_contract() fails
    """
    if not hasattr(algorithm, 'tensor_contract'):
        raise AttributeError(
            f"Algorithm of type {type(algorithm).__name__} does not have a "
            "'tensor_contract()' method. Cannot generate unitary matrix."
        )

    logger.verbose("Computing unitary matrix via tensor contraction...")
    try:
        return algorithm.tensor_contract()
    except Exception as e:
        logger.info(
            f"ERROR: Failed to compute unitary matrix: {e}\n"
            "This may indicate a bug in the algorithm's tensor_contract() implementation."
        )
        raise

# -------------------------------------------------------------------------------------------------

# TODO: Do I need this function?
def _compute_exact_matrix(hamiltonian, config_analysis):
    """
    Compute the exact matrix representation of the Hamiltonian.

    This computes the Hamiltonian matrix without any approximations (no Trotter,
    no double-factorization). The choice between dense and sparse/matrix-free
    representation is based on the memory threshold in config_analysis.

    Parameters:
        hamiltonian: The Hamiltonian object
        config_analysis: Analysis configuration with matrix_memory_threshold_gb

    Returns:
        Dense numpy array (small systems) or PauliStringOperator (large systems)

    Raises:
        Exception: If matrix computation fails
    """
    logger.verbose("Computing exact Hamiltonian matrix...")
    try:
        return hamiltonian.to_matrix(
            memory_threshold_gb=config_analysis.matrix_memory_threshold_gb
        )
    except Exception as e:
        logger.info(
            f"ERROR: Failed to compute exact Hamiltonian matrix: {e}\n"
            "This may indicate an issue with the Hamiltonian representation."
        )
        raise

# -------------------------------------------------------------------------------------------------
# Matrix output
# -------------------------------------------------------------------------------------------------

def output_unitary_matrix(filename, unitary_matrix) -> dict:
    """
    Save unitary matrix to a file

    Parameters:
        filename: Name of file to save matrix to
        unitary_matrix: The unitary matrix to save

    Returns:
        Dictionary with matrix metadata: shape, file, format, unitarity_error, norm
    """

    # Log basic properties
    logger.verbose(f"Matrix shape: {unitary_matrix.shape}")
    logger.verbose(f"Matrix dtype: {unitary_matrix.dtype}")

    # Compute unitarity check: ||U†U - I||_F
    try:
        matrix_norm = np.linalg.norm(unitary_matrix, ord='fro')
        unitarity_error = compute_unitarity_error(unitary_matrix)
        logger.verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        logger.verbose(f"Unitarity error ||U†U - I||_F: {unitarity_error:.6e}")
    except Exception as e:
        logger.info(f"WARNING: Could not compute unitarity check: {e}")
        matrix_norm = None
        unitarity_error = None

    # Save matrix to file (format auto-detected from extension)
    output_file = filename
    save_matrix(
        output_file, unitary_matrix,
        unitarity_error=unitarity_error,
        matrix_norm=matrix_norm
    )

    # Return metadata
    return {
        'matrix_shape': unitary_matrix.shape,
        'matrix_file': str(output_file),
        'matrix_format': Path(output_file).suffix,
        'unitarity_error': float(unitarity_error) if unitarity_error is not None else None,
        'matrix_norm': float(matrix_norm) if matrix_norm is not None else None
    }

# -------------------------------------------------------------------------------------------------
# Eigendecomposition computation
# -------------------------------------------------------------------------------------------------

def _eigendecompose_full(matrix, matrix_type):
    """
    Perform full eigendecomposition of a matrix.

    For exact (Hermitian) matrices: Uses eigh for efficiency (real eigenvalues).
    For approximate (unitary) matrices: Uses eig to get complex eigenvalues on unit circle.

    Args:
        matrix: Matrix to decompose (dense numpy array)
        matrix_type: String describing the matrix type ('exact' or 'approximate', for logging)

    Returns:
        tuple: (eigenvalues, eigenvectors)
            eigenvalues: All eigenvalues (unsorted, as returned by scipy)
            eigenvectors: All eigenvectors (columns correspond to eigenvalues)

    Raises:
        ValueError: If matrix is too large (matrix-free operators not supported for full decomposition)
    """
    import scipy.linalg
    from qhat.analysis.matrix_operations import PauliStringOperator

    dimension = matrix.shape[0]
    is_matrix_free = isinstance(matrix, PauliStringOperator) or hasattr(matrix, 'matvec')

    if is_matrix_free:
        raise ValueError(
            f"Full eigendecomposition not supported for matrix-free operators. "
            f"Matrix dimension {dimension} is too large for dense eigendecomposition."
        )

    logger.info(f"Computing full eigendecomposition for {matrix_type} matrix (dimension {dimension})")

    # Use appropriate eigendecomposition based on matrix type
    if matrix_type == 'exact':
        # Exact Hamiltonian is Hermitian: use eigh for real eigenvalues
        eigenvalues, eigenvectors = scipy.linalg.eigh(matrix)
        logger.verbose(f"  Eigenvalue range (unsorted): [{eigenvalues.min():.6e}, {eigenvalues.max():.6e}]")
    else:  # matrix_type == 'approximate'
        # Approximate unitary is not Hermitian: use eig for complex eigenvalues
        eigenvalues, eigenvectors = scipy.linalg.eig(matrix)
        # For unitary matrices, all eigenvalues should have magnitude 1
        magnitudes = np.abs(eigenvalues)
        logger.verbose(f"  Eigenvalue magnitude range: [{magnitudes.min():.6e}, {magnitudes.max():.6e}]")

    return eigenvalues, eigenvectors

# -------------------------------------------------------------------------------------------------

def _process_eigendecomposition(matrix, matrix_type, timestep=None, energy_shift=0.0):
    """
    Process full eigendecomposition: compute, convert (if unitary), sort, save, and return data.

    Args:
        matrix: Matrix to decompose (dense array, matrix-free not supported)
        matrix_type: String describing the matrix type ('exact' or 'approximate')
        timestep: Required for matrix_type='approximate' to convert phases to eigenenergies
        energy_shift: Energy shift to add back to eigenenergies (for comparison on same scale)

    Returns:
        Dictionary with eigenenergies, eigenvectors, file path, and metadata

    Raises:
        ValueError: If matrix is None or if timestep not provided for approximate
    """
    logger.info(f"Computing {matrix_type} matrix eigendecomposition (full spectrum)")

    if matrix is None:
        raise ValueError(
            f"{matrix_type}_matrix is required for {matrix_type} eigendecomposition. "
            "Compute the matrix before calling eigendecomposition_analysis()."
        )

    # Compute full eigendecomposition
    eigenvalues_raw, eigenvectors_raw = _eigendecompose_full(matrix, matrix_type)

    # Convert to eigenenergies and sort
    if matrix_type == 'approximate':
        if timestep is None:
            raise ValueError("timestep is required for approximate eigendecomposition to convert phases to eigenenergies")

        # Convert unitary eigenvalues to eigenenergies
        eigenenergies, eigenphases = convert_unitary_eigenvalues_to_eigenenergies(
            eigenvalues_raw, timestep
        )

        # Apply energy shift correction to restore original energy scale
        eigenenergies_corrected = eigenenergies + energy_shift

        # Sort by eigenenergy
        sort_indices = np.argsort(eigenenergies_corrected)
        eigenenergies_sorted = eigenenergies_corrected[sort_indices]
        eigenvectors_sorted = eigenvectors_raw[:, sort_indices]

        logger.info(f"Converted {len(eigenphases)} unitary eigenvalues to eigenenergies")
        if abs(energy_shift) > 1e-10:
            logger.info(f"  Applied energy shift correction: {energy_shift:.6e}")
        logger.info(f"  Eigenenergy range (sorted): [{eigenenergies_sorted[0]:.6e}, {eigenenergies_sorted[-1]:.6e}]")

        # Save sorted results with optional debugging data
        output_file = f'{matrix_type}_eigendecomposition.npz'
        save_eigendecomposition(
            output_file,
            eigenenergies=eigenenergies_sorted,
            eigenvectors=eigenvectors_sorted,
            matrix_type=matrix_type,
            timestep=timestep,
            unitary_eigenvalues=eigenvalues_raw,
            eigenphases=eigenphases
        )

    else:  # matrix_type == 'exact'
        # Eigenvalues are already eigenenergies, apply energy shift correction and sort
        eigenenergies_raw = np.real(eigenvalues_raw)  # Should be real for Hermitian

        # Apply energy shift correction to restore original energy scale
        eigenenergies_corrected = eigenenergies_raw + energy_shift

        sort_indices = np.argsort(eigenenergies_corrected)
        eigenenergies_sorted = eigenenergies_corrected[sort_indices]
        eigenvectors_sorted = eigenvectors_raw[:, sort_indices]

        if abs(energy_shift) > 1e-10:
            logger.info(f"Applied energy shift correction: {energy_shift:.6e}")
        logger.info(f"Eigenenergy range (sorted): [{eigenenergies_sorted[0]:.6e}, {eigenenergies_sorted[-1]:.6e}]")

        # Save sorted results
        output_file = f'{matrix_type}_eigendecomposition.npz'
        save_eigendecomposition(
            output_file,
            eigenenergies=eigenenergies_sorted,
            eigenvectors=eigenvectors_sorted,
            matrix_type=matrix_type
        )

    return {
        'file': output_file,
        'eigenenergies': eigenenergies_sorted,
        'num_eigenstates': len(eigenenergies_sorted),
        'eigenenergy_range': [float(eigenenergies_sorted[0]), float(eigenenergies_sorted[-1])]
    }

# -------------------------------------------------------------------------------------------------
# Eigendecomposition analysis
# -------------------------------------------------------------------------------------------------

# eigendecomposition analysis
# - selectively dispatches to _process_eigendecomposition (the selection appears to be kind of a
#   mess, so I'm not sure it's reliable and I haven't actually analyzed it in detail)
# - _process_eigendecomposition
#   - accepts a matrix
#   - eigendecomposes the matrix
#   - knows that approx = shifted time-evolution, exact = shifted Hamiltonian
#   - converts raw eigenvalues to unshifted eigenenergies
#   - sorts eigendecomposition by eigenenergies
#   - write (exact|approximate)_eigendecomposition.npz (with metadata)
#   - returns some summary metadata (including the list of eigenvalues)
# - generates *_eigendecomposition.npz, returns metadata including list of eigenvalues

def eigendecomposition_analysis(
        config_analysis: AnalysisConfiguration,
        timestep=None,
        energy_shift=0.0,
        exact_matrix=None,
        unitary_matrix=None,
        requires_exact_eigendecomposition_func=None,
        requires_approximate_eigendecomposition_func=None) -> dict:
    """
    Compute full eigendecompositions of exact and/or approximate matrices.

    This is the single decision point for computing eigendecompositions.
    Always computes full spectrum (all eigenstates), sorted by eigenenergy.

    Parameters:
        config_analysis: Analysis configuration with eigendecomposition settings
        timestep: Time evolution parameter (required for approximate eigendecomposition)
        energy_shift: Energy shift applied to Hamiltonian (added back to eigenenergies for comparison)
        exact_matrix: Pre-computed exact matrix (required if eigendecomposition_matrices is 'exact' or 'both')
        unitary_matrix: Pre-computed unitary matrix (required if eigendecomposition_matrices is 'approximate' or 'both')
        requires_exact_eigendecomposition_func: Function to check if exact eigendecomposition is needed
        requires_approximate_eigendecomposition_func: Function to check if approximate eigendecomposition is needed

    Returns:
        Dictionary with eigendecomposition data (eigenenergies, eigenvectors) and metadata

    Raises:
        ValueError: If required matrices or timestep are not provided
    """
    # Determine which eigendecompositions need to be computed
    # Use requires_* functions to check if eigendecompositions are needed
    # (handles both explicit requests and implicit needs like eigenvalue error analysis)
    if requires_exact_eigendecomposition_func is None:
        # Default behavior: check if matrices are provided
        need_exact = exact_matrix is not None
    else:
        need_exact = requires_exact_eigendecomposition_func(config_analysis)

    if requires_approximate_eigendecomposition_func is None:
        # Default behavior: check if matrices are provided
        need_approx = unitary_matrix is not None
    else:
        need_approx = requires_approximate_eigendecomposition_func(config_analysis)

    if not need_exact and not need_approx:
        logger.info("Eigendecomposition analysis not requested")
        return {}

    logger.info(f"Starting eigendecomposition analysis")
    if config_analysis.enable_eigenvalue_errors:
        logger.verbose(f"  enable_eigenvalue_errors: True (requires both eigendecompositions)")

    results = {}

    # Compute exact eigendecomposition if requested
    if need_exact:
        results['exact_eigendecomposition'] = _process_eigendecomposition(
            exact_matrix, 'exact', timestep=None, energy_shift=energy_shift
        )

    # Compute approximate eigendecomposition if requested
    if need_approx:
        results['approximate_eigendecomposition'] = _process_eigendecomposition(
            unitary_matrix, 'approximate', timestep=timestep, energy_shift=energy_shift
        )

    return results

# -------------------------------------------------------------------------------------------------
# Flexible operator output system
# -------------------------------------------------------------------------------------------------

def save_requested_operator_outputs(
        operator_output_requests,
        exact_op,
        approx_op) -> dict:
    """
    Save all requested operator forms using OperatorRepresentation.

    Returns
    -------
    dict
        Information about saved files
    """
    from qhat.analysis.operators import OperatorRepresentation

    results = {
        'matrix_outputs': [],
        'eigendecomposition_outputs': []
    }

    # Check if any requests exist
    if not operator_output_requests:
        return results

    # Create OperatorRepresentation wrappers if not provided
    if exact_op is None:
        raise ValueError("Missing exact_op in save_requested_operator_outputs.")
    if approx_op is None:
        raise ValueError("Missing approx_op in save_requested_operator_outputs.")

    operators = {
        'exact': exact_op,
        'approximate': approx_op
    }

    # Process operator output requests
    for request in operator_output_requests:
        op = operators[request['source']]
        energy_shifted = request['energy_shifted']
        representation = request['representation']

        shift_str = 'shifted' if energy_shifted else 'unshifted'
        logger.info(
            f"Saving {request['source']} {request['operator_type']} "
            f"({shift_str}) as {representation} to {request['filename']}"
        )

        if representation == 'matrix':
            # Get and save matrix representation
            matrix = op.get(
                operator_type=request['operator_type'],
                energy_shifted=energy_shifted,
                representation='dense_matrix'
            )

            save_matrix(request['filename'], matrix)

            results['matrix_outputs'].append({
                'filename': request['filename'],
                'source': request['source'],
                'operator_type': request['operator_type'],
                'energy_shifted': energy_shifted,
                'shape': matrix.shape
            })

        elif representation == 'eigendecomposition':
            # Get and save eigendecomposition representation
            eigendata = op.get(
                operator_type=request['operator_type'],
                energy_shifted=energy_shifted,
                representation='eigendecomposition'
            )

            # Sort by eigenvalues (ascending)
            sort_indices = np.argsort(eigendata['eigenvalues'].real)
            eigenvalues_sorted = eigendata['eigenvalues'][sort_indices]
            eigenvectors_sorted = eigendata['eigenvectors'][:, sort_indices]

            save_eigendecomposition(
                request['filename'],
                eigenenergies=eigenvalues_sorted,
                eigenvectors=eigenvectors_sorted,
                matrix_type=f"{request['source']}_{request['operator_type']}_{shift_str}",
                timestep=approx_op.timestep
            )

            results['eigendecomposition_outputs'].append({
                'filename': request['filename'],
                'source': request['source'],
                'operator_type': request['operator_type'],
                'energy_shifted': energy_shifted,
                'num_eigenvalues': len(eigenvalues_sorted),
                'eigenvalue_range': [float(eigenvalues_sorted[0].real), float(eigenvalues_sorted[-1].real)]
            })

    return results
