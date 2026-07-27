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

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------
# Matrix computation
# -------------------------------------------------------------------------------------------------

# TODO: Why is this prefixed with an underscore as if it's intended to be private?
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

# TODO: Why is this prefixed with an underscore as if it's intended to be private?
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

# TODO: Should this be in file_io.py?
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
                timestep=approx_op.tevol_hbar
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
