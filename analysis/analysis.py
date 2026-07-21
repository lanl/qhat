import numpy as np
from datetime import datetime
import logging
import os
from pathlib import Path

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount
from qualtran.cirq_interop.t_complexity_protocol import t_complexity

from qhat.analysis.config_types import AnalysisConfiguration, GeneralConfiguration
from qhat.analysis.file_io import save_matrix, load_state, save_state

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def convert_unitary_eigenvalues_to_eigenenergies(unitary_eigenvalues, timestep, hbar=1.0):
    """
    Convert unitary eigenvalues e^(-iφ) to Hamiltonian eigenenergies E.

    The time evolution operator is U = exp(-iHt/ℏ), so if H has eigenenergy E,
    then U has eigenvalue exp(-iEt/ℏ) = exp(-iφ) where φ = Et/ℏ is the eigenphase.

    Assumption: The Hamiltonian has been shifted and scaled (by existing code) such that
    all eigenenergies produce eigenphases in the range [0, 2π), preventing aliasing.

    Args:
        unitary_eigenvalues: Complex eigenvalues of unitary U = exp(-iHt/ℏ) (on unit circle)
        timestep: Time evolution parameter t (units: ℏ/energy, e.g., ℏ/Hartree)
        hbar: Value of ℏ (default: 1.0 for atomic units)

    Returns:
        tuple: (eigenenergies, eigenphases)
            eigenenergies: Real eigenvalues of Hamiltonian (shifted/scaled, units: energy)
            eigenphases: Phases φ ∈ [0, 2π) where φ = Et/ℏ
    """
    # Extract phases using np.angle, which returns [-π, π]
    # Since U = exp(-iφ), we have arg(U) = -φ (mod 2π)
    phases_neg_pi_to_pi = np.angle(unitary_eigenvalues)

    # Eigenphase φ = -arg(U), then map to [0, 2π) convention
    eigenphases_raw = -phases_neg_pi_to_pi
    eigenphases = np.where(eigenphases_raw < 0,
                           eigenphases_raw + 2*np.pi,
                           eigenphases_raw)

    # Convert to eigenenergies: E = φ * ℏ / t
    # These are shifted/scaled eigenenergies that correspond to the phases in [0, 2π)
    eigenenergies = eigenphases * hbar / timestep

    return eigenenergies, eigenphases

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    logger.verbose("Estimating resources with pyLIQTR.")

    # TODO: rotation error
    #       -- argument rotation_gate_precision sets the precision for a single rotation gate
    #       -- argument algorithm_precision sets the precision for the whole algorithm (i.e., it
    #          sets rotation_gate_precision to algorithm_precision / number of rotations)
    # TODO: profile?
    #       -- argument profile = True: keep rotations as a separate count
    #       -- argument profile = False: estimate rotations as Clifford+T
    resources = estimate_pyliqtr(algorithm)

    resource_dict = {
        "Clifford_count" : resources["Clifford"],
        "T_count"        : resources["T"],
        }
    if "LogicalQubits" in resources:
        resource_dict["qubit_count"] = resources["LogicalQubits"]
    else:
        get_cost_value(algorithm, QubitCount())
    return resource_dict

# -------------------------------------------------------------------------------------------------

def resource_estimation_qualtran(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    logger.verbose("Estimating resources with Qualtran.")

    t_cliff_resources = t_complexity(algorithm)
    qubits = get_cost_value(algorithm, QubitCount())

    # TODO:  add rotation count to this?
    resource_dict = {
        "Clifford_count" : t_cliff_resources.clifford,
        "T_count"        : t_cliff_resources.t,
        "qubit_count"    : qubits
        }
    return resource_dict

# -------------------------------------------------------------------------------------------------

def estimate_resources(
        config_analysis: AnalysisConfiguration,
        algorithm) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_analysis, algorithm)
    elif config_analysis.resource_estimator.lower() == "qualtran":
        return resource_estimation_qualtran(config_analysis, algorithm)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")

# -------------------------------------------------------------------------------------------------

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

def output_unitary_matrix(
        config_analysis: AnalysisConfiguration,
        algorithm,
        unitary_matrix) -> dict:
    """
    Generate and save the unitary matrix representation of the algorithm.

    Parameters:
        config_analysis: Analysis configuration with matrix_output_format and matrix_output_file
        algorithm: The algorithm bloq to analyze
        unitary_matrix: The unitary matrix to save (pre-computed)

    Returns:
        Dictionary with matrix metadata: shape, file, format, unitarity_error, norm
    """

    # Log basic properties
    logger.verbose(f"Matrix shape: {unitary_matrix.shape}")
    logger.verbose(f"Matrix dtype: {unitary_matrix.dtype}")

    # Compute unitarity check: ||U†U - I||_F
    try:
        matrix_norm = np.linalg.norm(unitary_matrix, ord='fro')
        U_dag_U = np.conj(unitary_matrix.T) @ unitary_matrix
        identity = np.eye(unitary_matrix.shape[0])
        unitarity_error = np.linalg.norm(U_dag_U - identity, ord='fro')
        logger.verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        logger.verbose(f"Unitarity error ||U†U - I||_F: {unitarity_error:.6e}")
    except Exception as e:
        logger.info(f"WARNING: Could not compute unitarity check: {e}")
        matrix_norm = None
        unitarity_error = None

    # Save matrix to file (format auto-detected from extension)
    output_file = config_analysis.matrix_output_file
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

def exact_matrix_output(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        exact_matrix) -> dict:
    """
    Save the exact Hamiltonian matrix representation.

    Parameters:
        config_analysis: Analysis configuration with exact_matrix_output_file
        hamiltonian: The Hamiltonian object
        exact_matrix: The exact matrix to save (pre-computed)

    Returns:
        Dictionary with matrix metadata: shape, file, format, hermiticity_error, norm

    Note:
        For large systems, exact_matrix may be a matrix-free operator rather than
        a dense array. In that case, certain properties (like saving to file) may
        not be supported or may require special handling.
    """
    from qhat.analysis.matrix_operations import PauliStringOperator

    # Check if this is a matrix-free operator
    is_matrix_free = isinstance(exact_matrix, PauliStringOperator)

    if is_matrix_free:
        logger.verbose(f"Matrix-free operator with shape: {exact_matrix.shape}")
        logger.info(
            "WARNING: Matrix-free operator cannot be directly saved to file. "
            "Skipping matrix output for large system."
        )
        return {
            'matrix_shape': exact_matrix.shape,
            'matrix_file': None,
            'matrix_format': None,
            'hermiticity_error': None,
            'matrix_norm': None,
            'matrix_free': True,
            'note': 'Matrix-free operator not saved (too large)'
        }

    # For dense matrices, proceed with normal output
    logger.verbose(f"Matrix shape: {exact_matrix.shape}")
    logger.verbose(f"Matrix dtype: {exact_matrix.dtype}")

    # Compute Hermiticity check: ||H - H†||_F
    try:
        matrix_norm = np.linalg.norm(exact_matrix, ord='fro')
        H_dag = np.conj(exact_matrix.T)
        hermiticity_error = np.linalg.norm(exact_matrix - H_dag, ord='fro')
        logger.verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        logger.verbose(f"Hermiticity error ||H - H†||_F: {hermiticity_error:.6e}")
    except Exception as e:
        logger.info(f"WARNING: Could not compute Hermiticity check: {e}")
        matrix_norm = None
        hermiticity_error = None

    # Save matrix to file (format auto-detected from extension)
    output_file = config_analysis.exact_matrix_output_file
    save_matrix(
        output_file, exact_matrix,
        hermiticity_error=hermiticity_error,
        matrix_norm=matrix_norm
    )

    # Return metadata
    return {
        'matrix_shape': exact_matrix.shape,
        'matrix_file': str(output_file),
        'matrix_format': Path(output_file).suffix,
        'hermiticity_error': float(hermiticity_error) if hermiticity_error is not None else None,
        'matrix_norm': float(matrix_norm) if matrix_norm is not None else None,
        'matrix_free': False
    }

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
    from qhat.analysis.file_io import save_eigendecomposition

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

def eigendecomposition_analysis(
        config_analysis: AnalysisConfiguration,
        timestep=None,
        energy_shift=0.0,
        exact_matrix=None,
        unitary_matrix=None) -> dict:
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

    Returns:
        Dictionary with eigendecomposition data (eigenenergies, eigenvectors) and metadata

    Raises:
        ValueError: If required matrices or timestep are not provided
    """
    which_matrices = config_analysis.eigendecomposition_matrices

    if which_matrices is None:
        logger.info("Eigendecomposition analysis not requested (eigendecomposition_matrices is None)")
        return {}

    logger.info(f"Starting eigendecomposition analysis")
    logger.verbose(f"  eigendecomposition_matrices: {which_matrices}")

    # Determine which matrices we need
    need_exact = which_matrices in ['exact', 'both']
    need_approx = which_matrices in ['approximate', 'both']

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

def error_analysis(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        algorithm,
        exact_matrix=None,
        unitary_matrix=None,
        exact_eigendecomp=None,
        approx_eigendecomp=None,
        timestep=None,
        energy_shift=0.0) -> dict:
    """
    Compute error metrics comparing exact and approximate representations.

    Three independent error types:
    1. Eigenvalue errors - Compares eigenenergies: λ(H_approx) vs λ(H_exact)
    2. Matrix norm errors - Compares time evolution operators: ||U_exact - U_approx||
    3. State-dependent errors - Compares time-evolved states: ||U_exact|ψ⟩ - U_approx|ψ⟩||

    Parameters:
        config_analysis: Analysis configuration with error analysis settings
        hamiltonian: Hamiltonian object
        algorithm: Algorithm bloq
        exact_matrix: Pre-computed exact matrix (optional)
        unitary_matrix: Pre-computed unitary matrix (optional)
        exact_eigendecomp: Pre-computed exact eigendecomposition (optional)
        approx_eigendecomp: Pre-computed approximate eigendecomposition (optional)
        timestep: Time evolution parameter t (required for matrix/state errors).
                  Used to compute U_exact = exp(-i * H_exact * t) from H_exact.
        energy_shift: Energy shift E applied to approximate Hamiltonian (default: 0.0).
                      Used to match global phases: U_exact_shifted = exp(i*E*t) * U_exact.

    Returns:
        Dictionary with error metrics
    """
    from qhat.analysis.matrix_operations import PauliStringOperator
    from qhat.analysis.file_io import load_eigendecomposition, load_state
    from qhat.analysis.operators import OperatorRepresentation
    import scipy.linalg

    logger.info("Starting error analysis")

    results = {}

    # =============================================================================================
    # 1. EIGENENERGY ERROR
    # =============================================================================================

    if config_analysis.enable_eigenvalue_errors:
        logger.info("Computing eigenenergy errors for all eigenstates")

        if exact_eigendecomp is None or approx_eigendecomp is None:
            raise ValueError(
                "Both eigendecompositions must be computed in order to compare eigenenergies. "
                "Ensure eigendecomposition_matrices is set to 'both' when enable_eigenvalue_errors is True."
            )

        # Get eigenenergies from both decompositions (both already sorted by energy)
        exact_eigenenergies = exact_eigendecomp['eigenenergies']
        approx_eigenenergies = approx_eigendecomp['eigenenergies']

        # Verify same dimension (should have all eigenvalues since full decomposition)
        if len(exact_eigenenergies) != len(approx_eigenenergies):
            raise ValueError(
                f"Dimension mismatch: exact has {len(exact_eigenenergies)} eigenstates, "
                f"approximate has {len(approx_eigenenergies)} eigenstates."
            )

        # Element-wise comparison (both already sorted by energy)
        absolute_errors = exact_eigenenergies - approx_eigenenergies
        relative_errors = absolute_errors / np.abs(exact_eigenenergies)

        num_eigenstates = len(exact_eigenenergies)
        logger.info(f"Computed errors for {num_eigenstates} eigenstates (sorted by energy)")
        logger.info(f"  Ground state: exact={exact_eigenenergies[0]:.6e}, approx={approx_eigenenergies[0]:.6e}, error={absolute_errors[0]:.6e}")
        logger.info(f"  Highest state: exact={exact_eigenenergies[-1]:.6e}, approx={approx_eigenenergies[-1]:.6e}, error={absolute_errors[-1]:.6e}")
        logger.info(f"  Max absolute error: {np.abs(absolute_errors).max():.6e}")
        logger.info(f"  Max relative error: {np.abs(relative_errors).max():.6e}")

        results['eigenenergy_errors'] = {
            'num_eigenstates': num_eigenstates,
            'absolute_errors': absolute_errors.tolist(),
            'relative_errors': relative_errors.tolist(),
            'max_absolute_error': float(np.abs(absolute_errors).max()),
            'max_relative_error': float(np.abs(relative_errors).max())
        }

    # =============================================================================================
    # OPERATOR CONVERSION: Wrap operators in OperatorRepresentation for clean handling
    # =============================================================================================
    # Available inputs:
    #   - exact_matrix: H_exact (unshifted Hamiltonian, dense matrix or PauliStringOperator)
    #   - unitary_matrix: U_s,approx (energy-shifted time-evolution operator)
    #
    # OperatorRepresentation provides unified interface for conversions:
    #   - H ↔ U (Hamiltonian ↔ time-evolution operator)
    #   - Shifted ↔ unshifted (energy shift application/removal)
    #   - Dense matrix ↔ eigendecomposition
    #
    # Phase 2: Dense matrices only; matrix-free support deferred to Phase 3

    # Wrap operators in OperatorRepresentation instances
    exact_op = None
    approx_op = None

    # Check if matrix/state errors are requested (these need unitary operators)
    needs_unitaries = (
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )

    if needs_unitaries:
        logger.info("Preparing operators for matrix/state error analysis")

        # Validate required inputs
        if exact_matrix is None:
            raise ValueError(
                "Matrix/state error analysis requires the exact Hamiltonian matrix. "
                "Ensure exact_matrix_output_file is set or eigendecomposition is enabled."
            )

        if unitary_matrix is None:
            raise ValueError(
                "Matrix/state error analysis requires the approximate unitary matrix. "
                "The algorithm must produce a unitary matrix representation."
            )

        if timestep is None:
            raise ValueError(
                "Matrix/state error analysis requires timestep parameter. "
                "Ensure the algorithm provides a timestep (e.g., time evolution algorithms)."
            )

        # Check if dense or matrix-free
        is_exact_dense = isinstance(exact_matrix, np.ndarray)
        is_approx_dense = isinstance(unitary_matrix, np.ndarray)

        if not (is_exact_dense and is_approx_dense):
            # Matrix-free case: not implemented in Phase 2
            raise NotImplementedError(
                "Matrix/state error analysis not yet implemented for matrix-free operators. "
                "Phase 2 supports only dense matrices (systems with n ≤ 15 qubits). "
                "Reduce system size to enable dense matrix computation, or disable these error types."
            )

        logger.verbose(f"Creating OperatorRepresentation instances")
        logger.verbose(f"  Timestep: t = {timestep}")
        logger.verbose(f"  Energy shift: E = {energy_shift}")

        # Wrap exact Hamiltonian in OperatorRepresentation
        # Note: exact_matrix is H' = H + E*I (shifted up by E to make eigenvalues positive)
        # unitary_matrix is U' = exp(-i*H'*t) (also uses the shifted Hamiltonian)
        # Both are on the shifted scale, and we want to unshift to the physical H for comparison.
        exact_op = OperatorRepresentation(
            data=exact_matrix,
            operator_type='hamiltonian',
            energy_shifted=True,  # Input IS shifted (H' = H + E*I)
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift  # Will unshift by subtracting E
        )
        logger.verbose(f"  Created exact operator representation (H', shifted)")

        # Wrap approximate time-evolution operator in OperatorRepresentation
        approx_op = OperatorRepresentation(
            data=unitary_matrix,
            operator_type='time_evolution',
            energy_shifted=True,  # Input IS shifted (U' from H')
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift  # Will unshift with phase factor
        )
        logger.verbose(f"  Created approx operator representation (U', shifted)")

        logger.info(f"Operator representations ready for conversion on demand")

    # =============================================================================================
    # 2. MATRIX NORM ERRORS: ||U_exact - U_approx||
    # =============================================================================================
    # Compares unshifted time evolution operators (both represent physical Hamiltonian).
    # Operators were converted in the section above.

    if config_analysis.error_matrix_norms is not None:
        # Note: error_matrix_norms is normalized to list during validation
        norms_to_compute = config_analysis.error_matrix_norms

        logger.info(f"Computing matrix norm errors: {norms_to_compute}")

        # Require matrices to be provided by caller
        if exact_matrix is None:
            raise ValueError(
                "Matrix norm error analysis requires the exact Hamiltonian matrix, but it was not computed. "
            )

        if unitary_matrix is None:
            raise ValueError(
                "Matrix norm error analysis requires the approximate/unitary matrix, but it was not computed. "
            )

        # Check if matrices are dense or matrix-free
        is_exact_dense = isinstance(exact_matrix, np.ndarray)
        is_approx_dense = isinstance(unitary_matrix, np.ndarray)
        is_dense = is_exact_dense and is_approx_dense

        dimension = exact_matrix.shape[0]
        num_qubits = int(np.log2(dimension))

        if is_dense:
            # Small systems: direct computation of ||U_exact - U_approx||
            logger.verbose(f"Using dense matrices for norm computation (dimension={dimension})")

            # Get unshifted time-evolution operators from OperatorRepresentation
            if exact_op is None or approx_op is None:
                raise RuntimeError(
                    "Operator representations should have been created but are None. "
                    "This is an internal error in the operator conversion logic."
                )

            logger.verbose(f"  Converting H' → U (unshifted, physical)")
            exact_unitary_matrix = exact_op.get(
                operator_type='time_evolution',
                energy_shifted=False,  # Get physical U from H
                representation='dense_matrix'
            )

            logger.verbose(f"  Converting U' → U (unshifted, physical)")
            approx_unitary_matrix = approx_op.get(
                operator_type='time_evolution',
                energy_shifted=False,  # Get physical U
                representation='dense_matrix'
            )

            # Verify unitarity
            identity = np.eye(dimension)
            exact_unitarity_error = np.linalg.norm(
                exact_unitary_matrix.conj().T @ exact_unitary_matrix - identity,
                'fro'
            )
            approx_unitarity_error = np.linalg.norm(
                approx_unitary_matrix.conj().T @ approx_unitary_matrix - identity,
                'fro'
            )
            logger.verbose(f"  U_exact unitarity check: ||U†U - I||_F = {exact_unitarity_error:.6e}")
            logger.verbose(f"  U_approx unitarity check: ||U†U - I||_F = {approx_unitarity_error:.6e}")

            if exact_unitarity_error > 1e-10:
                logger.warning(
                    f"WARNING: U_exact unitarity error {exact_unitarity_error:.6e} exceeds 1e-10."
                )
            if approx_unitarity_error > 1e-10:
                logger.warning(
                    f"WARNING: U_approx unitarity error {approx_unitarity_error:.6e} exceeds 1e-10."
                )

            diff_matrix = exact_unitary_matrix - approx_unitary_matrix

            for norm_type in norms_to_compute:
                if norm_type == 'frobenius':
                    error = np.linalg.norm(diff_matrix, 'fro')
                    logger.info(f"Frobenius norm ||U_exact - U_approx||_F: {error:.6e}")
                    results['matrix_frobenius_error'] = float(error)

                elif norm_type == 'spectral':
                    error = np.linalg.norm(diff_matrix, 2)
                    logger.info(f"Spectral norm ||U_exact - U_approx||_2: {error:.6e}")
                    results['matrix_spectral_error'] = float(error)

                else:
                    raise ValueError(f"Unknown matrix norm type: {norm_type}")

        else:
            # Matrix-free case: not implemented in Phase 2
            raise NotImplementedError(
                "Matrix-free matrix norm error analysis not yet implemented. "
                "Phase 2 supports only dense matrices (systems with n ≤ 15 qubits). "
                "This feature is planned for Phase 3 to support larger systems. "
                "Note: Matrix norm errors require O(N²) matrix-vector products in matrix-free mode, "
                "which may be impractical for very large systems anyway."
            )

    # =============================================================================================
    # 3. STATE-DEPENDENT ERRORS: ||U_exact|ψ⟩ - U_approx|ψ⟩||
    # =============================================================================================
    # Compares time evolution applied to specific quantum states.
    # Uses the same unshifted operators (U_exact, U_approx) as matrix norm errors.

    if config_analysis.error_state_inputs is not None:
        # Note: error_state_inputs is normalized to list during validation
        state_files = config_analysis.error_state_inputs

        logger.info(f"Computing state-dependent errors for {len(state_files)} state(s)")

        # Require matrices to be provided by caller
        if exact_matrix is None:
            raise ValueError(
                "State-dependent error analysis requires the exact Hamiltonian matrix, but it was not computed. "
            )

        if unitary_matrix is None:
            raise ValueError(
                "State-dependent error analysis requires the approximate/unitary matrix, but it was not computed. "
            )

        # Get unshifted time-evolution operators (reuse from matrix norm errors if already computed)
        if exact_op is None or approx_op is None:
            raise RuntimeError(
                "Operator representations should have been created but are None. "
                "This is an internal error in the operator conversion logic."
            )

        logger.verbose(f"  Getting U (unshifted, physical) for state evolution")
        exact_unitary_matrix = exact_op.get(
            operator_type='time_evolution',
            energy_shifted=False,  # Get physical U
            representation='dense_matrix'
        )

        approx_unitary_matrix = approx_op.get(
            operator_type='time_evolution',
            energy_shifted=False,  # Get physical U
            representation='dense_matrix'
        )

        state_errors = []

        for state_file in state_files:
            logger.verbose(f"Processing {state_file}")

            # Load state
            try:
                initial_state = load_state(state_file)
            except Exception as e:
                logger.info(f"ERROR: Failed to load state from {state_file}: {e}")
                raise

            # Apply exact time evolution operator: U_exact |ψ⟩
            # Note: For Phase 2 dense case, these are always numpy arrays (not matrix-free)
            exact_final = exact_unitary_matrix @ initial_state

            # Apply approximate time evolution operator: U_approx |ψ⟩
            approx_final = approx_unitary_matrix @ initial_state

            # Compute error: ||U_exact|ψ⟩ - U_approx|ψ⟩||
            diff = exact_final - approx_final
            absolute_error = np.linalg.norm(diff)
            relative_error = absolute_error / np.linalg.norm(exact_final)

            logger.info(f"  {state_file}: ||U_exact|ψ⟩ - U_approx|ψ⟩|| = {absolute_error:.6e} (relative: {relative_error:.6e})")

            state_errors.append({
                'input_file': state_file,
                'absolute_error': float(absolute_error),
                'relative_error': float(relative_error)
            })

        results['state_errors'] = state_errors

    # Save results to file
    output_file = 'error_analysis.npz'
    save_dict = {}

    if 'eigenenergy_errors' in results:
        save_dict['eigenenergy_absolute_errors'] = np.array(results['eigenenergy_errors']['absolute_errors'])
        save_dict['eigenenergy_relative_errors'] = np.array(results['eigenenergy_errors']['relative_errors'])
        save_dict['eigenenergy_num'] = results['eigenenergy_errors']['num_eigenstates']

    if 'matrix_frobenius_error' in results:
        save_dict['matrix_frobenius_error'] = results['matrix_frobenius_error']

    if 'matrix_spectral_error' in results:
        save_dict['matrix_spectral_error'] = results['matrix_spectral_error']

    if 'state_errors' in results:
        state_absolute = [s['absolute_error'] for s in results['state_errors']]
        state_relative = [s['relative_error'] for s in results['state_errors']]
        save_dict['state_absolute_errors'] = np.array(state_absolute)
        save_dict['state_relative_errors'] = np.array(state_relative)
        # Note: filenames are in results dict, not saved to npz

    if save_dict:
        np.savez(output_file, **save_dict)
        logger.info(f"Error analysis results saved to {output_file}")
        results['output_file'] = output_file

    return results

# -------------------------------------------------------------------------------------------------

def numerical_simulation(
        config_analysis: AnalysisConfiguration,
        algorithm,
        unitary_matrix) -> dict:
    """
    Perform numerical simulation by applying the unitary matrix to input state(s).

    Parameters:
        config_analysis: Analysis configuration with numerical_simulation_inputs
        algorithm: The algorithm bloq to analyze
        unitary_matrix: The unitary matrix to apply (pre-computed)

    Returns:
        Dictionary with simulation metadata: list of {input_file, output_file, output_norm}
    """

    # Log matrix properties
    logger.verbose(f"Matrix shape: {unitary_matrix.shape}")

    # Normalize inputs to list (in case this function is called directly without validation)
    raw_inputs = config_analysis.numerical_simulation_inputs
    if raw_inputs is None:
        raise ValueError("numerical_simulation_inputs is None")

    # Validate type before normalization
    if not isinstance(raw_inputs, (str, list)):
        raise ValueError(
            f"numerical_simulation_inputs must be a string or list of strings, "
            f"got {type(raw_inputs).__name__}"
        )

    input_files = _normalize_string_or_list_to_list(raw_inputs)

    logger.info(f"Running numerical simulation on {len(input_files)} input state(s)")

    results = []

    for input_file in input_files:
        logger.verbose(f"Processing {input_file}")

        # Load initial state
        try:
            initial_state = load_state(input_file)
        except Exception as e:
            logger.info(f"ERROR: Failed to load state from {input_file}: {e}")
            raise

        # Validate dimensions
        if initial_state.shape[0] != unitary_matrix.shape[1]:
            raise ValueError(
                f"Dimension mismatch: state vector has dimension {initial_state.shape[0]} "
                f"but matrix expects {unitary_matrix.shape[1]}"
            )

        # Perform matrix-vector multiplication
        logger.verbose("Performing matrix-vector multiplication")
        final_state = unitary_matrix @ initial_state

        # Compute norm
        final_norm = np.linalg.norm(final_state)
        logger.verbose(f"Final state norm: {final_norm:.6e}")

        # Generate output filename: input.npy -> input_final.npy
        input_path = Path(input_file)
        output_file = str(input_path.parent / f"{input_path.stem}_final{input_path.suffix}")

        # Save final state
        try:
            save_state(output_file, final_state)
        except Exception as e:
            logger.info(f"ERROR: Failed to save state to {output_file}: {e}")
            raise

        logger.info(f"Simulation complete: {input_file} -> {output_file}")

        results.append({
            'input_file': input_file,
            'output_file': output_file,
            'output_norm': float(final_norm)
        })

    return {'simulations': results}

# -------------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------------

def _normalize_string_or_list_to_list(value):
    """
    Normalize a configuration value that can be either a string or list into a list.

    This is a common pattern for config options that accept either a single item (string)
    or multiple items (list).

    Parameters:
        value: Either a string or a list of strings (or None)

    Returns:
        A list (or None if input was None)

    Examples:
        _normalize_string_or_list_to_list("item") -> ["item"]
        _normalize_string_or_list_to_list(["a", "b"]) -> ["a", "b"]
        _normalize_string_or_list_to_list(None) -> None
    """
    if value is None:
        return None
    elif isinstance(value, str):
        return [value]
    else:
        # Already a list (or list-like)
        return value

# -------------------------------------------------------------------------------------------------
# Functions to determine what expensive computations are required
# -------------------------------------------------------------------------------------------------

def requires_exact_eigendecomposition(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if exact eigendecomposition needs to be computed.

    Exact eigendecomposition is required for:
    - Eigendecomposition analysis with eigendecomposition_matrices = 'exact' or 'both'
    - Eigenvalue error analysis (always needs both eigendecompositions)

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if exact eigendecomposition computation is needed, False otherwise
    """
    # Eigendecomposition requested if eigendecomposition_matrices is not None
    eigendecomposition_requested = config_analysis.eigendecomposition_matrices is not None

    # Need exact eigendecomposition if:
    # 1. Eigendecomposition requested and matrices setting includes 'exact' or 'both'
    # 2. Eigenvalue error analysis is enabled (always needs both)
    return (
        (eigendecomposition_requested and
         config_analysis.eigendecomposition_matrices in ['exact', 'both']) or
        config_analysis.enable_eigenvalue_errors
    )


def requires_approximate_eigendecomposition(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if approximate eigendecomposition needs to be computed.

    Approximate eigendecomposition is required for:
    - Eigendecomposition analysis with eigendecomposition_matrices = 'approximate' or 'both'
    - Eigenvalue error analysis (always needs both eigendecompositions)

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if approximate eigendecomposition computation is needed, False otherwise
    """
    # Eigendecomposition requested if eigendecomposition_matrices is not None
    eigendecomposition_requested = config_analysis.eigendecomposition_matrices is not None

    # Need approximate eigendecomposition if:
    # 1. Eigendecomposition requested and matrices setting includes 'approximate' or 'both'
    # 2. Eigenvalue error analysis is enabled (always needs both)
    return (
        (eigendecomposition_requested and
         config_analysis.eigendecomposition_matrices in ['approximate', 'both']) or
        config_analysis.enable_eigenvalue_errors
    )


def requires_exact_matrix(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if the exact Hamiltonian matrix needs to be computed.

    The exact matrix is required for:
    - Exact matrix output to file
    - Exact eigendecomposition (which depends on the matrix)
    - Matrix norm error analysis
    - State-dependent error analysis

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if exact matrix computation is needed, False otherwise
    """
    return (
        config_analysis.exact_matrix_output_file is not None or
        requires_exact_eigendecomposition(config_analysis) or
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )


def requires_approximate_matrix(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if the approximate/unitary matrix needs to be computed.

    The approximate matrix is required for:
    - Matrix output to file
    - Numerical simulation
    - Approximate eigendecomposition (which depends on the matrix)
    - Matrix norm error analysis
    - State-dependent error analysis

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if approximate matrix computation is needed, False otherwise
    """
    return (
        config_analysis.matrix_output_file is not None or
        config_analysis.numerical_simulation_inputs is not None or
        requires_approximate_eigendecomposition(config_analysis) or
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )

# -------------------------------------------------------------------------------------------------

def validate_and_autocomplete_analysis_config(config_analysis: AnalysisConfiguration) -> None:
    """
    Validate configuration consistency and auto-enable dependent analyses where appropriate.

    This function is called early in driver.py, after loading configuration but before
    loading the Hamiltonian. This allows for fail-fast behavior if configuration is invalid.

    This function checks for:
    1. Missing dependencies (raises errors if configuration is needed)
    2. Opportunities to auto-enable analyses (logs when enabling)

    Parameters:
        config_analysis: Analysis configuration to validate and potentially modify

    Raises:
        ValueError: If configuration is inconsistent and cannot be auto-corrected

    Note:
        This function modifies the config_analysis object in-place when auto-enabling analyses.
    """

    # Check eigenvalue error analysis dependencies
    if config_analysis.enable_eigenvalue_errors:
        # Check if eigendecomposition is configured
        eigendecomposition_configured = config_analysis.eigendecomposition_matrices is not None

        if not eigendecomposition_configured:
            raise ValueError(
                "enable_eigenvalue_errors requires eigendecomposition. "
                "Set eigendecomposition_matrices to 'both' to enable eigenvalue error analysis."
            )

        # Must compute both eigendecompositions to compare
        if config_analysis.eigendecomposition_matrices != 'both':
            logger.info(
                "INFO: enable_eigenvalue_errors requires both exact and approximate eigendecompositions. "
                f"Auto-setting eigendecomposition_matrices from '{config_analysis.eigendecomposition_matrices}' to 'both'."
            )
            config_analysis.eigendecomposition_matrices = 'both'

    # Normalize string-or-list config values to always be lists
    # This allows downstream code to always assume list type
    config_analysis.error_matrix_norms = _normalize_string_or_list_to_list(
        config_analysis.error_matrix_norms
    )
    config_analysis.error_state_inputs = _normalize_string_or_list_to_list(
        config_analysis.error_state_inputs
    )
    config_analysis.numerical_simulation_inputs = _normalize_string_or_list_to_list(
        config_analysis.numerical_simulation_inputs
    )

    # Check if matrices will be computed and auto-enable output if not already set
    if requires_approximate_matrix(config_analysis):
        if config_analysis.matrix_output_file is None:
            default_filename = "unitary_matrix.npz"
            logger.info(
                f"INFO: Approximate/unitary matrix will be computed for requested analyses. "
                f"Auto-enabling matrix output to '{default_filename}' (essentially free)."
            )
            config_analysis.matrix_output_file = default_filename

    if requires_exact_matrix(config_analysis):
        if config_analysis.exact_matrix_output_file is None:
            default_filename = "exact_hamiltonian.npz"
            logger.info(
                f"INFO: Exact Hamiltonian matrix will be computed for requested analyses. "
                f"Auto-enabling exact matrix output to '{default_filename}' (essentially free)."
            )
            config_analysis.exact_matrix_output_file = default_filename

# -------------------------------------------------------------------------------------------------

def exact_numerical_simulation(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        exact_matrix) -> dict:
    """
    Perform exact numerical simulation by applying the exact Hamiltonian matrix to input state(s).

    This is similar to numerical_simulation() but uses the exact Hamiltonian matrix
    instead of the approximate unitary. Useful for comparing exact evolution with
    approximate algorithm results.

    Parameters:
        config_analysis: Analysis configuration with exact_simulation_inputs
        hamiltonian: Hamiltonian object (for metadata)
        exact_matrix: The exact Hamiltonian matrix to apply (pre-computed, can be matrix-free)

    Returns:
        Dictionary with simulation metadata: list of {input_file, output_file, output_norm}
    """

    # Log operator properties
    logger.verbose(f"Exact matrix shape: {exact_matrix.shape}")

    # Normalize input to list
    inputs = config_analysis.exact_simulation_inputs
    if isinstance(inputs, str):
        input_files = [inputs]
    elif isinstance(inputs, list):
        input_files = inputs
    else:
        raise ValueError(
            f"exact_simulation_inputs must be a string or list of strings, "
            f"got {type(inputs)}"
        )

    logger.info(f"Running exact numerical simulation on {len(input_files)} input state(s)")

    results = []

    for input_file in input_files:
        logger.verbose(f"Processing {input_file}")

        # Load initial state
        try:
            initial_state = load_state(input_file)
        except Exception as e:
            logger.info(f"ERROR: Failed to load state from {input_file}: {e}")
            raise

        # Validate dimensions
        if initial_state.shape[0] != exact_matrix.shape[1]:
            raise ValueError(
                f"Dimension mismatch: state vector has dimension {initial_state.shape[0]} "
                f"but exact matrix expects {exact_matrix.shape[1]}"
            )

        # Perform matrix-vector multiplication (or use matvec for matrix-free)
        logger.verbose("Performing exact matrix-vector multiplication")
        if hasattr(exact_matrix, 'matvec'):
            # Matrix-free operator
            final_state = exact_matrix.matvec(initial_state)
        else:
            # Dense matrix
            final_state = exact_matrix @ initial_state

        # Compute norm
        final_norm = np.linalg.norm(final_state)
        logger.verbose(f"Final state norm: {final_norm:.6e}")

        # Generate output filename: input.npy -> input_exact_final.npy
        input_path = Path(input_file)
        output_file = str(input_path.parent / f"{input_path.stem}_exact_final{input_path.suffix}")

        # Save final state
        try:
            save_state(output_file, final_state)
        except Exception as e:
            logger.info(f"ERROR: Failed to save state to {output_file}: {e}")
            raise

        logger.info(f"Exact simulation complete: {input_file} -> {output_file}")

        results.append({
            'input_file': input_file,
            'output_file': output_file,
            'output_norm': float(final_norm)
        })

    return {'exact_simulations': results}

# -------------------------------------------------------------------------------------------------

def save_requested_operator_outputs(
        config_analysis: AnalysisConfiguration,
        exact_matrix,
        unitary_matrix,
        timestep,
        energy_shift) -> dict:
    """
    Save all requested operator forms using OperatorRepresentation.

    Parameters
    ----------
    config_analysis : AnalysisConfiguration
        Configuration with output requests
    exact_matrix : ndarray
        Exact Hamiltonian matrix (shifted, H' = H + E*I)
    unitary_matrix : ndarray
        Approximate time-evolution operator (from shifted Hamiltonian)
    timestep : float
        Time evolution parameter
    energy_shift : float
        Energy shift value

    Returns
    -------
    dict
        Information about saved files
    """
    from qhat.analysis.operators import OperatorRepresentation
    from qhat.analysis.file_io import save_matrix, save_eigendecomposition
    import datetime

    results = {
        'matrix_outputs': [],
        'eigendecomposition_outputs': []
    }

    # Check if any requests exist
    if not config_analysis._matrix_output_requests and not config_analysis._eigendecomposition_output_requests:
        return results

    # Create OperatorRepresentation wrappers
    logger.verbose("Creating OperatorRepresentation instances for flexible output")
    logger.verbose(f"  Timestep: t = {timestep}")
    logger.verbose(f"  Energy shift: E = {energy_shift}")

    # Exact operator (from shifted Hamiltonian)
    exact_op = OperatorRepresentation(
        data=exact_matrix,
        operator_type='hamiltonian',
        energy_shifted=True,  # Input is H' = H + E*I
        representation='dense_matrix',
        timestep=timestep,
        energy_shift=energy_shift
    )
    logger.verbose("  Created exact operator representation (H', shifted)")

    # Approximate operator (from shifted Hamiltonian via Trotter/etc)
    approx_op = OperatorRepresentation(
        data=unitary_matrix,
        operator_type='time_evolution',
        energy_shifted=True,  # Input is U' from H'
        representation='dense_matrix',
        timestep=timestep,
        energy_shift=energy_shift
    )
    logger.verbose("  Created approximate operator representation (U', shifted)")

    operators = {
        'exact': exact_op,
        'approximate': approx_op
    }

    # Process matrix output requests
    for request in config_analysis._matrix_output_requests:
        op = operators[request['operator']]
        energy_shifted = (request['shift'] == 'shifted')

        logger.info(
            f"Saving {request['operator']} {request['form']} "
            f"({request['shift']}) matrix to {request['filename']}"
        )

        matrix = op.get(
            operator_type=request['form'],
            energy_shifted=energy_shifted,
            representation='dense_matrix'
        )

        # Save with metadata
        metadata = {
            'operator': request['operator'],
            'form': request['form'],
            'shift': request['shift'],
            'timestep': timestep,
            'energy_shift': energy_shift,
            'timestamp': datetime.datetime.now().isoformat(),
            'shape': matrix.shape
        }

        save_matrix(request['filename'], matrix, metadata=metadata)

        results['matrix_outputs'].append({
            'filename': request['filename'],
            'operator': request['operator'],
            'form': request['form'],
            'shift': request['shift'],
            'shape': matrix.shape
        })

    # Process eigendecomposition output requests
    for request in config_analysis._eigendecomposition_output_requests:
        op = operators[request['operator']]
        energy_shifted = (request['shift'] == 'shifted')

        logger.info(
            f"Saving {request['operator']} {request['form']} "
            f"({request['shift']}) eigendecomposition to {request['filename']}"
        )

        eigendata = op.get(
            operator_type=request['form'],
            energy_shifted=energy_shifted,
            representation='eigendecomposition'
        )

        # Sort by eigenvalues (ascending)
        sort_indices = np.argsort(eigendata['eigenvalues'].real)
        eigenvalues_sorted = eigendata['eigenvalues'][sort_indices]
        eigenvectors_sorted = eigendata['eigenvectors'][:, sort_indices]

        # Save with metadata
        metadata = {
            'operator': request['operator'],
            'form': request['form'],
            'shift': request['shift'],
            'timestep': timestep,
            'energy_shift': energy_shift,
            'timestamp': datetime.datetime.now().isoformat(),
            'dimension': len(eigenvalues_sorted)
        }

        save_eigendecomposition(
            request['filename'],
            eigenenergies=eigenvalues_sorted,
            eigenvectors=eigenvectors_sorted,
            matrix_type=f"{request['operator']}_{request['form']}_{request['shift']}",
            **metadata
        )

        results['eigendecomposition_outputs'].append({
            'filename': request['filename'],
            'operator': request['operator'],
            'form': request['form'],
            'shift': request['shift'],
            'num_eigenvalues': len(eigenvalues_sorted),
            'eigenvalue_range': [float(eigenvalues_sorted[0].real), float(eigenvalues_sorted[-1].real)]
        })

    return results

# -------------------------------------------------------------------------------------------------

def analyze_algorithm(
        config_analysis: AnalysisConfiguration,
        algorithm,
        hamiltonian=None,
        timestep=None,
        energy_shift=0.0) -> dict:
    """
    Analyze a quantum algorithm.

    Parameters:
        config_analysis: Analysis configuration
        algorithm: Algorithm bloq to analyze
        hamiltonian: Hamiltonian object (required for exact matrix/simulation)
        timestep: Time evolution parameter (required for approximate eigendecomposition)
                 If not provided, approximate eigendecomposition will fail.
        energy_shift: Energy shift applied to Hamiltonian (for correcting eigenvalue comparisons)

    Returns:
        Dictionary with analysis results
    """

    logger.info("Beginning algorithm analysis.")

    # Note: Configuration validation happens in driver.py before Hamiltonian is loaded

    # Check what analyses are requested
    eigendecomposition_requested = (
        requires_exact_eigendecomposition(config_analysis) or
        requires_approximate_eigendecomposition(config_analysis)
    )
    error_analysis_requested = (
        config_analysis.enable_eigenvalue_errors or
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )

    # Check if flexible output API is used
    flexible_outputs_requested = (
        config_analysis._matrix_output_requests or
        config_analysis._eigendecomposition_output_requests
    )

    # Validate at least one analysis requested
    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None and
        config_analysis.numerical_simulation_inputs is None and
        config_analysis.exact_matrix_output_file is None and
        not eigendecomposition_requested and
        not error_analysis_requested and
        config_analysis.exact_simulation_inputs is None and
        not flexible_outputs_requested):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')\n"
            "  - numerical_simulation_inputs (e.g., 'state.npy' or ['state1.npy', 'state2.npy'])\n"
            "  - exact_matrix_output_file (e.g., 'exact_hamiltonian.npz')\n"
            "  - eigendecomposition_matrices (e.g., 'exact', 'approximate', or 'both')\n"
            "  - enable_eigenvalue_errors (True to compute errors for all eigenvalues)\n"
            "  - error_matrix_norms (e.g., 'frobenius' or ['frobenius', 'spectral'])\n"
            "  - error_state_inputs (e.g., 'state.npy')\n"
            "  - exact_simulation_inputs (e.g., 'state.npy')\n"
            "  - analysis.save_matrix_to_file(...) (flexible matrix output API)\n"
            "  - analysis.save_eigendecomposition_to_file(...) (flexible eigendecomposition output API)"
        )

    results = {}

    # Determine what expensive computations are needed using shared functions
    needs_matrix = requires_approximate_matrix(config_analysis)
    needs_exact_matrix = requires_exact_matrix(config_analysis)

    # Compute matrices once if needed
    # Note: Opportunistic matrix output enabling happens during validation in driver.py
    unitary_matrix = None
    if needs_matrix:
        unitary_matrix = _compute_unitary_matrix(algorithm)

    exact_matrix = None
    if needs_exact_matrix:
        if hamiltonian is None:
            raise ValueError(
                "Exact matrix computation requires hamiltonian parameter. "
                "Pass hamiltonian to analyze_algorithm()."
            )
        exact_matrix = _compute_exact_matrix(hamiltonian, config_analysis)

    # Dispatch to requested analyses
    if config_analysis.resource_estimator is not None:
        logger.info(f"Performing resource estimation using {config_analysis.resource_estimator}.")
        results["resource_estimates"] = estimate_resources(config_analysis, algorithm)

    if config_analysis.matrix_output_file is not None:
        logger.info("Generating unitary matrix output.")
        results["matrix_output"] = output_unitary_matrix(config_analysis, algorithm, unitary_matrix)

    if config_analysis.exact_matrix_output_file is not None:
        logger.info("Generating exact Hamiltonian matrix output.")
        results["exact_matrix_output"] = exact_matrix_output(config_analysis, hamiltonian, exact_matrix)

    if config_analysis.numerical_simulation_inputs is not None:
        logger.info("Performing numerical simulation.")
        results["numerical_simulation"] = numerical_simulation(config_analysis, algorithm, unitary_matrix)

    # Eigendecomposition: single decision point for computing eigendecompositions
    exact_eigendecomp = None
    approx_eigendecomp = None
    if eigendecomposition_requested:
        logger.info("Performing eigendecomposition analysis.")
        eig_results = eigendecomposition_analysis(
            config_analysis,
            timestep=timestep,
            energy_shift=energy_shift,
            exact_matrix=exact_matrix,
            unitary_matrix=unitary_matrix
        )
        results["eigendecomposition"] = eig_results

        # Extract eigendecomposition data for use by error_analysis
        if 'exact_eigendecomposition' in eig_results:
            exact_eigendecomp = eig_results['exact_eigendecomposition']
        if 'approximate_eigendecomposition' in eig_results:
            approx_eigendecomp = eig_results['approximate_eigendecomposition']

    # Error analysis: receives eigendecomposition data, does not recompute
    if error_analysis_requested:
        logger.info("Performing error analysis.")
        results["error_analysis"] = error_analysis(
            config_analysis, hamiltonian, algorithm,
            exact_matrix=exact_matrix,
            unitary_matrix=unitary_matrix,
            exact_eigendecomp=exact_eigendecomp,
            approx_eigendecomp=approx_eigendecomp,
            timestep=timestep,
            energy_shift=energy_shift
        )

    if config_analysis.exact_simulation_inputs is not None:
        logger.info("Performing exact numerical simulation.")
        results["exact_simulation"] = exact_numerical_simulation(
            config_analysis, hamiltonian, exact_matrix
        )

    # Flexible operator outputs (new API)
    if config_analysis._matrix_output_requests or config_analysis._eigendecomposition_output_requests:
        logger.info("Saving requested operator outputs.")
        # Need both exact_matrix and unitary_matrix
        if exact_matrix is None:
            raise ValueError(
                "Flexible operator output requires exact_matrix. "
                "This should have been caught during validation."
            )
        if unitary_matrix is None:
            raise ValueError(
                "Flexible operator output requires unitary_matrix. "
                "This should have been caught during validation."
            )
        if timestep is None:
            raise ValueError(
                "Flexible operator output requires timestep parameter. "
                "Pass timestep to analyze_algorithm()."
            )

        results["flexible_operator_outputs"] = save_requested_operator_outputs(
            config_analysis,
            exact_matrix,
            unitary_matrix,
            timestep,
            energy_shift
        )

    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
