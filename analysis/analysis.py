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
    # Note: save_matrix expects 'unitarity_error' parameter, but for Hamiltonians
    # we use hermiticity_error. Pass it as unitarity_error for storage.
    save_matrix(
        output_file, exact_matrix,
        unitarity_error=hermiticity_error,
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

def _eigendecompose(matrix, matrix_type, num_eigenvalues, which_eigs):
    """
    Perform eigendecomposition of a matrix (full or partial).

    Parameters:
        matrix: Matrix to decompose (dense array or matrix-free operator)
        matrix_type: String describing the matrix type (for logging)
        num_eigenvalues: Number of eigenvalues to compute (int) or "all" for full decomposition
        which_eigs: Which eigenvalues to compute ('smallest', 'largest', or 'both')

    Returns:
        tuple: (eigenvalues, eigenvectors, num_computed)

    Raises:
        ValueError: If parameters are invalid or operation not supported
    """
    import scipy.linalg
    import scipy.sparse.linalg
    from qhat.analysis.matrix_operations import PauliStringOperator

    dimension = matrix.shape[0]
    is_matrix_free = isinstance(matrix, PauliStringOperator) or hasattr(matrix, 'matvec')

    # Determine if full or partial eigendecomposition
    is_full = isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"

    if is_full:
        # Full eigendecomposition
        if is_matrix_free:
            raise ValueError(
                f"Full eigendecomposition not supported for matrix-free operators. "
                f"Use num_eigenvalues=k for partial decomposition."
            )
        logger.verbose(f"Computing full eigendecomposition for {matrix_type} matrix")
        eigenvalues, eigenvectors = scipy.linalg.eigh(matrix)
        num_computed = len(eigenvalues)

    else:
        # Partial eigendecomposition using sparse methods
        k = int(num_eigenvalues)
        if k <= 0:
            raise ValueError(f"num_eigenvalues must be positive, got {k}")
        if k > dimension:
            raise ValueError(
                f"num_eigenvalues ({k}) must be less than or equal to dimension ({dimension})."
            )

        # Map user-friendly values to scipy's 'which' parameter
        which_map = {
            'smallest': 'SA',  # Smallest Algebraic (most negative)
            'largest': 'LA',   # Largest Algebraic (most positive)
        }

        if which_eigs == 'both':
            # Compute both smallest and largest
            logger.verbose(f"Computing {k} smallest and {k} largest eigenvalues for {matrix_type} matrix")
            eigs_small, vecs_small = scipy.sparse.linalg.eigsh(
                matrix, k=k, which='SA', return_eigenvectors=True
            )
            eigs_large, vecs_large = scipy.sparse.linalg.eigsh(
                matrix, k=k, which='LA', return_eigenvectors=True
            )
            # Concatenate results
            eigenvalues = np.concatenate([eigs_small, eigs_large])
            eigenvectors = np.concatenate([vecs_small, vecs_large], axis=1)
            num_computed = len(eigenvalues)
        else:
            # Compute only one set
            which_scipy = which_map.get(which_eigs)
            if which_scipy is None:
                raise ValueError(
                    f"which_eigenvalues must be 'smallest', 'largest', or 'both', "
                    f"got '{which_eigs}'"
                )
            logger.verbose(f"Computing {k} {which_eigs} eigenvalues for {matrix_type} matrix")
            eigenvalues, eigenvectors = scipy.sparse.linalg.eigsh(
                matrix, k=k, which=which_scipy, return_eigenvectors=True
            )
            num_computed = len(eigenvalues)

    logger.info(f"Computed {num_computed} eigenvalues for {matrix_type} matrix")
    logger.verbose(f"  Eigenvalue range: [{eigenvalues.min():.6e}, {eigenvalues.max():.6e}]")

    return eigenvalues, eigenvectors, num_computed

# -------------------------------------------------------------------------------------------------

def _process_eigendecomposition(matrix, matrix_type, num_eigenvalues, which_eigs):
    """
    Process eigendecomposition for a single matrix: compute, save, and return data.

    Parameters:
        matrix: Matrix to decompose (dense array or matrix-free operator)
        matrix_type: String describing the matrix type ('exact' or 'approximate')
        num_eigenvalues: Number of eigenvalues to compute (int) or "all"
        which_eigs: Which eigenvalues to compute ('smallest', 'largest', or 'both')

    Returns:
        Dictionary with file path, num_eigenvalues, eigenvalue_range, which, eigenvalues, eigenvectors

    Raises:
        ValueError: If matrix is None
    """
    from qhat.analysis.file_io import save_eigendecomposition

    logger.info(f"Computing {matrix_type} matrix eigendecomposition")

    if matrix is None:
        raise ValueError(
            f"{matrix_type}_matrix is required for {matrix_type} eigendecomposition. "
            "Compute the matrix before calling eigendecomposition_analysis()."
        )

    eigs, vecs, num_computed = _eigendecompose(
        matrix, matrix_type, num_eigenvalues, which_eigs
    )

    # Save to file
    output_file = f'{matrix_type}_eigendecomposition.npz'
    save_eigendecomposition(
        output_file, eigs, vecs, matrix_type,
        num_eigenvalues, which_eigs
    )

    return {
        'file': output_file,
        'num_eigenvalues': num_computed,
        'eigenvalue_range': [float(eigs.min()), float(eigs.max())],
        'which': which_eigs,
        'eigenvalues': eigs,
        'eigenvectors': vecs
    }

# -------------------------------------------------------------------------------------------------

def eigendecomposition_analysis(
        config_analysis: AnalysisConfiguration,
        exact_matrix=None,
        unitary_matrix=None) -> dict:
    """
    Compute eigendecompositions of exact and/or approximate matrices.

    Parameters:
        config_analysis: Analysis configuration with eigendecomposition settings
        exact_matrix: Pre-computed exact matrix (required if eigendecomposition_matrices is 'exact' or 'both')
        unitary_matrix: Pre-computed unitary matrix (required if eigendecomposition_matrices is 'approximate' or 'both')

    Returns:
        Dictionary with eigendecomposition results and file paths

    Raises:
        ValueError: If required matrices are not provided
    """
    from qhat.analysis.file_io import save_eigendecomposition
    from qhat.analysis.matrix_operations import PauliStringOperator
    import scipy.linalg
    import scipy.sparse.linalg

    num_eigenvalues = config_analysis.num_eigenvalues
    which_matrices = config_analysis.eigendecomposition_matrices
    which_eigs = config_analysis.which_eigenvalues

    logger.info(f"Starting eigendecomposition analysis")
    logger.verbose(f"  num_eigenvalues: {num_eigenvalues}")
    logger.verbose(f"  eigendecomposition_matrices: {which_matrices}")
    logger.verbose(f"  which_eigenvalues: {which_eigs}")

    # Determine if full or partial eigendecomposition
    is_full = isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"

    # Determine which matrices we need
    need_exact = which_matrices in ['exact', 'both']
    need_approx = which_matrices in ['approximate', 'both']

    results = {}

    # Compute exact eigendecomposition if requested
    if need_exact:
        results['exact_eigendecomposition'] = _process_eigendecomposition(
            exact_matrix, 'exact', num_eigenvalues, which_eigs
        )

    # Compute approximate eigendecomposition if requested
    if need_approx:
        results['approximate_eigendecomposition'] = _process_eigendecomposition(
            unitary_matrix, 'approximate', num_eigenvalues, which_eigs
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
        approx_eigendecomp=None) -> dict:
    """
    Compute error metrics comparing exact and approximate representations.

    Three independent error types:
    1. Eigenvalue errors (if enable_eigenvalue_errors is True)
    2. Matrix norm errors (if error_matrix_norms is not None)
    3. State-dependent errors (if error_state_inputs is not None)

    Parameters:
        config_analysis: Analysis configuration with error analysis settings
        hamiltonian: Hamiltonian object
        algorithm: Algorithm bloq
        exact_matrix: Pre-computed exact matrix (optional)
        unitary_matrix: Pre-computed unitary matrix (optional)
        exact_eigendecomp: Pre-computed exact eigendecomposition (optional)
        approx_eigendecomp: Pre-computed approximate eigendecomposition (optional)

    Returns:
        Dictionary with error metrics
    """
    from qhat.analysis.matrix_operations import PauliStringOperator
    from qhat.analysis.file_io import load_eigendecomposition, load_state
    import scipy.linalg

    logger.info("Starting error analysis")

    results = {}

    # =============================================================================================
    # 1. EIGENVALUE ERROR
    # =============================================================================================

    if config_analysis.enable_eigenvalue_errors:
        logger.info("Computing eigenvalue errors for all eigenvalues from eigendecomposition")

        if exact_eigendecomp is None or approx_eigendecomp is None:
            raise ValueError(
                "Both eigendecompositions must be computed in order to compare eigenvalues. "
                "Ensure eigendecomposition_matrices is set to 'both' when enable_eigenvalue_errors is True."
            )

        # Get all eigenvalues from the eigendecompositions
        exact_eigs = exact_eigendecomp['eigenvalues']
        approx_eigs = approx_eigendecomp['eigenvalues']

        # Verify that the same number of eigenvalues were computed
        if len(exact_eigs) != len(approx_eigs):
            raise ValueError(
                f"Mismatch in number of eigenvalues: exact has {len(exact_eigs)}, "
                f"approximate has {len(approx_eigs)}. Both eigendecompositions must compute "
                "the same number of eigenvalues."
            )

        # Compare all eigenvalues
        absolute_errors = exact_eigs - approx_eigs
        relative_errors = absolute_errors / np.abs(exact_eigs)

        num_eigenvalues = len(exact_eigs)
        logger.info(f"Computed errors for {num_eigenvalues} eigenvalues")
        logger.info(f"Eigenvalue absolute error range: [{absolute_errors.min():.6e}, {absolute_errors.max():.6e}]")
        logger.info(f"Eigenvalue relative error range: [{relative_errors.min():.6e}, {relative_errors.max():.6e}]")

        results['eigenvalue_errors'] = {
            'num_eigenvalues': num_eigenvalues,
            'absolute_errors': absolute_errors.tolist(),
            'relative_errors': relative_errors.tolist(),
            'max_absolute_error': float(np.abs(absolute_errors).max()),
            'max_relative_error': float(np.abs(relative_errors).max())
        }

    # =============================================================================================
    # 2. MATRIX NORM ERRORS
    # =============================================================================================

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
            # Small systems: direct computation
            logger.verbose(f"Using dense matrices for norm computation (dimension={dimension})")
            diff_matrix = exact_matrix - unitary_matrix

            for norm_type in norms_to_compute:
                if norm_type == 'frobenius':
                    error = np.linalg.norm(diff_matrix, 'fro')
                    logger.info(f"Frobenius norm error: {error:.6e}")
                    results['matrix_frobenius_error'] = float(error)

                elif norm_type == 'spectral':
                    error = np.linalg.norm(diff_matrix, 2)
                    logger.info(f"Spectral norm error: {error:.6e}")
                    results['matrix_spectral_error'] = float(error)

                else:
                    raise ValueError(f"Unknown matrix norm type: {norm_type}")

        else:
            # Large systems: matrix-free computation
            logger.info(f"WARNING: Matrix-free norm computation for {num_qubits} qubits")
            logger.info(f"  This requires 2^{num_qubits} = {dimension:,} matrix-vector products")
            logger.info(f"  Estimated time: {'<1 minute' if dimension < 10000 else '1-10 minutes' if dimension < 100000 else '10+ minutes'}")

            for norm_type in norms_to_compute:
                if norm_type == 'frobenius':
                    # Compute ||A||_F^2 = sum_i ||A e_i||^2
                    logger.verbose("Computing Frobenius norm via matrix-vector products")
                    frobenius_squared = 0.0
                    for i in range(dimension):
                        if i % max(1, dimension // 10) == 0:
                            logger.verbose(f"  Progress: {i}/{dimension} ({100*i//dimension}%)")
                        # Create basis vector e_i
                        e_i = np.zeros(dimension, dtype=complex)
                        e_i[i] = 1.0
                        # Compute difference: (H - U) e_i
                        if hasattr(exact_matrix, 'matvec'):
                            exact_result = exact_matrix.matvec(e_i)
                        else:
                            exact_result = exact_matrix @ e_i
                        if hasattr(unitary_matrix, 'matvec'):
                            approx_result = unitary_matrix.matvec(e_i)
                        else:
                            approx_result = unitary_matrix @ e_i
                        diff_i = exact_result - approx_result
                        frobenius_squared += np.linalg.norm(diff_i) ** 2

                    error = np.sqrt(frobenius_squared)
                    logger.info(f"Frobenius norm error (matrix-free): {error:.6e}")
                    results['matrix_frobenius_error'] = float(error)

                elif norm_type == 'spectral':
                    # Compute ||A||_2 = largest singular value via power iteration
                    logger.verbose("Computing spectral norm via power iteration")

                    # Random starting vector
                    v = np.random.randn(dimension) + 1j * np.random.randn(dimension)
                    v = v / np.linalg.norm(v)

                    max_iterations = 100
                    tolerance = 1e-6

                    for iteration in range(max_iterations):
                        # Apply (H - U)† (H - U) to v
                        # First: (H - U) v
                        if hasattr(exact_matrix, 'matvec'):
                            exact_result = exact_matrix.matvec(v)
                        else:
                            exact_result = exact_matrix @ v
                        if hasattr(unitary_matrix, 'matvec'):
                            approx_result = unitary_matrix.matvec(v)
                        else:
                            approx_result = unitary_matrix @ v
                        diff_v = exact_result - approx_result

                        # Then: (H - U)† diff_v
                        if hasattr(exact_matrix, 'rmatvec'):
                            exact_adjoint = exact_matrix.rmatvec(diff_v)
                        else:
                            exact_adjoint = exact_matrix.conj().T @ diff_v
                        if hasattr(unitary_matrix, 'rmatvec'):
                            approx_adjoint = unitary_matrix.rmatvec(diff_v)
                        else:
                            approx_adjoint = unitary_matrix.conj().T @ diff_v
                        result = exact_adjoint - approx_adjoint

                        # Normalize
                        v_new = result / np.linalg.norm(result)

                        # Check convergence
                        if np.linalg.norm(v_new - v) < tolerance:
                            logger.verbose(f"  Converged after {iteration+1} iterations")
                            break

                        v = v_new

                        if (iteration + 1) % 10 == 0:
                            logger.verbose(f"  Iteration {iteration+1}/{max_iterations}")

                    # Rayleigh quotient to get eigenvalue (= squared singular value)
                    if hasattr(exact_matrix, 'matvec'):
                        exact_result = exact_matrix.matvec(v)
                    else:
                        exact_result = exact_matrix @ v
                    if hasattr(unitary_matrix, 'matvec'):
                        approx_result = unitary_matrix.matvec(v)
                    else:
                        approx_result = unitary_matrix @ v
                    diff_v = exact_result - approx_result

                    if hasattr(exact_matrix, 'rmatvec'):
                        exact_adjoint = exact_matrix.rmatvec(diff_v)
                    else:
                        exact_adjoint = exact_matrix.conj().T @ diff_v
                    if hasattr(unitary_matrix, 'rmatvec'):
                        approx_adjoint = unitary_matrix.rmatvec(diff_v)
                    else:
                        approx_adjoint = unitary_matrix.conj().T @ diff_v
                    result = exact_adjoint - approx_adjoint

                    eigenvalue = np.vdot(v, result)
                    error = np.sqrt(np.abs(eigenvalue))

                    logger.info(f"Spectral norm error (matrix-free, power iteration): {error:.6e}")
                    results['matrix_spectral_error'] = float(error)

                else:
                    raise ValueError(f"Unknown matrix norm type: {norm_type}")

    # =============================================================================================
    # 3. STATE-DEPENDENT ERRORS
    # =============================================================================================

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

        state_errors = []

        for state_file in state_files:
            logger.verbose(f"Processing {state_file}")

            # Load state
            try:
                initial_state = load_state(state_file)
            except Exception as e:
                logger.info(f"ERROR: Failed to load state from {state_file}: {e}")
                raise

            # Apply exact operator
            if hasattr(exact_matrix, 'matvec'):
                exact_final = exact_matrix.matvec(initial_state)
            else:
                exact_final = exact_matrix @ initial_state

            # Apply approximate operator
            if hasattr(unitary_matrix, 'matvec'):
                approx_final = unitary_matrix.matvec(initial_state)
            else:
                approx_final = unitary_matrix @ initial_state

            # Compute error
            diff = exact_final - approx_final
            absolute_error = np.linalg.norm(diff)
            relative_error = absolute_error / np.linalg.norm(exact_final)

            logger.info(f"  {state_file}: absolute error = {absolute_error:.6e}, relative error = {relative_error:.6e}")

            state_errors.append({
                'input_file': state_file,
                'absolute_error': float(absolute_error),
                'relative_error': float(relative_error)
            })

        results['state_errors'] = state_errors

    # Save results to file
    output_file = 'error_analysis.npz'
    save_dict = {}

    if 'eigenvalue_errors' in results:
        save_dict['eigenvalue_absolute_errors'] = np.array(results['eigenvalue_errors']['absolute_errors'])
        save_dict['eigenvalue_relative_errors'] = np.array(results['eigenvalue_errors']['relative_errors'])
        save_dict['eigenvalue_num'] = results['eigenvalue_errors']['num_eigenvalues']

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
    # Check if eigendecomposition is requested at all
    num_eigenvalues = config_analysis.num_eigenvalues
    eigendecomposition_requested = (
        isinstance(num_eigenvalues, int) and num_eigenvalues > 0
    ) or (
        isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"
    )

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
    # Check if eigendecomposition is requested at all
    num_eigenvalues = config_analysis.num_eigenvalues
    eigendecomposition_requested = (
        isinstance(num_eigenvalues, int) and num_eigenvalues > 0
    ) or (
        isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"
    )

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
        # Check if num_eigenvalues is configured
        num_eigenvalues = config_analysis.num_eigenvalues
        eigendecomposition_configured = (
            isinstance(num_eigenvalues, int) and num_eigenvalues > 0
        ) or (
            isinstance(num_eigenvalues, str) and num_eigenvalues.lower() == "all"
        )

        if not eigendecomposition_configured:
            raise ValueError(
                "enable_eigenvalue_errors requires eigendecomposition. "
                "Set num_eigenvalues to a positive integer or 'all'."
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

def analyze_algorithm(
        config_analysis: AnalysisConfiguration,
        algorithm,
        hamiltonian=None) -> dict:

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

    # Validate at least one analysis requested
    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None and
        config_analysis.numerical_simulation_inputs is None and
        config_analysis.exact_matrix_output_file is None and
        not eigendecomposition_requested and
        not error_analysis_requested and
        config_analysis.exact_simulation_inputs is None):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')\n"
            "  - numerical_simulation_inputs (e.g., 'state.npy' or ['state1.npy', 'state2.npy'])\n"
            "  - exact_matrix_output_file (e.g., 'exact_hamiltonian.npz')\n"
            "  - num_eigenvalues (e.g., 5 or 'all')\n"
            "  - enable_eigenvalue_errors (True to compute errors for all eigenvalues)\n"
            "  - error_matrix_norms (e.g., 'frobenius' or ['frobenius', 'spectral'])\n"
            "  - error_state_inputs (e.g., 'state.npy')\n"
            "  - exact_simulation_inputs (e.g., 'state.npy')"
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

    if eigendecomposition_requested:
        logger.info("Performing eigendecomposition analysis.")
        eig_results = eigendecomposition_analysis(
            config_analysis,
            exact_matrix=exact_matrix,
            unitary_matrix=unitary_matrix
        )
        # Keep eigenvalue/eigenvector data for error analysis, but remove it before storing
        exact_eigendecomp = None
        approx_eigendecomp = None
        if 'exact_eigendecomposition' in eig_results:
            exact_eigendecomp = eig_results['exact_eigendecomposition']
        if 'approximate_eigendecomposition' in eig_results:
            approx_eigendecomp = eig_results['approximate_eigendecomposition']

        # Remove numpy arrays from results dict (they're already saved to files)
        eig_results_for_storage = {}
        for key, value in eig_results.items():
            filtered_value = {k: v for k, v in value.items() if k not in ['eigenvalues', 'eigenvectors']}
            eig_results_for_storage[key] = filtered_value
        results["eigendecomposition"] = eig_results_for_storage

    if error_analysis_requested:
        logger.info("Performing error analysis.")
        # Use eigendecomposition data if available (from above)
        if not eigendecomposition_requested:
            exact_eigendecomp = None
            approx_eigendecomp = None
        results["error_analysis"] = error_analysis(
            config_analysis, hamiltonian, algorithm,
            exact_matrix=exact_matrix,
            unitary_matrix=unitary_matrix,
            exact_eigendecomp=exact_eigendecomp,
            approx_eigendecomp=approx_eigendecomp
        )

    if config_analysis.exact_simulation_inputs is not None:
        logger.info("Performing exact numerical simulation.")
        results["exact_simulation"] = exact_numerical_simulation(
            config_analysis, hamiltonian, exact_matrix
        )

    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
