import numpy as np
import logging

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import load_state
from qhat.analysis.matrix_operations import compute_unitarity_error
from qhat.analysis.utils import normalize_string_or_list_to_list

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def _compute_eigenvalue_errors(exact_op, approx_op) -> dict:
    """
    Compute eigenenergy errors between exact and approximate Hamiltonians.

    Parameters:
        exact_op: OperatorRepresentation for exact Hamiltonian
        approx_op: OperatorRepresentation for approximate Hamiltonian

    Returns:
        Dictionary with eigenenergy error metrics
    """
    logger.info("Computing eigenenergy errors for all eigenstates")

    if exact_op is None:
        raise ValueError("missing exact operator")
    if approx_op is None:
        raise ValueError("missing approximate operator")

    # Get eigenenergies from both decompositions (both already sorted by energy)
    def get_eigenvalues(op):
        return np.sort(op.get(operator_type="Hamiltonian",
                      energy_shifted=False,
                      representation="eigendecomposition")["eigenvalues"])
    exact_eigenenergies = get_eigenvalues(exact_op)
    approx_eigenenergies = get_eigenvalues(approx_op)

    # Verify same dimension (should have all eigenvalues since full decomposition)
    if len(exact_eigenenergies) != len(approx_eigenenergies):
        raise ValueError(
            f"Dimension mismatch: exact operator has {len(exact_eigenenergies)} eigenstates, "
            f"approximate operator has {len(approx_eigenenergies)} eigenstates."
        )

    # Element-wise comparison (both already sorted by energy)
    absolute_errors = exact_eigenenergies - approx_eigenenergies
    relative_errors = absolute_errors / np.abs(exact_eigenenergies)

    num_eigenstates = len(exact_eigenenergies)
    logger.info(f"Computed errors for {num_eigenstates} eigenstates (sorted by energy)")
    logger.verbose(f"  Ground state: exact={exact_eigenenergies[0]:.6e}, approx={approx_eigenenergies[0]:.6e}, error={absolute_errors[0]:.6e}")
    logger.verbose(f"  Highest state: exact={exact_eigenenergies[-1]:.6e}, approx={approx_eigenenergies[-1]:.6e}, error={absolute_errors[-1]:.6e}")
    logger.verbose(f"  Max absolute error: {np.abs(absolute_errors).max():.6e}")
    logger.verbose(f"  Max relative error: {np.abs(relative_errors).max():.6e}")

    return {
        'eigenenergy_errors': {
            'num_eigenstates': num_eigenstates,
            'absolute_errors': absolute_errors.tolist(),
            'relative_errors': relative_errors.tolist(),
            'max_absolute_error': float(np.abs(absolute_errors).max()),
            'max_relative_error': float(np.abs(relative_errors).max())
        }
    }

# -------------------------------------------------------------------------------------------------

def _compute_matrix_norm_errors(exact_unitary_matrix, approx_unitary_matrix, norms_to_compute) -> dict:
    """
    Compute matrix norm errors between exact and approximate time evolution operators.

    Parameters:
        exact_unitary_matrix: Dense matrix for exact time evolution operator
        approx_unitary_matrix: Dense matrix for approximate time evolution operator
        norms_to_compute: List of norm types to compute ('frobenius', 'spectral')

    Returns:
        Dictionary with matrix norm error metrics
    """
    logger.info(f"Computing matrix norm errors: {norms_to_compute}")

    # Verify unitarity
    exact_unitarity_error = compute_unitarity_error(exact_unitary_matrix)
    approx_unitarity_error = compute_unitarity_error(approx_unitary_matrix)
    logger.verbose(f"  U_exact unitarity check: ||U†U - I||_F = {exact_unitarity_error:.6e}")
    logger.verbose(f"  U_approx unitarity check: ||U†U - I||_F = {approx_unitarity_error:.6e}")

    if exact_unitarity_error > 1e-10:
        logger.warning(
            f"WARNING: U_exact unitarity error {exact_unitarity_error:.6e} > 1e-10."
        )
    if approx_unitarity_error > 1e-10:
        logger.warning(
            f"WARNING: U_approx unitarity error {approx_unitarity_error:.6e} > 1e-10."
        )

    diff_matrix = exact_unitary_matrix - approx_unitary_matrix

    results = {}
    for norm_type in normalize_string_or_list_to_list(norms_to_compute):
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

    return results

# -------------------------------------------------------------------------------------------------

def _compute_state_dependent_errors(exact_unitary_matrix, approx_unitary_matrix, state_files) -> dict:
    """
    Compute state-dependent errors for time evolution of specific quantum states.

    Parameters:
        exact_unitary_matrix: Dense matrix for exact time evolution operator
        approx_unitary_matrix: Dense matrix for approximate time evolution operator
        state_files: List of file paths containing quantum states to evolve

    Returns:
        Dictionary with state-dependent error metrics
    """
    logger.info(f"Computing state-dependent errors for {len(state_files)} state(s)")

    state_errors = []

    for state_file in normalize_string_or_list_to_list(state_files):
        logger.verbose(f"Processing {state_file}")

        # Load state
        try:
            initial_state = load_state(state_file)
        except Exception as e:
            logger.error(f"ERROR: Failed to load state from {state_file}: {e}")
            raise

        # Apply exact time evolution operator: U_exact |ψ⟩
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

    return {'state_errors': state_errors}

# -------------------------------------------------------------------------------------------------

def _save_error_results(results: dict, config_general=None) -> str:
    """
    Save error analysis results to npz file.

    Parameters:
        results: Dictionary containing error analysis results
        config_general: General configuration for output directory handling

    Returns:
        Path to output file
    """
    output_file = 'error_analysis.npz'
    if config_general:
        output_file = config_general.get_output_path(output_file)
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

    return output_file

# -------------------------------------------------------------------------------------------------

def error_analysis(
        config_analysis: AnalysisConfiguration,
        hamiltonian,
        algorithm,
        exact_op,
        approx_op,
        timestep=None,
        energy_shift=0.0,
        config_general=None) -> dict:
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
        exact_eigendecomp: Pre-computed exact eigendecomposition (optional)
        approx_eigendecomp: Pre-computed approximate eigendecomposition (optional)
        timestep: Time evolution parameter t/ħ (required for matrix/state errors).
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

    # Determine if we need to extract unitary matrices (shared by matrix norm and state errors)
    need_unitaries = (
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )

    # Extract unitary matrices once if needed by multiple error types
    exact_unitary_matrix = None
    approx_unitary_matrix = None
    if need_unitaries:
        if exact_op is None or approx_op is None:
            raise RuntimeError(
                "Operator representations should have been created but are None. "
                "This is an internal error in the operator conversion logic."
            )

        logger.verbose("Converting shifted H to unshifted U, constructing matrix")
        exact_unitary_matrix = exact_op.get(
            operator_type='time_evolution',
            energy_shifted=False,
            representation='dense_matrix'
        )
        logger.verbose("Converting shifted U to unshifted U, constructing matrix")
        approx_unitary_matrix = approx_op.get(
            operator_type='time_evolution',
            energy_shifted=False,
            representation='dense_matrix'
        )

    # Dispatch to independent error computation functions
    if config_analysis.enable_eigenvalue_errors:
        results.update(_compute_eigenvalue_errors(exact_op, approx_op))

    if config_analysis.error_matrix_norms is not None:
        results.update(_compute_matrix_norm_errors(
            exact_unitary_matrix,
            approx_unitary_matrix,
            config_analysis.error_matrix_norms
        ))

    if config_analysis.error_state_inputs is not None:
        results.update(_compute_state_dependent_errors(
            exact_unitary_matrix,
            approx_unitary_matrix,
            config_analysis.error_state_inputs
        ))

    # Save results to file
    if results:
        output_file = _save_error_results(results, config_general)
        results['output_file'] = output_file

    return results
