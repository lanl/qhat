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

    Strategy: Convert H_exact → U_exact, then compare U eigenvalues directly.
    This avoids the ill-conditioned inverse problem (U → H via logarithm) which
    amplifies numerical errors by 1/t.

    The ratio λ_U_approx / λ_U_exact = exp(-i*ΔE*t/ℏ) directly encodes the
    energy error ΔE, which we recover via: ΔE = -angle(ratio) * ℏ/t

    Parameters:
        exact_op: OperatorRepresentation for exact Hamiltonian
        approx_op: OperatorRepresentation for approximate time evolution operator

    Returns:
        Dictionary with eigenenergy error metrics
    """
    logger.info("Computing eigenenergy errors for all eigenstates")

    if exact_op is None:
        raise ValueError("missing exact operator")
    if approx_op is None:
        raise ValueError("missing approximate operator")

    # Get exact H eigendecomposition first (native form for exact_op)
    exact_H_data = exact_op.get(operator_type="Hamiltonian",
                                 energy_shifted=False,
                                 representation="eigendecomposition")
    exact_H_eigenvalues = exact_H_data["eigenvalues"]
    exact_H_eigenvectors = exact_H_data["eigenvectors"]

    # Get exact U eigendecomposition (converted from H, but well-conditioned H→U)
    # H and U share the same eigenvectors (they commute)
    exact_U_data = exact_op.get(operator_type="time_evolution",
                                 energy_shifted=False,
                                 representation="eigendecomposition")
    exact_U_eigenvalues = exact_U_data["eigenvalues"]
    exact_U_eigenvectors = exact_U_data["eigenvectors"]

    # Get approximate U eigendecomposition (native form for approx_op - NO conversion!)
    approx_U_data = approx_op.get(operator_type="time_evolution",
                                   energy_shifted=False,
                                   representation="eigendecomposition")
    approx_U_eigenvalues = approx_U_data["eigenvalues"]
    approx_U_eigenvectors = approx_U_data["eigenvectors"]

    # Match approximate U eigenvectors to exact U eigenvectors by overlap
    # For each exact eigenvector, find the best matching approximate eigenvector
    overlaps = np.abs(exact_U_eigenvectors.conj().T @ approx_U_eigenvectors)  # [n_exact, n_approx]
    approx_indices = np.argmax(overlaps, axis=1)  # Best match for each exact state

    # Reorder approximate U eigenvalues to match exact ordering
    approx_U_eigenvalues = approx_U_eigenvalues[approx_indices]

    logger.debug(f"Eigenvector matching: min overlap = {overlaps.max(axis=1).min():.6f}, "
                 f"mean overlap = {overlaps.max(axis=1).mean():.6f}")

    # Sort everything by exact H eigenvalue for nice reporting
    energy_sort_indices = np.argsort(exact_H_eigenvalues)
    exact_H_eigenvalues = exact_H_eigenvalues[energy_sort_indices]
    exact_U_eigenvalues = exact_U_eigenvalues[energy_sort_indices]
    approx_U_eigenvalues = approx_U_eigenvalues[energy_sort_indices]

    # Verify same dimension
    if len(exact_U_eigenvalues) != len(approx_U_eigenvalues):
        raise ValueError(
            f"Dimension mismatch: exact operator has {len(exact_U_eigenvalues)} eigenstates, "
            f"approximate operator has {len(approx_U_eigenvalues)} eigenstates."
        )

    # Compute energy errors from U eigenvalue ratios
    # ratio = λ_U_approx / λ_U_exact = exp(-i*ΔE*t/ℏ)
    # Therefore: ΔE = -angle(ratio) * ℏ/t
    ratio = approx_U_eigenvalues / exact_U_eigenvalues
    phase_error = np.angle(ratio)  # ∈ (-π, π]

    # Convert phase error to energy error
    if exact_op.tevol_hbar is None:
        raise ValueError("tevol_hbar is required for eigenvalue error calculation")
    energy_errors = -phase_error / exact_op.tevol_hbar

    # Relative errors (use exact H eigenvalues as reference)
    relative_errors = energy_errors / np.abs(exact_H_eigenvalues)

    num_eigenstates = len(exact_H_eigenvalues)
    logger.info(f"Computed errors for {num_eigenstates} eigenstates (sorted by energy)")
    logger.verbose(f"  Ground state: exact={exact_H_eigenvalues[0]:.6e}, error={energy_errors[0]:.6e}, rel_err={relative_errors[0]:.6e}")
    logger.verbose(f"  Highest state: exact={exact_H_eigenvalues[-1]:.6e}, error={energy_errors[-1]:.6e}, rel_err={relative_errors[-1]:.6e}")
    logger.info(f"  Max absolute energy error: {np.abs(energy_errors).max():.6e}")
    logger.info(f"  Max relative energy error: {np.abs(relative_errors).max():.6e}")

    # Also log phase errors for comparison with matrix norms
    max_phase_error = np.abs(phase_error).max()
    logger.info(f"  Max phase error (for comparison with matrix norm): {max_phase_error:.6e}")

    return {
        'eigenenergy_errors': {
            'num_eigenstates': num_eigenstates,
            'absolute_errors': energy_errors.tolist(),
            'relative_errors': relative_errors.tolist(),
            'max_absolute_error': float(np.abs(energy_errors).max()),
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
            logger.error(f"Failed to load state from {state_file}: {e}")
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
