import logging

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.resource_estimation import estimate_resources
from qhat.analysis.error_analysis import error_analysis
from qhat.analysis.numerical_simulation import numerical_simulation, exact_numerical_simulation
from qhat.analysis.matrix_eigendecomposition import (
    _compute_unitary_matrix,
    _compute_exact_matrix,
    output_unitary_matrix,
    eigendecomposition_analysis,
    save_requested_operator_outputs
)
from qhat.analysis.utils import normalize_string_or_list_to_list

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------
# Functions to determine what expensive computations are required
# -------------------------------------------------------------------------------------------------

def requires_exact_eigendecomposition(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if exact eigendecomposition needs to be computed.

    Exact eigendecomposition is required for:
    - Flexible API requests exact eigendecompositions
    - Eigenvalue error analysis (always needs both eigendecompositions)

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if exact eigendecomposition computation is needed, False otherwise
    """
    # Check if flexible API requests exact operator eigendecompositions
    has_exact_eigendecomp_request = any(
        req['operator'] == 'exact'
        for req in config_analysis._eigendecomposition_output_requests
    )

    # Need exact eigendecomposition if:
    # 1. Flexible API requests it
    # 2. Eigenvalue error analysis is enabled (always needs both)
    return (
        has_exact_eigendecomp_request or
        config_analysis.enable_eigenvalue_errors
    )


def requires_approximate_eigendecomposition(config_analysis: AnalysisConfiguration) -> bool:
    """
    Determine if approximate eigendecomposition needs to be computed.

    Approximate eigendecomposition is required for:
    - Flexible API requests approximate eigendecompositions
    - Eigenvalue error analysis (always needs both eigendecompositions)

    Parameters:
        config_analysis: Analysis configuration

    Returns:
        True if approximate eigendecomposition computation is needed, False otherwise
    """
    # Check if flexible API requests approximate operator eigendecompositions
    has_approx_eigendecomp_request = any(
        req['operator'] == 'approximate'
        for req in config_analysis._eigendecomposition_output_requests
    )

    # Need approximate eigendecomposition if:
    # 1. Flexible API requests it
    # 2. Eigenvalue error analysis is enabled (always needs both)
    return (
        has_approx_eigendecomp_request or
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
    # Check if flexible API requests exact operator matrices
    has_exact_flexible_output = any(
        req['operator'] == 'exact'
        for req in config_analysis._matrix_output_requests
    )

    return (
        has_exact_flexible_output or
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
        config_analysis.algorithm_matrix_output_file is not None or
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
    # Note: If enable_eigenvalue_errors is True, the requires_*_eigendecomposition() functions
    # will return True, ensuring eigendecompositions are computed even if not explicitly requested
    # for output. No validation error needed - the system auto-enables what's required.

    # Normalize string-or-list config values to always be lists
    # This allows downstream code to always assume list type
    config_analysis.error_matrix_norms = normalize_string_or_list_to_list(
        config_analysis.error_matrix_norms
    )
    config_analysis.error_state_inputs = normalize_string_or_list_to_list(
        config_analysis.error_state_inputs
    )
    config_analysis.numerical_simulation_inputs = normalize_string_or_list_to_list(
        config_analysis.numerical_simulation_inputs
    )

    # Check if matrices will be computed and auto-enable output if not already set
    if requires_approximate_matrix(config_analysis):
        if config_analysis.algorithm_matrix_output_file is None:
            default_filename = "unitary_matrix.npz"
            logger.info(
                f"INFO: Approximate/unitary matrix will be computed for requested analyses. "
                f"Auto-enabling matrix output to '{default_filename}' (essentially free)."
            )
            config_analysis.algorithm_matrix_output_file = default_filename

    # Note: No auto-enabling for exact matrix - users should use flexible API
    # analysis.save_matrix_to_file(operator='exact', form='hamiltonian', ...)

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
        config_analysis.algorithm_matrix_output_file is None and
        config_analysis.numerical_simulation_inputs is None and
        not eigendecomposition_requested and
        not error_analysis_requested and
        config_analysis.exact_simulation_inputs is None and
        not flexible_outputs_requested):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - algorithm_matrix_output_file (e.g., 'matrix.npz') - full algorithm circuit matrix\n"
            "  - numerical_simulation_inputs (e.g., 'state.npy' or ['state1.npy', 'state2.npy'])\n"
            "  - enable_eigenvalue_errors (True to compute errors for all eigenvalues)\n"
            "  - error_matrix_norms (e.g., 'frobenius' or ['frobenius', 'spectral'])\n"
            "  - error_state_inputs (e.g., 'state.npy')\n"
            "  - exact_simulation_inputs (e.g., 'state.npy')\n"
            "  - analysis.save_matrix_to_file(...) - flexible matrix output API\n"
            "  - analysis.save_eigendecomposition_to_file(...) - flexible eigendecomposition output API"
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

    if config_analysis.algorithm_matrix_output_file is not None:
        logger.info("Generating algorithm matrix output.")
        results["matrix_output"] = output_unitary_matrix(config_analysis, algorithm, unitary_matrix)

    # exact_matrix_output removed - use flexible API instead

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
            unitary_matrix=unitary_matrix,
            requires_exact_eigendecomposition_func=requires_exact_eigendecomposition,
            requires_approximate_eigendecomposition_func=requires_approximate_eigendecomposition
        )
        results["eigendecomposition"] = eig_results

        # Extract eigendecomposition data for use by error_analysis
        if 'exact_eigendecomposition' in eig_results:
            exact_eigendecomp = eig_results['exact_eigendecomposition']
        if 'approximate_eigendecomposition' in eig_results:
            approx_eigendecomp = eig_results['approximate_eigendecomposition']

    # Create OperatorRepresentation instances once for reuse across analyses
    # This avoids redundant eigendecompositions when both error analysis and flexible outputs are requested
    exact_op = None
    approx_op = None
    needs_operators = (
        error_analysis_requested or
        len(config_analysis._matrix_output_requests) > 0 or
        len(config_analysis._eigendecomposition_output_requests) > 0
    )
    if needs_operators:
        from qhat.analysis.operators import OperatorRepresentation

        logger.info("Creating shared OperatorRepresentation instances")
        logger.verbose(f"  Timestep: t = {timestep}")
        logger.verbose(f"  Energy shift: E = {energy_shift}")

        # Exact operator (from shifted Hamiltonian H' = H + E*I)
        exact_op = OperatorRepresentation(
            data=exact_matrix,
            operator_type='hamiltonian',
            energy_shifted=True,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )
        logger.verbose("  Created exact operator representation (H', shifted)")

        # Approximate operator (from shifted time-evolution U' = exp(-i*H'*t))
        approx_op = OperatorRepresentation(
            data=unitary_matrix,
            operator_type='time_evolution',
            energy_shifted=True,
            representation='dense_matrix',
            timestep=timestep,
            energy_shift=energy_shift
        )
        logger.verbose("  Created approximate operator representation (U', shifted)")

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
            energy_shift=energy_shift,
            exact_op=exact_op,
            approx_op=approx_op
        )

    if config_analysis.exact_simulation_inputs is not None:
        logger.info("Performing exact numerical simulation.")
        results["exact_simulation"] = exact_numerical_simulation(
            config_analysis, hamiltonian, exact_matrix
        )

    # Flexible operator outputs
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
            energy_shift,
            exact_op=exact_op,
            approx_op=approx_op
        )

    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    logger.info("Algorithm analysis complete.")

    return results
