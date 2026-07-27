import logging

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.resource_estimation import estimate_resources
from qhat.analysis.error_analysis import error_analysis
from qhat.analysis.numerical_simulation import numerical_simulation
from qhat.analysis.matrix_eigendecomposition import (
    _compute_unitary_matrix,
    _compute_exact_matrix,
    output_unitary_matrix,
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
        req['source'] == 'exact' and req['representation'] == 'eigendecomposition'
        for req in config_analysis._operator_output_requests
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
        req['source'] == 'approximate' and req['representation'] == 'eigendecomposition'
        for req in config_analysis._operator_output_requests
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
        req['source'] == 'exact' and req['representation'] == 'matrix'
        for req in config_analysis._operator_output_requests
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
    # analysis.save_operator_to_file(source='exact', operator_type='hamiltonian', energy_shifted=False, representation='matrix', ...)

# -------------------------------------------------------------------------------------------------

def analyze_algorithm(
        config_analysis: AnalysisConfiguration,
        algorithm,
        unitary_encoding,
        approximate_time_evolution=None,
        exact_hamiltonian=None,
        timestep=None,
        energy_shift=0.0) -> dict:
    """
    Analyze a quantum algorithm.

    Parameters:
        config_analysis: Analysis configuration
        algorithm: Algorithm bloq to analyze
        unitary_encoding: Method used to encode the Hamiltonian into a unitary operator
        exact_hamiltonian: Hamiltonian object (required for exact matrix/simulation)
        timestep: Time evolution parameter (required for approximate eigendecomposition)
                 If not provided, approximate eigendecomposition will fail.
        energy_shift: Energy shift applied to Hamiltonian (for correcting eigenvalue comparisons)

    Returns:
        Dictionary with analysis results
    """

    logger.info("Beginning algorithm analysis.")

    # Note: Configuration validation happens in driver.py before Hamiltonian is loaded

    # Aggregate options to help summarize options _________________________________________________

    # TODO: Should some of these be methods of config_analysis?

    # Is the unitary encoding message Trotter?  Some analyses only apply to Trotter.
    is_trotter = unitary_encoding in ["ramped trotter"]

    # Are any error analyses turned on?
    error_analysis_requested = (
        config_analysis.enable_eigenvalue_errors or
        config_analysis.error_matrix_norms is not None or
        config_analysis.error_state_inputs is not None
    )

    # Consistency check: currently cannot do error analysis on non-Trotter methods
    if error_analysis_requested and not is_trotter:
        logger.warning("Error analysis was requested, but currently not available for non-Trotter "
                       "unitary encodings.  Turning off error analyses.")
        error_analysis_requested = False

    # Are exact and/or approximate OperatorRepresentation needed?
    exact_op_requests = [d for d in config_analysis._operator_output_requests
                         if d["source"] == "exact"]
    approx_op_requests = [d for d in config_analysis._operator_output_requests
                          if d["source"] == "approximate"]
    exact_op_needed = error_analysis_requested or len(exact_op_requests) > 0
    approximate_op_needed = error_analysis_requested or len(approx_op_requests) > 0

    # Consistency check: approximate operator cannot be constructed for non-Trotter methods
    # - We've already rules out error_analysis_requested, so that means that the user has asked
    #     for a Hamiltonian or time-evolution operator to be written to a file
    if approximate_op_needed and not is_trotter:
        logger.warning("Saving the approximate operator is currently unavailable for non-Trotter "
                       "unitary encodings.  Turning off these requests.")
        approx_op_requests = []

    # Is full algorithm matrix needed?
    algorithm_matrix_needed = (
        config_analysis.algorithm_matrix_output_file is not None or
        config_analysis.numerical_simulation_inputs is not None
    )

    # Validate at least one analysis requested
    if (not algorithm_matrix_needed and
        not error_analysis_requested and
        not exact_op_needed and
        not approximate_op_needed and 
        config_analysis.resource_estimator is None):
        raise ValueError("No analyses requested. Turn on at least one analysis.\n")

    logger.verbose("\n".join(["Analysis summary flags:",
        f"   error analysis requested:                  {error_analysis_requested}",
        f"   exact OperatorRepresentation needed:       {exact_op_needed}",
        f"   approximate OperatorRepresentation needed: {approximate_op_needed}",
        f"   full algorithm unitary matrix needed:      {algorithm_matrix_needed}",
        f"   number of exact operator output requests:  {len(exact_op_requests)}",
        f"   number of approx operator output requests: {len(approx_op_requests)}",
    ]))

    # Preliminary calculations ____________________________________________________________________

    # TODO: These blocks look like they could and perhaps should be functions

    # Construct the unitary matrix equivalent to the full algorithm
    algorithm_mat = None
    if algorithm_matrix_needed:
        logger.verbose("Constructing the unitary matrix equivalent to the full algorithm")
        algorithm_mat = _compute_unitary_matrix(algorithm)
        logger.info("Created full algorithm matrix")

    # Construct the "exact" (from the true Hamiltonian) OperatorRepresentation if necessary
    exact_op = None
    if exact_op_needed:
        if exact_hamiltonian is None:
            raise ValueError("Exact matrix/eigendecomposition computation requires exact_hamiltonian parameter.")
        logger.verbose("Constructing the exact energy-shifted Hamiltonian matrix")
        exact_matrix = _compute_exact_matrix(exact_hamiltonian, config_analysis.matrix_memory_threshold_gb)
        logger.verbose(
            "Constructing the exact OperatorRepresentation instance"
            f" (t = {timestep}, ΔE = {energy_shift})"
        )
        from qhat.analysis.operators import OperatorRepresentation
        exact_op = OperatorRepresentation(
            data=exact_matrix,
            operator_type='hamiltonian',
            energy_shifted=True,
            representation='dense_matrix',
            tevol_hbar=timestep,
            energy_shift=energy_shift
        )
        logger.info("Created exact operator representation")

    # Construct the "approximate" (from the unitary Hamiltonian encoding) OperatorRepresentation if
    # necessary
    approx_op = None
    if approximate_op_needed:
        if approximate_time_evolution is None:
            raise ValueError("Approximate matrix/eigendecomposition computation requires approximate_time_evolution parameter.")
        logger.verbose("Constructing the approximate energy-shifted time-evolution matrix")
        approx_matrix = _compute_unitary_matrix(approximate_time_evolution)
        logger.verbose(
            "Constructing the approximate OperatorRepresentation instance"
            f" (t = {timestep}, ΔE = {energy_shift})"
        )
        from qhat.analysis.operators import OperatorRepresentation
        approx_op = OperatorRepresentation(
            data=approx_matrix,
            operator_type='time_evolution',
            energy_shifted=True,
            representation='dense_matrix',
            tevol_hbar=timestep,
            energy_shift=energy_shift
        )
        logger.info("Created approximate operator representation")

    # Storage in which to aggregate the results ___________________________________________________

    results = {}

    # Perform analyses ____________________________________________________________________________

    # Analysis Category : resource estimation _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?
    if config_analysis.resource_estimator is not None:
        # TODO: Modify to allow resource estimation with multiple approaches
        logger.info(f"Performing resource estimation using {config_analysis.resource_estimator}.")
        results["resource_estimates"] = estimate_resources(config_analysis.resource_estimator, algorithm)
        print("keys = ", results["resource_estimates"].keys())

    # Analysis Category : error analysis  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    if error_analysis_requested:
        logger.info("Performing error analysis.")
        results["error_analysis"] = error_analysis(
            config_analysis, exact_hamiltonian, algorithm,
            exact_op=exact_op,
            approx_op=approx_op,
            timestep=timestep,
            energy_shift=energy_shift
        )

    # Analysis Category : matrices and eigendecompositions  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    # TODO: Why are we extracting a private data value?  Is it private just to "hide" it from users?
    if len(config_analysis._operator_output_requests) > 0:
        logger.info("Generating matrix and/or eigendecomposition outputs and saving to file(s).")
        me_results = dict()
        me_results.update(save_requested_operator_outputs(exact_op_requests, exact_op, "exact"))
        me_results.update(save_requested_operator_outputs(approx_op_requests, approx_op, "approx"))
        results["matrices_and_eigendecompositions"] = me_results
    if config_analysis.algorithm_matrix_output_file is not None:
        logger.info(f"Generating algorithm matrix output and saving to file {config_analysis.algorithm_matrix_output_file}.")
        results["matrix_output"] = output_unitary_matrix(config_analysis.algorithm_matrix_output_file, algorithm_mat)

    # Analysis Category : numerical simulation  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    if config_analysis.numerical_simulation_inputs is not None:
        logger.info("Performing numerical simulation of the full algorithm.")
        results["numerical_simulation"] = \
            numerical_simulation(config_analysis.numerical_simulation_inputs, algorithm, algorithm_mat)

    # Epiloque ____________________________________________________________________________________

    logger.info("Algorithm analysis complete.")

    return results
