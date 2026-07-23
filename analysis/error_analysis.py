import numpy as np
import logging

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import load_state
from qhat.analysis.matrix_operations import compute_unitarity_error

logger = logging.getLogger(__name__)

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
        energy_shift=0.0,
        exact_op=None,
        approx_op=None) -> dict:
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
                "Use the flexible API to save both eigendecompositions when enable_eigenvalue_errors is True."
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
                "Use the flexible API: analysis.save_matrix_to_file(operator='exact', ...) "
                "or enable eigendecomposition."
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
            # Matrix-free case: not implemented yet
            raise NotImplementedError(
                "Matrix/state error analysis not yet implemented for matrix-free operators."
            )

        # Create OperatorRepresentation instances if not provided
        if exact_op is None or approx_op is None:
            from qhat.analysis.operators import OperatorRepresentation

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
        else:
            logger.verbose("Using pre-created OperatorRepresentation instances (shared with other analyses)")

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
            exact_unitarity_error = compute_unitarity_error(exact_unitary_matrix)
            approx_unitarity_error = compute_unitarity_error(approx_unitary_matrix)
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
            # Matrix-free case: not implemented yet
            raise NotImplementedError(
                "Matrix-free matrix norm error analysis not yet implemented. "
                "Note: Matrix norm errors require O(N²) matrix-vector products in matrix-free "
                "mode, which may be impractical for very large systems anyway."
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
            # Note: For dense case, these are always numpy arrays (not matrix-free)
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
