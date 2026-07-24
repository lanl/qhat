import numpy as np
import logging
from pathlib import Path

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.file_io import load_state, save_state
from qhat.analysis.utils import normalize_string_or_list_to_list

logger = logging.getLogger(__name__)

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

    input_files = normalize_string_or_list_to_list(raw_inputs)

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

