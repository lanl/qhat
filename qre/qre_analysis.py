from qre_types import GeneralConfiguration, AnalysisConfiguration

from pyLIQTR.utils.resource_analysis import estimate_resources as estimate_pyliqtr
from qualtran.resource_counting import get_cost_value, QubitCount

import numpy as np
from datetime import datetime
from pathlib import Path

# -------------------------------------------------------------------------------------------------

def resource_estimation_cirq(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:
    raise NotImplementedError

# -------------------------------------------------------------------------------------------------

def resource_estimation_pyliqtr(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    config_general.log_verbose("Estimating resources with pyLIQTR.")

    # TODO: rotation error
    #       -- argument rotation_gate_precision sets the precision for a single rotation gate
    #       -- argument circuit_precision sets the precision for the whole circuit (i.e., it sets
    #          rotation_gate_precision to circuit_precision / number of rotations)
    # TODO: profile?
    #       -- argument profile = True: keep rotations as a separate count
    #       -- argument profile = False: estimate rotations as Clifford+T
    resources = estimate_pyliqtr(qpe_circuit)

    resource_dict = {
        "Clifford_count" : resources["Clifford"],
        "T_count"        : resources["T"],
        }
    if "LogicalQubits" in resources:
        resource_dict["qubit_count"] = resources["LogicalQubits"]
    else:
        get_cost_value(qpe_circuit, QubitCount())
    return resource_dict

# -------------------------------------------------------------------------------------------------

def estimate_resources(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    if config_analysis.resource_estimator.lower() == "pyliqtr":
        return resource_estimation_pyliqtr(config_general, config_analysis, qpe_circuit)
    elif config_analysis.resource_estimator.lower() == "cirq":
        return resource_estimation_cirq(config_general, config_analysis, qpe_circuit)
    else:
        raise ValueError(
                f"Invalid resource estimator method \"{config_analysis.resource_estimator}\".")

# -------------------------------------------------------------------------------------------------
# Matrix file format handlers
# -------------------------------------------------------------------------------------------------

def _save_matrix_numpy(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm):
    """Save matrix in NumPy .npz format with metadata."""
    np.savez(
        output_path,
        matrix=unitary_matrix,
        shape=unitary_matrix.shape,
        git_hash=git_hash,
        timestamp=datetime.now().isoformat(),
        unitarity_error=unitarity_error,
        matrix_norm=matrix_norm
    )

def _save_matrix_hdf5(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm):
    """Save matrix in HDF5 format with metadata as attributes."""
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required for HDF5 output")

    with h5py.File(output_path, 'w') as f:
        dset = f.create_dataset('matrix', data=unitary_matrix, compression='gzip')
        dset.attrs['shape'] = unitary_matrix.shape
        dset.attrs['git_hash'] = git_hash
        dset.attrs['timestamp'] = datetime.now().isoformat()
        if unitarity_error is not None:
            dset.attrs['unitarity_error'] = float(unitarity_error)
        if matrix_norm is not None:
            dset.attrs['matrix_norm'] = float(matrix_norm)

def _save_matrix_text(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm, config_general):
    """Save matrix in human-readable sparse text format (coordinate format).

    Only non-zero entries are written to save space. This is especially beneficial
    for large sparse matrices common in quantum chemistry.
    """
    # Threshold for considering an entry "zero"
    # For float64 (machine precision ~2e-16), entries below 1e-15 are likely accumulated
    # rounding errors from tensor contraction rather than physically meaningful values.
    # For unitary matrices with ||U||_F = sqrt(N) and entries bounded by 1, this threshold
    # is conservative: well above numerical noise but below any realistic quantum amplitude.
    threshold = 1e-15

    # Count non-zero entries for logging
    nnz = np.sum(np.abs(unitary_matrix) > threshold)
    sparsity = 1.0 - (nnz / unitary_matrix.size)

    config_general.log_verbose(
        f"Matrix sparsity: {sparsity*100:.2f}% ({nnz} non-zero of {unitary_matrix.size} total)"
    )

    with open(output_path, 'w') as f:
        f.write(f"# Unitary Matrix (Sparse Coordinate Format)\n")
        f.write(f"# Shape: {unitary_matrix.shape}\n")
        f.write(f"# Non-zero entries: {nnz}\n")
        f.write(f"# Git hash: {git_hash}\n")
        f.write(f"# Timestamp: {datetime.now().isoformat()}\n")
        if unitarity_error is not None:
            f.write(f"# Unitarity error: {unitarity_error:.6e}\n")
        if matrix_norm is not None:
            f.write(f"# Matrix norm: {matrix_norm:.6e}\n")
        f.write(f"# Format: row,col,real,imag (only non-zero entries, |value| > {threshold})\n")
        f.write(f"#\n")

        # Write only non-zero entries
        for i in range(unitary_matrix.shape[0]):
            for j in range(unitary_matrix.shape[1]):
                val = unitary_matrix[i, j]
                if abs(val) > threshold:
                    f.write(f"{i},{j},{val.real:.16e},{val.imag:.16e}\n")

def _get_format_from_extension(filename):
    """
    Infer matrix format from file extension.

    Returns: (format_name, extension) tuple
    Raises: ValueError if extension is not recognized
    """
    extension_map = {
        '.npz': 'numpy',
        '.h5': 'hdf5',
        '.hdf5': 'hdf5',
        '.txt': 'text',
        '.dat': 'text',
        '.csv': 'text',
    }

    # Get extension from filename
    ext = Path(filename).suffix.lower()

    if ext not in extension_map:
        raise ValueError(
            f"Cannot determine format from extension '{ext}'. "
            f"Supported extensions: {', '.join(extension_map.keys())}"
        )

    return extension_map[ext]

def _save_matrix_file(output_path, matrix_format, unitary_matrix, git_hash,
                      unitarity_error, matrix_norm, config_general):
    """
    Save matrix to file in the specified format.

    Dispatches to format-specific save functions.
    """
    format_handlers = {
        'numpy': _save_matrix_numpy,
        'hdf5': _save_matrix_hdf5,
        'text': _save_matrix_text,
    }

    handler = format_handlers.get(matrix_format)
    if handler is None:
        raise ValueError(
            f"Invalid matrix output format: {matrix_format}. "
            f"Valid options are: {', '.join(repr(k) for k in format_handlers.keys())}"
        )

    # Text format needs config_general for warnings, others don't
    if matrix_format == 'text':
        handler(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm, config_general)
    else:
        handler(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm)

    config_general.log(f"Matrix saved to {output_path}")

# -------------------------------------------------------------------------------------------------

def output_unitary_matrix(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:
    """
    Generate and save the unitary matrix representation of the circuit.

    Parameters:
        config_general: General configuration
        config_analysis: Analysis configuration with matrix_output_format and matrix_output_file
        qpe_circuit: The circuit bloq to analyze

    Returns:
        Dictionary with matrix metadata: shape, file, format, unitarity_error, norm
    """

    # Check if the circuit has tensor_contract method
    if not hasattr(qpe_circuit, 'tensor_contract'):
        raise AttributeError(
            f"Circuit of type {type(qpe_circuit).__name__} does not have a "
            "'tensor_contract()' method. Cannot generate unitary matrix."
        )

    # Extract the unitary matrix
    config_general.log_verbose("Computing unitary matrix via tensor contraction...")
    try:
        unitary_matrix = qpe_circuit.tensor_contract()
    except Exception as e:
        config_general.log(
            f"ERROR: Failed to compute unitary matrix: {e}\n"
            "This may indicate a bug in the circuit's tensor_contract() implementation."
        )
        raise

    # Log basic properties
    config_general.log_verbose(f"Matrix shape: {unitary_matrix.shape}")
    config_general.log_verbose(f"Matrix dtype: {unitary_matrix.dtype}")

    # Compute unitarity check: ||U†U - I||_F
    try:
        matrix_norm = np.linalg.norm(unitary_matrix, ord='fro')
        U_dag_U = np.conj(unitary_matrix.T) @ unitary_matrix
        identity = np.eye(unitary_matrix.shape[0])
        unitarity_error = np.linalg.norm(U_dag_U - identity, ord='fro')
        config_general.log_verbose(f"Matrix Frobenius norm: {matrix_norm:.6e}")
        config_general.log_verbose(f"Unitarity error ||U†U - I||_F: {unitarity_error:.6e}")
    except Exception as e:
        config_general.log(f"WARNING: Could not compute unitarity check: {e}")
        matrix_norm = None
        unitarity_error = None

    # Determine output filename and format from extension
    output_file = config_analysis.matrix_output_file
    matrix_format = _get_format_from_extension(output_file)

    output_path = Path(output_file)

    # Save matrix to file
    config_general.log(f"Saving matrix to {output_path} in {matrix_format} format...")
    _save_matrix_file(
        output_path, matrix_format, unitary_matrix,
        config_general.git_hash, unitarity_error, matrix_norm,
        config_general
    )

    # Return metadata
    return {
        'matrix_shape': unitary_matrix.shape,
        'matrix_file': str(output_path),
        'matrix_format': matrix_format,
        'unitarity_error': float(unitarity_error) if unitarity_error is not None else None,
        'matrix_norm': float(matrix_norm) if matrix_norm is not None else None
    }

# -------------------------------------------------------------------------------------------------

def analyze_circuit(
        config_general: GeneralConfiguration,
        config_analysis: AnalysisConfiguration,
        qpe_circuit) -> dict:

    config_general.log("Beginning circuit analysis.")

    # Validate at least one analysis requested
    if (config_analysis.resource_estimator is None and
        config_analysis.matrix_output_file is None):
        raise ValueError(
            "No analyses requested. Set at least one of:\n"
            "  - resource_estimator (e.g., 'pyliqtr', 'cirq')\n"
            "  - matrix_output_file (e.g., 'matrix.npz', 'matrix.h5', 'matrix.txt')"
        )

    results = {}

    # Dispatch to requested analyses
    if config_analysis.resource_estimator is not None:
        config_general.log(f"Performing resource estimation using {config_analysis.resource_estimator}.")
        results["resource_estimates"] = estimate_resources(
            config_general, config_analysis, qpe_circuit)

    if config_analysis.matrix_output_file is not None:
        config_general.log("Generating unitary matrix output.")
        results["matrix_output"] = output_unitary_matrix(
            config_general, config_analysis, qpe_circuit)

    # TODO: Add error estimation
    # TODO: Add an option for detailed error analysis (explicitly compute the eigenvalues of the
    #       original Hamiltonian and the final unitary, compute ground state energy from both,
    #       compare the results; will only work for small systems)
    # TODO: Add gate parallelism / gate depth analysis
    # TODO: Would it be useful to analyze in terms of a different basis (e.g., Toffoli gates)?

    config_general.log("Circuit analysis complete.")

    return results
