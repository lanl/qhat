"""
File I/O utilities for quantum algorithm analysis.

This module provides clean interfaces for reading and writing:
- Unitary matrices (NumPy .npz, HDF5 .h5/.hdf5, text .txt/.dat/.csv)
- Quantum state vectors (NumPy .npy)

Format is automatically detected from file extension.
"""

import logging
import numpy as np
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)


# =================================================================================================
# Matrix I/O Functions
# =================================================================================================

def _save_matrix_numpy(output_path, unitary_matrix, unitarity_error, matrix_norm):
    """Save matrix in NumPy .npz format with metadata."""
    np.savez(
        output_path,
        matrix=unitary_matrix,
        shape=unitary_matrix.shape,
        timestamp=datetime.now().isoformat(),
        unitarity_error=unitarity_error,
        matrix_norm=matrix_norm
    )


def _save_matrix_hdf5(output_path, unitary_matrix, unitarity_error, matrix_norm):
    """Save matrix in HDF5 format with metadata as attributes."""
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required for HDF5 output")

    with h5py.File(output_path, 'w') as f:
        dset = f.create_dataset('matrix', data=unitary_matrix, compression='gzip')
        dset.attrs['shape'] = unitary_matrix.shape
        dset.attrs['timestamp'] = datetime.now().isoformat()
        if unitarity_error is not None:
            dset.attrs['unitarity_error'] = float(unitarity_error)
        if matrix_norm is not None:
            dset.attrs['matrix_norm'] = float(matrix_norm)


def _save_matrix_text(output_path, unitary_matrix, unitarity_error, matrix_norm):
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

    # Count non-zero entries
    nnz = np.sum(np.abs(unitary_matrix) > threshold)
    sparsity = 1.0 - (nnz / unitary_matrix.size)

    with open(output_path, 'w') as f:
        f.write(f"# Unitary Matrix (Sparse Coordinate Format)\n")
        f.write(f"# Shape: {unitary_matrix.shape}\n")
        f.write(f"# Non-zero entries: {nnz}\n")
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


def _get_matrix_format_from_extension(filename):
    """
    Infer matrix format from file extension.

    Returns: format_name
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


def save_matrix(output_path, unitary_matrix, unitarity_error=None, matrix_norm=None):
    """
    Save unitary matrix to file with automatic format detection.

    Parameters:
        output_path: Output file path (format inferred from extension)
        unitary_matrix: The matrix to save (numpy array)
        unitarity_error: Optional unitarity error for metadata
        matrix_norm: Optional matrix norm for metadata

    Supported formats:
        - .npz: NumPy compressed format
        - .h5, .hdf5: HDF5 format with compression
        - .txt, .dat, .csv: Human-readable sparse text format

    Raises:
        ValueError: If format not recognized
        ImportError: If h5py not available for HDF5 format
    """
    # Determine format from extension
    matrix_format = _get_matrix_format_from_extension(output_path)

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

    handler(output_path, unitary_matrix, unitarity_error, matrix_norm)
    logger.info(f"Matrix saved to {output_path}")


# =================================================================================================
# State Vector I/O Functions
# =================================================================================================

def _load_state_numpy(path):
    """Load state vector from NumPy .npy format."""
    return np.load(path, allow_pickle=False)


def _save_state_numpy(path, vec):
    """Save state vector to NumPy .npy format."""
    np.save(path, vec, allow_pickle=False)


def _get_state_format_from_extension(filename):
    """
    Infer state vector format from file extension.

    Returns: format_name
    Raises: ValueError if extension is not recognized
    """
    extension_map = {
        '.npy': 'numpy',
    }

    # Get extension from filename
    ext = Path(filename).suffix.lower()

    if ext not in extension_map:
        raise ValueError(
            f"Cannot determine format from extension '{ext}'. "
            f"Supported extensions: {', '.join(extension_map.keys())}"
        )

    return extension_map[ext]


# =================================================================================================
# Eigendecomposition I/O Functions
# =================================================================================================

def save_eigendecomposition(output_path, eigenvalues, eigenvectors, matrix_type,
                            num_eigenvalues, which_eigenvalues, git_hash=None):
    """
    Save eigendecomposition results to NumPy .npz format.

    Parameters:
        output_path: Path to output file (must end in .npz)
        eigenvalues: 1D array of eigenvalues
        eigenvectors: 2D array where column i is the eigenvector for eigenvalue i
        matrix_type: 'exact' or 'approximate'
        num_eigenvalues: Number of eigenvalues computed (int or "all")
        which_eigenvalues: 'smallest', 'largest', or 'both'
        git_hash: Optional git hash for provenance
    """
    if not output_path.endswith('.npz'):
        raise ValueError(f"Output path must end with .npz, got: {output_path}")

    metadata = {
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'matrix_type': matrix_type,
        'num_eigenvalues': num_eigenvalues if isinstance(num_eigenvalues, int) else str(num_eigenvalues),
        'which_eigenvalues': which_eigenvalues,
        'timestamp': datetime.now().isoformat(),
    }

    if git_hash is not None:
        metadata['git_hash'] = git_hash

    np.savez(output_path, **metadata)
    logger.info(f"Eigendecomposition saved to {output_path}")


def load_eigendecomposition(path):
    """
    Load eigendecomposition results from NumPy .npz format.

    Parameters:
        path: Path to .npz file

    Returns:
        Dictionary with keys: eigenvalues, eigenvectors, matrix_type,
        num_eigenvalues, which_eigenvalues, timestamp, and optionally git_hash
    """
    if not path.endswith('.npz'):
        raise ValueError(f"Path must end with .npz, got: {path}")

    data = np.load(path, allow_pickle=False)
    result = {
        'eigenvalues': data['eigenvalues'],
        'eigenvectors': data['eigenvectors'],
        'matrix_type': str(data['matrix_type']),
        'num_eigenvalues': data['num_eigenvalues'],
        'which_eigenvalues': str(data['which_eigenvalues']),
        'timestamp': str(data['timestamp']),
    }

    if 'git_hash' in data:
        result['git_hash'] = str(data['git_hash'])

    logger.info(f"Eigendecomposition loaded from {path}")
    return result


def load_state(path):
    """
    Load quantum state vector from file with automatic format detection.

    Parameters:
        path: Path to state vector file (format inferred from extension)

    Returns:
        numpy array containing the state vector (1D complex array)

    Supported formats:
        - .npy: NumPy binary format

    Validation:
        - Converts to complex dtype if needed
        - Checks dimension is power of 2
        - Ensures 1D array

    Raises:
        ValueError: If format not recognized or vector invalid
        FileNotFoundError: If file doesn't exist
    """
    format_handlers = {
        'numpy': _load_state_numpy,
    }

    # Determine format from extension
    state_format = _get_state_format_from_extension(path)

    handler = format_handlers.get(state_format)
    if handler is None:
        raise ValueError(
            f"Invalid state format: {state_format}. "
            f"Valid options are: {', '.join(repr(k) for k in format_handlers.keys())}"
        )

    # Load the state
    logger.verbose(f"Loading state vector from {path}")
    vec = handler(path)

    # Validate
    if not isinstance(vec, np.ndarray):
        raise ValueError(f"Loaded state is not a numpy array (got {type(vec)})")
    if vec.ndim != 1:
        raise ValueError(f"State vector must be 1-dimensional (got shape {vec.shape})")
    if not np.iscomplexobj(vec):
        vec = vec.astype(complex)
        logger.verbose("Converted state vector to complex dtype")

    # Check if dimension is a power of 2
    n = vec.shape[0]
    if n == 0 or (n & (n - 1)) != 0:
        raise ValueError(f"State vector dimension {n} is not a power of 2")

    logger.verbose(f"Loaded state vector with dimension {n} (n_qubits={int(np.log2(n))})")

    return vec


def save_state(path, vec):
    """
    Save quantum state vector to file with automatic format detection.

    Parameters:
        path: Output file path (format inferred from extension)
        vec: State vector (numpy array)

    Supported formats:
        - .npy: NumPy binary format

    Raises:
        ValueError: If format not recognized
    """
    format_handlers = {
        'numpy': _save_state_numpy,
    }

    # Determine format from extension
    state_format = _get_state_format_from_extension(path)

    handler = format_handlers.get(state_format)
    if handler is None:
        raise ValueError(
            f"Invalid state format: {state_format}. "
            f"Valid options are: {', '.join(repr(k) for k in format_handlers.keys())}"
        )

    logger.verbose(f"Saving state vector to {path}")
    handler(path, vec)
    logger.verbose(f"State vector saved successfully")
