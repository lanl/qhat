"""
File I/O utilities for quantum algorithm analysis.

This module provides interfaces for writing unitary matrices in multiple formats.
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


def _save_matrix_text(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm):
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


def save_matrix(output_path, unitary_matrix, git_hash=None,
                unitarity_error=None, matrix_norm=None):
    """
    Save unitary matrix to file with automatic format detection.

    Parameters:
        output_path: Output file path (format inferred from extension)
        unitary_matrix: The matrix to save (numpy array)
        git_hash: Optional git hash for metadata
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

    handler(output_path, unitary_matrix, git_hash, unitarity_error, matrix_norm)
    logger.info(f"Matrix saved to {output_path}")
