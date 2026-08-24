"""Configuration utilities shared across QHAT tools."""

import hashlib


def parse_key_value_params(param_strings):
    """
    Parse KEY=VALUE parameter strings from command line.

    Parameters
    ----------
    param_strings : list of str
        List of "KEY=VALUE" strings

    Returns
    -------
    dict
        Dictionary mapping keys to values, with type conversion applied

    Raises
    ------
    ValueError
        If a parameter string is not in KEY=VALUE format

    Notes
    -----
    - Attempts to evaluate values as Python literals (int, float, list, etc.)
    - Falls back to string if evaluation fails
    - Example: "distance=1.5" → {"distance": 1.5}
    - Example: "name=hello" → {"name": "hello"}
    """
    config_params = {}
    for param in param_strings:
        if '=' not in param:
            raise ValueError(f"Parameter must be in KEY=VALUE format, got: {param}")
        key, value = param.split('=', 1)
        # Try to evaluate as Python literal (numbers, lists, etc.)
        try:
            config_params[key] = eval(value)
        except:
            # If evaluation fails, treat as string
            config_params[key] = value
    return config_params


def meV_to_Hartree(meV):
    """
    Convert millielectronvolts (meV) to Hartree atomic units.

    Parameters
    ----------
    meV : float
        Energy in millielectronvolts

    Returns
    -------
    float
        Energy in Hartree atomic units

    Notes
    -----
    Conversion factor: 1 meV = 3.67493221757e-5 Hartree
    """
    return 3.67493221757e-5 * meV


def string_to_seed(s):
    """
    Convert a string to a deterministic integer seed for random number generators.

    Parameters
    ----------
    s : str
        Input string

    Returns
    -------
    int
        Deterministic seed derived from the string (fits in 64-bit integer)

    Notes
    -----
    - Uses SHA-256 hash for deterministic, uniform distribution
    - Same string always produces same seed
    - Different strings produce very different seeds
    """
    hash_bytes = hashlib.sha256(s.encode('utf-8')).digest()
    # Take first 8 bytes and convert to integer (fits in 64-bit)
    seed = int.from_bytes(hash_bytes[:8], byteorder='big')
    return seed


def get_standard_exec_namespace():
    """
    Get a dictionary of standard utility functions for config file execution.

    Returns
    -------
    dict
        Dictionary of utility function names to functions

    Notes
    -----
    This provides the standard utilities available in config files:
    - meV_to_Hartree: Energy unit conversion
    - string_to_seed: Deterministic seed generation

    Tools should add their own config objects (general, hamiltonian, etc.)
    and params dictionary to this namespace before executing config files.
    """
    return {
        'meV_to_Hartree': meV_to_Hartree,
        'string_to_seed': string_to_seed,
    }
