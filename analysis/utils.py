"""
General-purpose utility functions for the analysis module.

This module contains helper functions that don't belong to any specific
domain like file I/O, configuration, or mathematical operations.
"""


def value(x, default):
    """
    Return x if it is not None, otherwise return default.

    This is a null-coalescing utility function that provides a concise way
    to handle optional values with defaults.

    Parameters:
        x: Value to check (can be None)
        default: Default value to return if x is None

    Returns:
        x if x is not None, otherwise default

    Examples:
        value(5, 10) -> 5
        value(None, 10) -> 10
        value(0, 10) -> 0  # Note: 0 is not None, so returns 0
    """
    if x is None:
        return default
    else:
        return x


def normalize_string_or_list_to_list(value):
    """
    Normalize a configuration value that can be either a string or list into a list.

    This is a common pattern for config options that accept either a single item (string)
    or multiple items (list).

    Parameters:
        value: Either a string or a list of strings (or None)

    Returns:
        A list (or None if input was None)

    Examples:
        normalize_string_or_list_to_list("item") -> ["item"]
        normalize_string_or_list_to_list(["a", "b"]) -> ["a", "b"]
        normalize_string_or_list_to_list(None) -> None
    """
    if value is None:
        return None
    elif isinstance(value, str):
        return [value]
    else:
        # Already a list (or list-like)
        return value
