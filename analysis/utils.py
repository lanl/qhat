"""
General-purpose utility functions for the analysis module.

This module contains helper functions that don't belong to any specific
domain like file I/O, configuration, or mathematical operations.
"""


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
