"""Resource estimation for ground-state phase estimation with partially
randomized product formulas.  The public entry point is :func:`estimate_resources`.
"""

from __future__ import annotations

from .api import (
    ResourceEstimate,
    SplitEstimate,
    estimate_resources,
)

__all__ = [
    "estimate_resources",
    "ResourceEstimate",
    "SplitEstimate",
]
