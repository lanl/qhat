"""QHAT-owned spatial-to-spin conversion with a configurable cutoff."""

import math
import numbers

import numpy as np


DEFAULT_COEFFICIENT_THRESHOLD = 1.0e-8


def validate_coefficient_threshold(value) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError("coefficient threshold must be a real number")
    threshold = float(value)
    if not math.isfinite(threshold) or threshold < 0.0:
        raise ValueError("coefficient threshold must be finite and non-negative")
    return threshold


def spinorb_from_spatial(
    one_body_integrals,
    two_body_integrals,
    coefficient_threshold=DEFAULT_COEFFICIENT_THRESHOLD,
):
    """Expand spatial integrals into spin-orbital tensors.

    This is the narrow piece of ``MolecularData.get_molecular_hamiltonian``
    QHAT needs to own in order to expose the otherwise fixed OpenFermion
    cutoff.  It preserves OpenFermion's strict ``abs(value) < threshold``
    boundary and supports ``0.0`` to disable magnitude truncation.
    """

    threshold = validate_coefficient_threshold(coefficient_threshold)
    n_qubits = 2 * one_body_integrals.shape[0]

    one_spin = np.zeros((n_qubits, n_qubits))
    two_spin = np.zeros((n_qubits,) * 4)

    one_spin[0::2, 0::2] = one_body_integrals
    one_spin[1::2, 1::2] = one_body_integrals

    two_spin[0::2, 1::2, 1::2, 0::2] = two_body_integrals
    two_spin[1::2, 0::2, 0::2, 1::2] = two_body_integrals
    two_spin[0::2, 0::2, 0::2, 0::2] = two_body_integrals
    two_spin[1::2, 1::2, 1::2, 1::2] = two_body_integrals

    one_spin[np.abs(one_spin) < threshold] = 0.0
    two_spin[np.abs(two_spin) < threshold] = 0.0

    return one_spin, two_spin
