"""Small QHAT-owned hotfix for OpenFermion's fixed coefficient cutoff."""

from collections import defaultdict
import math
import numbers

import numpy as np
from openfermion import (
    FermionOperator,
    InteractionOperator,
    QubitOperator,
    bravyi_kitaev,
    jordan_wigner,
)


DEFAULT_COEFFICIENT_THRESHOLD = 1.0e-8


def validate_coefficient_threshold(value) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError("coefficient threshold must be a real number")
    threshold = float(value)
    if not math.isfinite(threshold) or threshold < 0.0:
        raise ValueError("coefficient threshold must be finite and non-negative")
    return threshold


def _validate_spatial_integrals(one_body_integrals, two_body_integrals):
    one_body = np.asarray(one_body_integrals)
    two_body = np.asarray(two_body_integrals)

    if one_body.ndim != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError("one-body spatial integrals must be a square rank-2 array")
    expected_two_body_shape = (one_body.shape[0],) * 4
    if two_body.shape != expected_two_body_shape:
        raise ValueError(
            "two-body spatial integrals must have shape "
            f"{expected_two_body_shape}, got {two_body.shape}"
        )
    return one_body, two_body


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
    one_spatial, two_spatial = _validate_spatial_integrals(
        one_body_integrals, two_body_integrals
    )
    n_qubits = 2 * one_spatial.shape[0]

    one_spin = np.zeros((n_qubits, n_qubits), dtype=one_spatial.dtype)
    two_spin = np.zeros((n_qubits,) * 4, dtype=two_spatial.dtype)

    one_spin[0::2, 0::2] = one_spatial
    one_spin[1::2, 1::2] = one_spatial

    two_spin[0::2, 1::2, 1::2, 0::2] = two_spatial
    two_spin[1::2, 0::2, 0::2, 1::2] = two_spatial
    two_spin[0::2, 0::2, 0::2, 0::2] = two_spatial
    two_spin[1::2, 1::2, 1::2, 1::2] = two_spatial

    if threshold > 0.0:
        one_spin[np.abs(one_spin) < threshold] = 0.0
        two_spin[np.abs(two_spin) < threshold] = 0.0

    return one_spin, two_spin


def molecular_hamiltonian_from_spatial(
    constant,
    one_body_integrals,
    two_body_integrals,
    coefficient_threshold=DEFAULT_COEFFICIENT_THRESHOLD,
):
    """Construct an ``InteractionOperator`` and record cutoff provenance."""

    threshold = validate_coefficient_threshold(coefficient_threshold)
    one_spatial, two_spatial = _validate_spatial_integrals(
        one_body_integrals, two_body_integrals
    )
    one_before = 2 * int(np.count_nonzero(one_spatial))
    two_before = 4 * int(np.count_nonzero(two_spatial))
    one_spin, two_spin = spinorb_from_spatial(one_spatial, two_spatial, threshold)

    operator = InteractionOperator(constant, one_spin, 0.5 * two_spin)
    operator.coefficient_threshold = threshold
    operator.threshold_stats = {
        "one_body_before": one_before,
        "one_body_after": int(np.count_nonzero(one_spin)),
        "one_body_removed": one_before - int(np.count_nonzero(one_spin)),
        "two_body_before": two_before,
        "two_body_after": int(np.count_nonzero(two_spin)),
        "two_body_removed": two_before - int(np.count_nonzero(two_spin)),
    }
    return operator


def get_molecular_hamiltonian(
    molecule,
    occupied_indices=None,
    active_indices=None,
    coefficient_threshold=DEFAULT_COEFFICIENT_THRESHOLD,
):
    """QHAT-local counterpart of ``MolecularData.get_molecular_hamiltonian``."""

    if occupied_indices is None and active_indices is None:
        one_body_integrals, two_body_integrals = molecule.get_integrals()
        constant = molecule.nuclear_repulsion
    else:
        core_adjustment, one_body_integrals, two_body_integrals = (
            molecule.get_active_space_integrals(occupied_indices, active_indices)
        )
        constant = molecule.nuclear_repulsion + core_adjustment

    return molecular_hamiltonian_from_spatial(
        constant,
        one_body_integrals,
        two_body_integrals,
        coefficient_threshold,
    )


def map_interaction_operator(
    operator,
    mapping_name,
    coefficient_threshold=DEFAULT_COEFFICIENT_THRESHOLD,
):
    """Map without changing OpenFermion's global ``EQ_TOLERANCE``.

    Each unit-coefficient monomial safely survives OpenFermion's internal
    cutoff.  QHAT restores the real coefficient, sums matching Pauli terms,
    and applies the requested cutoff once at the end.
    """

    threshold = validate_coefficient_threshold(coefficient_threshold)
    if mapping_name not in ("jordan-wigner", "bravyi-kitaev"):
        raise ValueError(f"invalid fermion-to-qubit mapping {mapping_name!r}")

    contributions = defaultdict(list)
    for fermion_term in operator:
        coefficient = operator[fermion_term]
        unit_term = FermionOperator(fermion_term, 1.0)
        if mapping_name == "jordan-wigner":
            mapped_unit = jordan_wigner(unit_term)
        else:
            mapped_unit = bravyi_kitaev(unit_term, n_qubits=operator.n_qubits)
        for pauli_term, unit_coefficient in mapped_unit.terms.items():
            contributions[pauli_term].append(coefficient * unit_coefficient)

    terms = {}
    for pauli_term, values in contributions.items():
        coefficient = complex(
            math.fsum(complex(value).real for value in values),
            math.fsum(complex(value).imag for value in values),
        )
        # Molecular Hamiltonians are Hermitian.  Permit only machine-scale
        # imaginary residue from floating-point tensor symmetries.
        scale = math.fsum(abs(value) for value in values)
        imag_tolerance = 128.0 * np.finfo(float).eps * max(1.0, scale)
        if abs(coefficient.imag) > imag_tolerance:
            raise ValueError(f"non-Hermitian Pauli coefficient {coefficient!r}")
        coefficient = coefficient.real
        if coefficient != 0.0 and abs(coefficient) >= threshold:
            terms[pauli_term] = coefficient

    result = QubitOperator()
    # Do not use += here; it would reintroduce OpenFermion's fixed cutoff.
    result.terms = terms
    return result
