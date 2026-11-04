"""Majorana expansion utilities for OpenFermion ``FermionOperator`` objects."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from openfermion import FermionOperator


def single_majorana_as_fermion_operator(k: int) -> FermionOperator:
    """Return Majorana ``gamma_k`` as an OpenFermion ``FermionOperator``."""
    p, bit = divmod(k, 2)

    if bit == 0:
        return FermionOperator(((p, 0),), 1.0) + FermionOperator(((p, 1),), 1.0)

    return FermionOperator(((p, 1),), 1.0j) + FermionOperator(((p, 0),), -1.0j)


def majorana_monomial_as_fermion_operator(label: Iterable[int]) -> FermionOperator:
    """Convert a canonical Majorana label into a ``FermionOperator``."""
    op = FermionOperator((), 1.0)

    for k in label:
        op *= single_majorana_as_fermion_operator(k)

    op.compress()
    return op


def multiply_majorana_label_by_generator(label: Iterable[int], k: int):
    """Multiply a canonical Majorana monomial by ``gamma_k`` on the right."""
    label = tuple(label)
    num_greater = sum(1 for j in label if j > k)
    sign = -1 if num_greater % 2 else 1

    if k in label:
        new_label = tuple(j for j in label if j != k)
    else:
        new_label = tuple(sorted(label + (k,)))

    return new_label, sign


def fermion_term_to_majorana_terms(term, coeff: complex) -> dict[tuple[int, ...], complex]:
    """Expand one ladder-operator monomial into Majorana monomials."""
    acc = {(): complex(coeff)}

    for p, action in term:
        if action == 1:
            choices = ((2 * p, 0.5), (2 * p + 1, -0.5j))
        elif action == 0:
            choices = ((2 * p, 0.5), (2 * p + 1, 0.5j))
        else:
            raise ValueError(f"Unexpected ladder action {action}; expected 0 or 1.")

        new_acc = defaultdict(complex)

        for label, c0 in acc.items():
            for gamma_index, c1 in choices:
                new_label, sign = multiply_majorana_label_by_generator(label, gamma_index)
                new_acc[new_label] += c0 * c1 * sign

        acc = dict(new_acc)

    return {label: c for label, c in acc.items() if abs(c) > 0.0}


def fermion_operator_to_majorana_operator_dict(
    fermion_op: FermionOperator,
    atol: float = 1e-12,
) -> dict[tuple[int, ...], complex]:
    """Expand a ``FermionOperator`` into canonical Majorana monomials."""
    majorana_terms = defaultdict(complex)

    for term, coeff in fermion_op.terms.items():
        for label, c in fermion_term_to_majorana_terms(term, coeff).items():
            majorana_terms[label] += c

    return {label: coeff for label, coeff in majorana_terms.items() if abs(coeff) > atol}
