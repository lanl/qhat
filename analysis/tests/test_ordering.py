"""
Tests for Pauli-string ordering behavior.

Tests include:
- Default ordering preserves input order and returns a dictionary.
- Empty input returns an empty dictionary.
- group_evolve_xyz correctly handles Y-only Pauli terms.
- group_evolve_xyz safely handles identity-only Pauli terms.
- group_evolve_xyz rejects mixed Pauli terms.
"""

import pytest
from analysis.ordering import reorder_paulis


def test_none_ordering_returns_dict_and_preserves_order():
    pauli_strings = {"ZI": 1.0, "XI": 2.0, "IZ": 3.0}

    ordered = reorder_paulis(pauli_strings, None)

    assert isinstance(ordered, dict)
    assert list(ordered.items()) == list(pauli_strings.items())


def test_empty_ordering_input_returns_empty_dict():
    ordered = reorder_paulis({}, None)

    assert isinstance(ordered, dict)
    assert ordered == {}


def test_group_evolve_xyz_handles_y_only_terms():
    pauli_strings = {"ZII": 3.0, "YII": 2.0, "XII": 1.0}

    ordered = reorder_paulis(pauli_strings, "group_evolve_xyz")

    assert isinstance(ordered, dict)
    assert list(ordered.keys()) == ["XII", "YII", "ZII"]


def test_group_evolve_xyz_handles_identity_only_terms():
    pauli_strings = {"III": 0.5, "YII": 1.0, "ZII": 2.0}

    ordered = reorder_paulis(pauli_strings, "group_evolve_xyz")

    assert list(ordered.keys()) == ["III", "YII", "ZII"]


def test_group_evolve_xyz_rejects_mixed_pauli_terms():
    pauli_strings = {"XYZ": 1.0}

    with pytest.raises(Exception, match="group_evolve_xyz can only be used"):
        reorder_paulis(pauli_strings, "group_evolve_xyz")
