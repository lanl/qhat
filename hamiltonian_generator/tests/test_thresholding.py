"""Regression tests for QHAT coefficient thresholds."""

import pickle

import numpy as np
import pytest
from openfermion import InteractionOperator, MolecularData, QubitOperator
from openfermion.chem.molecular_data import (
    spinorb_from_spatial as openfermion_spinorb_from_spatial,
)

from qhat.hamiltonian_generator import hamgen
from qhat.hamiltonian_generator.hamgen_types import HamiltonianConfiguration
from qhat.hamiltonian_generator.thresholding import (
    DEFAULT_COEFFICIENT_THRESHOLD,
    spinorb_from_spatial,
)


class _ActiveSpaceState:
    def __init__(self, threshold=DEFAULT_COEFFICIENT_THRESHOLD):
        self.config_hamiltonian = HamiltonianConfiguration()
        self.config_hamiltonian.num_active_occupied = 4
        self.config_hamiltonian.num_active_vacant = 6
        self.config_hamiltonian.coefficient_threshold = threshold

    def log(self, *args, **kwargs):
        pass

    def log_verbose(self, *args, **kwargs):
        pass


class _PassthroughPaths:
    def get_cache_path(self, filename):
        return filename

    def get_output_path(self, filename):
        return filename


class _CacheState:
    def __init__(self, cache_path, threshold):
        self.cache_path = cache_path
        self.config_general = _PassthroughPaths()
        self.config_hamiltonian = HamiltonianConfiguration()
        self.config_hamiltonian.coefficient_threshold = threshold

    def filename_ham1(self):
        return str(self.cache_path.with_name("hartree-fock.pickle"))

    def filename_ham2(self):
        return str(self.cache_path)

    def log(self, *args, **kwargs):
        pass


def _interaction_operator(constant, threshold=None):
    operator = InteractionOperator(
        constant,
        np.diag([constant, 0.0]),
        np.zeros((2, 2, 2, 2)),
    )

    operator.hf_energy = -1.0
    operator.asmeta = {
        "n_frz_occ_so": 0,
        "n_act_occ_so": 1,
        "n_act_vac_so": 1,
        "n_frz_vac_so": 0,
    }
    operator.basis = "test-basis"
    operator.separation = 1.0
    operator.hf_time = 0.0
    operator.as_time = 0.0

    if threshold is not None:
        operator.coefficient_threshold = threshold

    return operator


def _synthetic_active_space_molecule():
    """
    Build a deterministic MolecularData object without PySCF or pickle files.

    The synthetic molecule has:

    - 6 electrons
    - 6 spatial orbitals
    - 12 spin orbitals
    - 1 frozen occupied spatial orbital
    - 5 active spatial orbitals

    This reproduces the active-space ranges used by the original lithium
    fixture without requiring a generated binary artifact.
    """
    threshold = DEFAULT_COEFFICIENT_THRESHOLD
    n_spatial_orbitals = 6

    molecule = MolecularData(
        geometry=[
            ("H", (0.0, 0.0, float(index)))
            for index in range(n_spatial_orbitals)
        ],
        basis="sto-3g",
        multiplicity=1,
        charge=0,
    )

    molecule.n_orbitals = n_spatial_orbitals
    molecule.n_qubits = 2 * n_spatial_orbitals
    molecule.nuclear_repulsion = 1.0
    molecule.hf_energy = -1.0

    molecule.basis = "sto-3g"
    molecule.separation = None
    molecule.hf_time = 0.0

    one_body_integrals = np.zeros(
        (n_spatial_orbitals, n_spatial_orbitals),
    )

    # All nonzero entries are inside the active spatial-orbital range [1, 6).
    #
    # Include coefficients below, exactly at, and above the default threshold.
    one_body_integrals[1, 1] = threshold / 2.0
    one_body_integrals[1, 2] = threshold
    one_body_integrals[2, 1] = threshold
    one_body_integrals[2, 2] = -2.0 * threshold

    two_body_integrals = np.zeros(
        (
            n_spatial_orbitals,
            n_spatial_orbitals,
            n_spatial_orbitals,
            n_spatial_orbitals,
        ),
    )

    two_body_integrals[1, 1, 1, 1] = threshold / 2.0
    two_body_integrals[1, 1, 1, 2] = threshold
    two_body_integrals[2, 1, 1, 2] = -2.0 * threshold

    molecule.one_body_integrals = one_body_integrals
    molecule.two_body_integrals = two_body_integrals

    return molecule


def test_configuration_uses_openfermion_default_threshold():
    config = HamiltonianConfiguration()

    assert (
        config.coefficient_threshold
        == DEFAULT_COEFFICIENT_THRESHOLD
        == 1.0e-8
    )


@pytest.mark.parametrize(
    "value",
    [
        0,
        1.0e-12,
        np.float64(1.0e-8),
    ],
)
def test_configuration_accepts_non_negative_real_thresholds(value):
    config = HamiltonianConfiguration()
    config.coefficient_threshold = value

    assert config.coefficient_threshold == float(value)


@pytest.mark.parametrize(
    ("value", "exception"),
    [
        (-1.0e-8, ValueError),
        (float("nan"), ValueError),
        (float("inf"), ValueError),
        (-float("inf"), ValueError),
        (True, TypeError),
        ("1e-8", TypeError),
        (None, TypeError),
    ],
)
def test_configuration_rejects_invalid_thresholds(value, exception):
    config = HamiltonianConfiguration()

    with pytest.raises(exception):
        config.coefficient_threshold = value


def test_spinorb_default_parity_strict_boundary_and_zero_threshold():
    threshold = DEFAULT_COEFFICIENT_THRESHOLD
    below_threshold = threshold / 2.0

    one_body = np.array(
        [
            [below_threshold, threshold],
            [-2.0 * threshold, 0.0],
        ]
    )

    two_body = np.zeros((2, 2, 2, 2))
    two_body[0, 0, 0, 0] = below_threshold
    two_body[0, 0, 0, 1] = threshold
    two_body[1, 0, 0, 1] = -2.0 * threshold

    expected_one, expected_two = openfermion_spinorb_from_spatial(
        one_body,
        two_body,
    )

    default_one, default_two = spinorb_from_spatial(
        one_body,
        two_body,
    )

    zero_one, zero_two = spinorb_from_spatial(
        one_body,
        two_body,
        0.0,
    )

    np.testing.assert_array_equal(
        default_one,
        expected_one,
    )

    np.testing.assert_array_equal(
        default_two,
        expected_two,
    )

    # Coefficients strictly below 1e-8 are removed.
    assert default_one[0, 0] == 0.0
    assert default_two[0, 1, 1, 0] == 0.0

    # Coefficients exactly equal to 1e-8 are retained.
    assert default_one[0, 2] == threshold
    assert default_two[0, 1, 1, 2] == threshold

    # A zero threshold retains the smaller coefficients.
    assert zero_one[0, 0] == below_threshold
    assert zero_two[0, 1, 1, 0] == below_threshold


def test_active_space_matches_openfermion_and_respects_zero_threshold():
    molecule = _synthetic_active_space_molecule()

    expected_operator = molecule.get_molecular_hamiltonian(
        occupied_indices=range(1),
        active_indices=range(1, 6),
    )

    default_operator = hamgen.apply_active_space(
        _ActiveSpaceState(),
        molecule,
    )

    assert default_operator.constant == expected_operator.constant

    np.testing.assert_array_equal(
        default_operator.one_body_tensor,
        expected_operator.one_body_tensor,
    )

    np.testing.assert_array_equal(
        default_operator.two_body_tensor,
        expected_operator.two_body_tensor,
    )

    assert (
        default_operator.coefficient_threshold
        == DEFAULT_COEFFICIENT_THRESHOLD
    )

    zero_operator = hamgen.apply_active_space(
        _ActiveSpaceState(0.0),
        molecule,
    )

    assert zero_operator.coefficient_threshold == 0.0

    assert np.count_nonzero(
        zero_operator.one_body_tensor
    ) > np.count_nonzero(
        default_operator.one_body_tensor
    )

    assert np.count_nonzero(
        zero_operator.two_body_tensor
    ) > np.count_nonzero(
        default_operator.two_body_tensor
    )


def test_mapping_and_metadata_record_configured_threshold():
    state = _ActiveSpaceState(0.0)
    state.metadata = {}

    active_operator = _interaction_operator(
        1.0,
        threshold=0.0,
    )

    mapped_operator = hamgen.map_fermions_to_qubits(
        state,
        active_operator,
    )

    hamgen.compute_metadata(
        state,
        mapped_operator,
    )

    assert (
        state.metadata[
            "spin-orbital coefficient threshold (Hartrees)"
        ]
        == 0.0
    )

    assert "Pauli-string coefficient threshold (Hartrees)" not in state.metadata


@pytest.mark.parametrize(
    ("mapping_setting", "mapping_name", "mapper_name"),
    [
        (
            "Jordan-Wigner",
            "jordan-wigner",
            "jordan_wigner",
        ),
        (
            "Bravyi-Kitaev",
            "bravyi-kitaev",
            "bravyi_kitaev",
        ),
    ],
)
def test_mapping_does_not_apply_additional_pauli_string_threshold(
    monkeypatch,
    mapping_setting,
    mapping_name,
    mapper_name,
):
    threshold = DEFAULT_COEFFICIENT_THRESHOLD

    # Assign terms directly so OpenFermion does not remove the synthetic
    # below-threshold coefficient before QHAT receives the mapped operator.
    mapped_operator = QubitOperator()
    mapped_operator.terms = {
        tuple(): 1.0,
        ((0, "X"),): 0.5 * threshold,
        ((0, "Y"),): threshold,
        ((0, "Z"),): 2.0 * threshold,
    }

    def fake_mapping(active_operator):
        return mapped_operator

    monkeypatch.setattr(
        hamgen,
        mapper_name,
        fake_mapping,
    )

    state = _ActiveSpaceState()
    state.metadata = {}
    state.config_hamiltonian.f2q_mapping = mapping_setting

    active_operator = _interaction_operator(
        1.0,
        threshold=DEFAULT_COEFFICIENT_THRESHOLD,
    )

    result = hamgen.map_fermions_to_qubits(
        state,
        active_operator,
    )

    # QHAT preserves all terms returned by OpenFermion, including coefficients
    # below the tensor-construction threshold.
    assert result.terms[((0, "X"),)] == 0.5 * threshold

    # A coefficient exactly equal to the threshold is retained.
    assert result.terms[((0, "Y"),)] == threshold

    # A coefficient above the threshold is retained unchanged.
    assert result.terms[((0, "Z"),)] == 2.0 * threshold

    # The identity term is also retained unchanged.
    assert result.terms[tuple()] == 1.0

    assert result.f2q_mapping == mapping_name
    assert not hasattr(result, "pauli_coefficient_threshold")
    assert len(result.terms) == 4

    hamgen.compute_metadata(
        state,
        result,
    )

    assert "Pauli-string coefficient threshold (Hartrees)" not in state.metadata

    assert (
        state.metadata[
            "number of terms in sum of Pauli strings"
        ]
        == 4
    )


def test_active_space_cache_rejects_mismatch_and_rebuilds_after_removal(
    tmp_path,
    monkeypatch,
):
    cache_path = tmp_path / "active-space.pickle"

    legacy_operator = _interaction_operator(1.0)

    with open(cache_path, "wb") as file:
        pickle.dump(legacy_operator, file)

    calls = []

    recomputed_operator = _interaction_operator(
        2.0,
        threshold=0.0,
    )

    def fake_get_ham1(state):
        calls.append("get_ham1")
        return object()

    def fake_apply_active_space(state, molecule):
        calls.append("apply_active_space")
        return recomputed_operator

    monkeypatch.setattr(
        hamgen,
        "get_ham1",
        fake_get_ham1,
    )

    monkeypatch.setattr(
        hamgen,
        "apply_active_space",
        fake_apply_active_space,
    )

    default_state = _CacheState(
        cache_path,
        DEFAULT_COEFFICIENT_THRESHOLD,
    )

    loaded_legacy = hamgen.get_ham2(default_state)

    assert loaded_legacy.constant == legacy_operator.constant
    assert calls == []

    zero_state = _CacheState(
        cache_path,
        0.0,
    )

    with pytest.raises(
        ValueError,
        match="remove the cache and rerun",
    ):
        hamgen.get_ham2(zero_state)

    assert calls == []

    cache_path.unlink()

    rebuilt = hamgen.get_ham2(zero_state)

    assert rebuilt is recomputed_operator
    assert calls == [
        "get_ham1",
        "apply_active_space",
    ]

    with open(cache_path, "rb") as file:
        persisted = pickle.load(file)

    assert persisted.coefficient_threshold == 0.0
    assert persisted.constant == recomputed_operator.constant

    tensor_path = tmp_path / "active-space.tensors.npz"

    with np.load(tensor_path) as tensors:
        assert set(tensors.files) == {
            "constant",
            "one_body",
            "two_body",
        }

        assert tensors["constant"] == recomputed_operator.constant

        np.testing.assert_array_equal(
            tensors["one_body"],
            recomputed_operator.one_body_tensor,
        )

        np.testing.assert_array_equal(
            tensors["two_body"],
            recomputed_operator.two_body_tensor,
        )

    loaded_match = hamgen.get_ham2(zero_state)

    assert loaded_match.constant == recomputed_operator.constant
    assert loaded_match.coefficient_threshold == 0.0

    assert calls == [
        "get_ham1",
        "apply_active_space",
    ]
