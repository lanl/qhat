"""
Light unit tests for the QHAT -> prqpe adapter (qhat.analysis.prqpe_adapter).

These tests deliberately import ONLY the adapter, the configuration types,
openfermion, and the prqpe core.  They do NOT import hamiltonian.py or
analysis.py, so they run without the heavy qhat stack (pyLIQTR, qualtran,
pyscf, ...).  A small duck-typed stub stands in for a real QHAT Hamiltonian.

Coverage:
- _qubit_operator_from_qhat drops the identity term and keeps the rest.
- resource_estimation_prqpe returns a well-formed, positive, TOML-serializable
  result dict.
- _resolve_C_gs tiering: explicit value > cgs_rule > small-system auto >
  ValueError when the system is too large and nothing is configured.
- A missing Hamiltonian raises ValueError.
"""

import pytest
import tomlkit
from types import SimpleNamespace

from openfermion import QubitOperator

from qhat.analysis.prqpe_adapter import (
    _qubit_operator_from_qhat,
    _resolve_C_gs,
    resource_estimation_prqpe,
)
from qhat.prqpe.representation import build_representation, sorted_weights
from qhat.prqpe.cgs import cgs_rule


# ==================================================================================
# Shared stub + fixtures
# ==================================================================================

class _StubHam:
    """Duck-typed stand-in for a QHAT Hamiltonian (avoids the heavy stack)."""

    def __init__(self, terms, nq):
        self._t, self._nq = terms, nq

    def get_all_pauli_strings(self, return_as="tuples"):
        assert return_as == "tuples"
        return self._t

    def num_qubits(self):
        return self._nq


# A small Hermitian two-qubit example with an explicit identity term.
SAMPLE_TERMS = {
    (): 0.7,
    ((0, 'Z'),): 1.0,
    ((1, 'Z'),): -0.5,
    ((0, 'X'), (1, 'X')): 0.3,
}


@pytest.fixture
def sample_ham():
    """A small two-qubit Hermitian stub Hamiltonian."""
    return _StubHam(dict(SAMPLE_TERMS), 2)


@pytest.fixture
def config_analysis():
    """A minimal analysis config selecting prqpe with small-system auto C_gs."""
    return SimpleNamespace(
        resource_estimator="prqpe",
        prqpe_C_gs=None,
        prqpe_cgs_rule=None,
        prqpe_cgs_max_qubits=14,
        prqpe_target_precision=0.0016,
    )


@pytest.fixture
def algorithm():
    """A minimal prqpe algorithm carrier (duck-typed)."""
    return SimpleNamespace(
        method="partially_randomized",
        overlap=1.0,
        xi=None,
        randomizer="rte",
        commuting_group_size=None,
    )


# ==================================================================================
# Test: _qubit_operator_from_qhat
# ==================================================================================

class TestQubitOperatorFromQhat:
    """Test the QHAT-Hamiltonian -> openfermion QubitOperator conversion."""

    def test_drops_identity_keeps_rest(self, sample_ham):
        """The identity (empty-tuple) term is dropped; all others survive."""
        qop = _qubit_operator_from_qhat(sample_ham)

        assert isinstance(qop, QubitOperator)
        # Identity term must be absent (it carries no resource cost).
        assert () not in qop.terms
        # The three non-identity terms must be present with their coefficients.
        assert len(qop.terms) == 3
        assert qop.terms[((0, 'Z'),)] == 1.0
        assert qop.terms[((1, 'Z'),)] == -0.5
        assert qop.terms[((0, 'X'), (1, 'X'))] == 0.3

    def test_no_identity_term_is_noop(self):
        """A Hamiltonian without an identity term converts all terms."""
        ham = _StubHam({((0, 'X'),): 0.5, ((1, 'Y'),): 0.25}, 2)
        qop = _qubit_operator_from_qhat(ham)

        assert () not in qop.terms
        assert len(qop.terms) == 2


# ==================================================================================
# Test: resource_estimation_prqpe
# ==================================================================================

class TestResourceEstimationPrqpe:
    """Test the top-level adapter entry point."""

    def test_returns_well_formed_dict(self, config_analysis, algorithm, sample_ham):
        """A successful estimate is a flat dict with positive, typed fields."""
        out = resource_estimation_prqpe(config_analysis, algorithm, hamiltonian=sample_ham)

        assert out["estimator"] == "prqpe"

        assert isinstance(out["toffoli_count"], float)
        assert out["toffoli_count"] > 0

        assert isinstance(out["logical_qubits"], int)
        assert out["logical_qubits"] > 0

        assert isinstance(out["C_gs"], float)
        assert out["C_gs"] > 0

        assert isinstance(out["metadata"], dict)

    def test_result_is_toml_serializable(self, config_analysis, algorithm, sample_ham):
        """The result dict round-trips through tomlkit (wrapped in a table)."""
        out = resource_estimation_prqpe(config_analysis, algorithm, hamiltonian=sample_ham)

        document = tomlkit.document()
        document["resource_estimates"] = out
        rendered = tomlkit.dumps(document)
        reloaded = tomlkit.loads(rendered)

        assert reloaded["resource_estimates"]["estimator"] == "prqpe"
        assert reloaded["resource_estimates"]["logical_qubits"] == out["logical_qubits"]
        assert reloaded["resource_estimates"]["toffoli_count"] == pytest.approx(
            out["toffoli_count"]
        )

    def test_missing_hamiltonian_raises(self, config_analysis, algorithm):
        """A missing Hamiltonian is a usage error."""
        with pytest.raises(ValueError, match="requires the Hamiltonian"):
            resource_estimation_prqpe(config_analysis, algorithm, hamiltonian=None)


# ==================================================================================
# Test: _resolve_C_gs tiering
# ==================================================================================

class TestResolveCgs:
    """Test the tiered ground-state Trotter constant resolution."""

    def test_explicit_value_wins(self):
        """An explicit prqpe_C_gs takes precedence over everything else."""
        config = SimpleNamespace(
            prqpe_C_gs=3.14,
            prqpe_cgs_rule=(0.5, 2.0),  # present, but must be ignored
            prqpe_cgs_max_qubits=14,
        )
        # pauli_ham is unused on this path, so None is fine.
        assert _resolve_C_gs(None, 2, 1.0, config) == 3.14

    def test_cgs_rule_path(self):
        """With no explicit value, the (A, b) rule yields cgs_rule(lam, A, b)."""
        A, b, lam = 0.5, 2.0, 3.0
        config = SimpleNamespace(
            prqpe_C_gs=None,
            prqpe_cgs_rule=(A, b),
            prqpe_cgs_max_qubits=14,
        )
        result = _resolve_C_gs(None, 2, lam, config)
        assert result == pytest.approx(cgs_rule(lam, A, b))
        assert result == pytest.approx(A * lam ** b)

    def test_small_system_auto_path(self, sample_ham):
        """A small system with nothing set diagonalizes for a positive C_gs."""
        pauli_ham = build_representation(_qubit_operator_from_qhat(sample_ham))
        lam = float(sorted_weights(pauli_ham).sum())
        config = SimpleNamespace(
            prqpe_C_gs=None,
            prqpe_cgs_rule=None,
            prqpe_cgs_max_qubits=14,
        )
        result = _resolve_C_gs(pauli_ham, sample_ham.num_qubits(), lam, config)
        assert isinstance(result, float)
        assert result > 0

    def test_too_large_without_config_raises(self, sample_ham):
        """Beyond prqpe_cgs_max_qubits, with nothing else set, raises ValueError."""
        pauli_ham = build_representation(_qubit_operator_from_qhat(sample_ham))
        config = SimpleNamespace(
            prqpe_C_gs=None,
            prqpe_cgs_rule=None,
            prqpe_cgs_max_qubits=1,  # below the 2-qubit system size
        )
        with pytest.raises(ValueError, match="C_gs required"):
            _resolve_C_gs(pauli_ham, 2, 1.0, config)
