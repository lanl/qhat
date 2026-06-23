import numpy as np
import pytest

openfermion = pytest.importorskip("openfermion")

from qhat.trotter_isomorphism.comparison import compare_jw_bk_operator_norm_errors
from qhat.trotter_isomorphism.jw_bk import build_clifford_from_z_map, build_jw_bk_analysis


def test_build_jw_bk_analysis_for_two_mode_hopping_operator():
    FermionOperator = openfermion.FermionOperator
    fermion_op = (
        FermionOperator(((0, 1), (1, 0)), 1.0)
        + FermionOperator(((1, 1), (0, 0)), 1.0)
        + FermionOperator(((0, 1), (0, 0)), 0.25)
        + FermionOperator(((1, 1), (1, 0)), -0.75)
    )

    analysis = build_jw_bk_analysis(
        fermion_op,
        n_qubits=2,
        atol=1e-10,
        remove_identity=True,
        build_matrices=True,
    )

    assert analysis.n_qubits == 2
    assert len(analysis.jw_terms) == len(analysis.bk_terms)
    assert () not in analysis.jw_order

    C_jw_to_bk, _, _ = build_clifford_from_z_map(analysis.jw_to_bk, analysis.n_qubits)
    assert np.allclose(C_jw_to_bk.conj().T @ C_jw_to_bk, np.eye(4))
    assert np.allclose(
        C_jw_to_bk @ analysis.H_jw @ C_jw_to_bk.conj().T,
        analysis.H_bk,
    )


def test_matched_operator_norm_errors_agree():
    FermionOperator = openfermion.FermionOperator
    fermion_op = (
        FermionOperator(((0, 1), (1, 0)), 1.0)
        + FermionOperator(((1, 1), (0, 0)), 1.0)
        + FermionOperator(((0, 1), (0, 0)), 0.25)
        + FermionOperator(((1, 1), (1, 0)), -0.75)
    )

    analysis = build_jw_bk_analysis(
        fermion_op,
        n_qubits=2,
        atol=1e-10,
        remove_identity=True,
        build_matrices=True,
    )
    result = compare_jw_bk_operator_norm_errors(
        H_jw=analysis.H_jw,
        H_bk=analysis.H_bk,
        jw_terms=analysis.jw_terms,
        bk_terms=analysis.bk_terms,
        time=0.4,
        num_steps=2,
        n_qubits=2,
        method=2,
    )

    assert result["error_difference"] < 1e-10
