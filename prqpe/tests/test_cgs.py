"""First-principles tests for the ground-state Trotter constant estimator.

The benchmark is a small transverse-field Ising chain, whose coupling and field
terms do not commute, so the symmetric second-order Trotter step carries a
nonzero ground-state bias.  The expected behavior is derived from the formula
order, not from any tabulated value: the bias scales as ``delta^2``, so the
fitted exponent must be close to two and the fitted constant positive and
finite.  Everything here is ``numpy``-only.
"""

import numpy as np
import pytest
import scipy.linalg as sla

from qhat.prqpe.cgs import cgs_rule, estimate_cgs, fit_cgs
from qhat.prqpe.cgs.trotter_matrix import (
    effective_ground_energy,
    ground_state,
    trot2_matrix,
)
from qhat.prqpe.representation.pauli import PauliHamiltonian


def transverse_field_ising(n_qubits, J=1.0, g=0.6):
    """A transverse-field Ising chain ``H = -J sum ZZ - g sum X`` (open chain).

    The ``ZZ`` couplings and the ``X`` fields do not commute, so the symmetric
    Trotter step has a nonzero ground-state bias.
    """
    xs, zs, cs = [], [], []
    for q in range(n_qubits - 1):
        zz = (1 << (n_qubits - 1 - q)) | (1 << (n_qubits - 1 - (q + 1)))
        xs.append(0)
        zs.append(zz)
        cs.append(-J)
    for q in range(n_qubits):
        xs.append(1 << (n_qubits - 1 - q))
        zs.append(0)
        cs.append(-g)
    return PauliHamiltonian(n_qubits, xs, zs, cs)


# ----------------------------------------------------------------------------
# The bias power law: second-order Trotter gives p ~= 2.
# ----------------------------------------------------------------------------

def test_estimated_exponent_is_about_two():
    # On a clean step-size grid the fitted exponent of |E0 - E0~| = C_gs*delta^p
    # must be close to the second-order value p = 2.
    ham = transverse_field_ising(4)
    e0, psi0 = ground_state(ham)
    deltas = np.geomspace(0.2, 0.02, 7)
    errors = np.array(
        [abs(e0 - effective_ground_energy(trot2_matrix(ham, d), d, psi0)[0])
         for d in deltas]
    )
    _c_gs, p, _res = fit_cgs(errors, deltas)
    assert abs(p - 2.0) < 0.15


def test_estimate_cgs_is_positive_and_finite():
    c_gs = estimate_cgs(transverse_field_ising(4))
    assert np.isfinite(c_gs)
    assert c_gs > 0.0


def test_estimate_cgs_matches_manual_fit():
    # The convenience entry point reproduces the explicit fit on its own grid.
    ham = transverse_field_ising(3)
    e0, psi0 = ground_state(ham)
    from qhat.prqpe.cgs import DEFAULT_DELTAS

    errors = np.array(
        [abs(e0 - effective_ground_energy(trot2_matrix(ham, d), d, psi0)[0])
         for d in DEFAULT_DELTAS]
    )
    c_gs, _p, _res = fit_cgs(errors, DEFAULT_DELTAS)
    assert estimate_cgs(ham) == pytest.approx(c_gs)


# ----------------------------------------------------------------------------
# Trotter-unitary and effective-energy sanity.
# ----------------------------------------------------------------------------

def test_trot2_matrix_is_unitary():
    ham = transverse_field_ising(3)
    U = trot2_matrix(ham, 0.1)
    np.testing.assert_allclose(U.conj().T @ U, np.eye(ham.dim), atol=1e-12)


def test_trot2_matrix_is_time_reversal_symmetric():
    # The symmetric Strang step obeys S_2(-delta) = S_2(delta)^H.
    ham = transverse_field_ising(3)
    for delta in (0.2, 0.1, 0.05):
        forward = trot2_matrix(ham, delta)
        backward = trot2_matrix(ham, -delta)
        np.testing.assert_allclose(forward, backward.conj().T, atol=1e-12)


def test_trot2_matrix_local_error_is_third_order():
    # ||S_2(delta) - exp(-i delta H)||_2 scales as delta^3 for the symmetric step.
    ham = transverse_field_ising(3)
    H = ham.to_dense()
    deltas = np.geomspace(0.2, 0.025, 5)
    errors = np.array(
        [np.linalg.norm(trot2_matrix(ham, d) - sla.expm(-1j * d * H), 2)
         for d in deltas]
    )
    slope = np.polyfit(np.log(deltas), np.log(errors), 1)[0]
    assert abs(slope - 3.0) < 0.2


def test_effective_energy_converges_to_exact_as_step_shrinks():
    # The Trotter bias vanishes with the step, so E0~ -> E0 as delta -> 0 and
    # the selected eigenstate's overlap with the exact ground state -> 1.
    ham = transverse_field_ising(3)
    e0, psi0 = ground_state(ham)
    e_coarse, ov_coarse = effective_ground_energy(trot2_matrix(ham, 0.2), 0.2, psi0)
    e_fine, ov_fine = effective_ground_energy(trot2_matrix(ham, 0.01), 0.01, psi0)
    assert abs(e_fine - e0) < abs(e_coarse - e0)
    assert abs(e_fine - e0) < 1e-4
    assert ov_fine > ov_coarse
    assert ov_fine == pytest.approx(1.0, abs=1e-6)


def test_ground_state_matches_dense_diagonalization():
    ham = transverse_field_ising(4)
    e0, psi0 = ground_state(ham)
    evals = np.linalg.eigvalsh(ham.to_dense())
    assert e0 == pytest.approx(evals[0])
    # psi0 is a unit eigenvector at energy E0.
    Hpsi = ham.to_dense() @ psi0
    assert np.vdot(psi0, psi0) == pytest.approx(1.0)
    np.testing.assert_allclose(Hpsi, e0 * psi0, atol=1e-8)


# ----------------------------------------------------------------------------
# Extrapolation rule.
# ----------------------------------------------------------------------------

def test_cgs_rule_is_the_power_law():
    assert cgs_rule(10.0, 0.5, 1.3) == pytest.approx(0.5 * 10.0 ** 1.3)
    assert cgs_rule(64.0, 0.02, 1.0) == pytest.approx(0.02 * 64.0)


def test_cgs_rule_scales_with_lambda():
    # Doubling lambda multiplies C_gs by 2**b.
    base = cgs_rule(50.0, 0.3, 1.2)
    doubled = cgs_rule(100.0, 0.3, 1.2)
    assert doubled / base == pytest.approx(2.0 ** 1.2)
