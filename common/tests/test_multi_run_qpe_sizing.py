"""
Tests for the window-independent multi-run sizing (`qhat.common.multi_run_qpe_sizing`), in the
regime of textbook QPE: per-run failures up to 0.2 and quantiles up to about N/2. The tapered
regime (small u, K of 1 to 6) is covered by `test_tapered_qpe_sizing.py`.
"""

import math

import numpy as np
import pytest

from qhat.common import tapered_qpe_sizing, textbook_qpe_sizing
from qhat.common.tapered_window_state import (
    BandwidthOptimizedWindowState, DPSSWindowState, TruncatedDPSSWindowState)
from qhat.common.trotter_flattened import Trotterization

from qhat.common import multi_run_qpe_sizing
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.multi_run_qpe_sizing import (
    check_targets,
    controlled_unitary_cost,
    max_runs_within,
    fewest_runs_for_failure,
    worst_failure,
)


def direct_failure(n, k, under, over):
    """P(the k-th smallest of n runs is an under- or overestimate), summed term by term."""
    def at_least(j_min, q):
        return sum(math.comb(n, j) * q**j * (1 - q)**(n - j) for j in range(j_min, n + 1))
    return at_least(k, under) + at_least(n - k + 1, over)


def brute_force(u, hit, target, force, max_runs):
    for n in range(1, max_runs + 1):
        quantiles = [1] if force else range(1, n + 1)
        failures = {k: worst_failure(n, k, u, hit) for k in quantiles}
        if min(failures.values()) <= target:
            return n, min(failures, key=failures.get)
    return None


class TestWorstFailure:
    @pytest.mark.parametrize("n, k", [(1, 1), (9, 5), (38, 9), (60, 31)])
    @pytest.mark.parametrize("u, overlap", [(0.06, 0.5), (0.14, 1.0), (0.2, 0.3)])
    def test_matches_direct_sum(self, n, k, u, overlap):
        hit = overlap * (1 - u)
        expected = max(direct_failure(n, k, 0.0, 1 - hit), direct_failure(n, k, u, 1 - u - hit))
        assert worst_failure(n, k, u, hit) == pytest.approx(expected, rel=1e-10, abs=1e-15)

    @pytest.mark.parametrize("n, k", [(13, 7), (38, 9), (82, 18)])
    def test_worst_split_is_an_endpoint(self, n, k):
        u, overlap = 0.13, 0.6
        hit = overlap * (1 - u)
        bound = worst_failure(n, k, u, hit)
        for under in np.linspace(0, u, 41):
            assert direct_failure(n, k, under, 1 - hit - under) <= bound * (1 + 1e-12)


class TestFewestRunsForFailure:
    # Answers range from N = 1 to about 330 and K up to about N/2; in some, forcing K = 1 is
    # impossible.
    @pytest.mark.parametrize("u", [1e-3, 1e-2, 0.05, 0.13, 0.2])
    @pytest.mark.parametrize("overlap", [1.0, 0.7, 0.35])
    @pytest.mark.parametrize("target", [0.1, 1e-3, 1e-6])
    @pytest.mark.parametrize("force", [False, True])
    def test_matches_exhaustive_search(self, u, overlap, target, force):
        hit = overlap * (1 - u)
        max_runs = 400
        expected = brute_force(u, hit, target, force, max_runs)
        got = fewest_runs_for_failure(u, hit, target, force, max_runs)
        if expected is None:
            assert got is None
            return
        assert got[0] == expected[0]
        assert worst_failure(*got, u, hit) == pytest.approx(
            worst_failure(*expected, u, hit), rel=1e-9)
        assert worst_failure(*got, u, hit) <= target

    @pytest.mark.parametrize("u, overlap", [(0.3, 0.4), (0.5, 0.8), (0.25, 0.25 / 0.75)])
    def test_underestimates_as_likely_as_hits_never_work(self, u, overlap):
        hit = overlap * (1 - u)
        assert u >= hit
        # One of the two extreme splits fails with probability at least 1/2 for every N and K.
        for n in [1, 7, 40, 150]:
            for k in range(1, n + 1, max(1, n // 10)):
                assert worst_failure(n, k, u, hit) >= 0.5 - 1e-12
        assert fewest_runs_for_failure(u, hit, 0.49, False, 10**9) is None

    def test_huge_quantile_is_fast(self):
        # Underestimates nearly as likely as hits: the only plans have about 4e8 runs and K about
        # 4e5. Stepping K one at a time took 90 s; the search must still finish quickly.
        u = 1 / (2 * (2**9 - 2))
        hit = 1e-3 * (1 - u)
        got = fewest_runs_for_failure(u, hit, 1e-9, False, 4 * 10**8)
        assert got is None or worst_failure(*got, u, hit) <= 1e-9

    @pytest.mark.parametrize("u, hit, target, max_runs, expected", [
        # K = 801 works, 802 does not, 803 does: the K search must not skip past 801.
        (0.3563216679642024, 0.4384013708679965, 8.105632705620942e-05, 17705, (2017, 801)),
    ])
    def test_quantile_jitter(self, u, hit, target, max_runs, expected):
        assert fewest_runs_for_failure(u, hit, target, False, max_runs)[0] == expected[0]

    def test_respects_max_runs(self):
        u, hit = 0.05, 0.3
        unlimited = fewest_runs_for_failure(u, hit, 1e-6, False, 10**6)
        assert fewest_runs_for_failure(u, hit, 1e-6, False, unlimited[0]) == unlimited
        assert fewest_runs_for_failure(u, hit, 1e-6, False, unlimited[0] - 1) is None


class TestCheckTargets:
    @pytest.mark.parametrize("probability_of_failure", [0.0, 1.0, -1e-3])
    def test_rejects_probability_of_failure(self, probability_of_failure):
        with pytest.raises(ValueError, match="probability of failure"):
            check_targets(probability_of_failure, 0.5)

    @pytest.mark.parametrize("overlap", [0.0, -0.1, 1.5])
    def test_rejects_overlap(self, overlap):
        with pytest.raises(ValueError, match="overlap"):
            check_targets(1e-3, overlap)


class TestControlledUnitaryCost:
    """The plan cost, N (2^m - 1) controlled-U invocations, defined in one place."""

    @pytest.mark.parametrize("num_phase_qubits, num_runs, cost",
                             [(1, 1, 1), (5, 38, 1178), (7, 14, 1778), (10, 3, 3069)])
    def test_values(self, num_phase_qubits, num_runs, cost):
        assert controlled_unitary_cost(num_phase_qubits, num_runs) == cost

    @pytest.mark.parametrize("window", [
        None, DPSSWindowState(2, 3), TruncatedDPSSWindowState(2, 3),
        BandwidthOptimizedWindowState(2, 3)], ids=["rectangular", "dpss", "truncated", "bo"])
    def test_matches_flexible_qpe_call_graph(self, window):
        # The Trotterization fast-forwards U^(2^j) as 2^j times the steps, so the call graph's
        # controlled Trotterizations add up to the number of controlled U.
        u = Trotterization.from_method(pauli_terms=[("XZ", 0.5), ("ZY", 0.3)],
                                       method="first order", time=0.4, num_steps=3)
        qpe = FlexibleQPE(u, 5) if window is None else FlexibleQPE(u, 5, window)
        total = sum(count * bloq.subbloq.num_steps // u.num_steps
                    for bloq, count in qpe.build_call_graph(None)
                    if isinstance(getattr(bloq, "subbloq", None), Trotterization))
        assert total == controlled_unitary_cost(5)

    @pytest.mark.parametrize("limit, num_phase_qubits, runs",
                             [(1178, 5, 38), (1177, 5, 37), (30, 5, 0),
                              (1e18, 10, 977517106549364)])
    def test_max_runs_within(self, limit, num_phase_qubits, runs):
        assert max_runs_within(limit, num_phase_qubits) == runs

    def test_sizing_uses_the_one_cost_function(self, monkeypatch):
        # Make one extra phase qubit cost almost nothing: if every sizing path uses
        # controlled_unitary_cost, both searches then pick the most qubits allowed.
        def flat_cost(num_phase_qubits, num_runs=1):
            return num_runs * (1000 + num_phase_qubits)

        textbook = textbook_qpe_sizing.size_textbook_qpe(0.2 / 1.7, 1e-3, 0.5)
        tapered = tapered_qpe_sizing.size_tapered_qpe(DPSSWindowState, 0.2 / 1.7, 1e-3, 0.5)
        monkeypatch.setattr(multi_run_qpe_sizing, "controlled_unitary_cost", flat_cost)
        textbook_flat = textbook_qpe_sizing.size_textbook_qpe(0.2 / 1.7, 1e-3, 0.5)
        tapered_flat = tapered_qpe_sizing.size_tapered_qpe(DPSSWindowState, 0.2 / 1.7, 1e-3, 0.5)
        assert textbook_flat.num_phase_qubits > textbook.num_phase_qubits
        assert textbook_flat.controlled_unitary_cost == flat_cost(textbook_flat.num_phase_qubits,
                                                                  textbook_flat.num_runs)
        assert tapered_flat.num_confidence_qubits > tapered.num_confidence_qubits
        assert tapered_flat.controlled_unitary_cost == flat_cost(tapered_flat.num_phase_qubits,
                                                                 tapered_flat.num_runs)
