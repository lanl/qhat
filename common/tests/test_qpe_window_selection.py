"""
Tests for ordering QPE windows by cost (`qhat.common.qpe_window_selection`).
"""

import itertools

import attrs
import pytest

from qhat.common.qpe_window_selection import rank_windows, TIE_BREAK_ORDER


@attrs.frozen
class Plan:
    controlled_unitary_cost: int


def names(ranked):
    return [name for name, _ in ranked]


class TestRankWindows:

    def test_cheapest_first(self):
        plans = {"dpss": Plan(1778), "rectangular": Plan(1178),
                 "truncated dpss": Plan(1905), "bandwidth-optimized": Plan(1905)}
        assert names(rank_windows(plans)) == [
            "rectangular", "dpss", "bandwidth-optimized", "truncated dpss"]

    @pytest.mark.parametrize("order", list(itertools.permutations(TIE_BREAK_ORDER)))
    def test_ties_follow_preparation_cost(self, order):
        # Whatever order the candidates arrive in, equal costs are ordered by TIE_BREAK_ORDER.
        plans = {name: Plan(100) for name in order}
        assert names(rank_windows(plans)) == list(TIE_BREAK_ORDER)

    def test_tie_below_cheaper_windows(self):
        plans = {"rectangular": Plan(500), "dpss": Plan(300), "bandwidth-optimized": Plan(300),
                 "truncated dpss": Plan(300)}
        assert names(rank_windows(plans)) == [
            "bandwidth-optimized", "truncated dpss", "dpss", "rectangular"]

    def test_no_plan_last_in_tie_break_order(self):
        plans = {"dpss": None, "rectangular": None, "truncated dpss": Plan(7),
                 "bandwidth-optimized": Plan(9)}
        ranked = rank_windows(plans)
        assert names(ranked) == ["truncated dpss", "bandwidth-optimized", "rectangular", "dpss"]
        assert ranked[2][1] is None and ranked[3][1] is None

    def test_subset_of_windows(self):
        assert names(rank_windows({"dpss": Plan(3), "truncated dpss": Plan(3)})) == [
            "truncated dpss", "dpss"]

    def test_unknown_window_raises(self):
        with pytest.raises(ValueError, match="Unknown window names"):
            rank_windows({"kaiser": Plan(1)})
