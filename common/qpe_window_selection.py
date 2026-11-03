"""
Choosing the QPE window state by cost: order candidate windows' multi-run plans by
`controlled_unitary_cost` (see `qhat.common.multi_run_qpe_sizing`), cheapest first.

Ties go to the window that is cheaper to prepare, in `TIE_BREAK_ORDER`:

- "rectangular": one Hadamard per qubit;
- "bandwidth-optimized" and "truncated dpss": the same small circuit (2K + 1 rotations and a QFT),
  with the bandwidth-optimized taper failing no more often;
- "dpss": an exact state preparation on all phase qubits.

Windows with no plan come last, in the same order.

Example: plans costing rectangular 1178, dpss 1778, bandwidth-optimized 1905 and truncated dpss
1905 are ordered rectangular, dpss, bandwidth-optimized, truncated dpss.
"""

from typing import Any, Mapping, Optional, Sequence, Tuple

TIE_BREAK_ORDER = ("rectangular", "bandwidth-optimized", "truncated dpss", "dpss")


def rank_windows(plans: Mapping[str, Optional[Any]]) -> Sequence[Tuple[str, Optional[Any]]]:
    """`(window name, plan)` pairs, cheapest first; `plans` maps names in `TIE_BREAK_ORDER` to a
    sizing result with a `controlled_unitary_cost`, or None if that window has no plan."""
    unknown = set(plans) - set(TIE_BREAK_ORDER)
    if unknown:
        raise ValueError(f"Unknown window names {sorted(unknown)}; expected names from "
                         f"{TIE_BREAK_ORDER}.")

    def key(item):
        name, plan = item
        cost = float("inf") if plan is None else plan.controlled_unitary_cost
        return (plan is None, cost, TIE_BREAK_ORDER.index(name))

    return sorted(plans.items(), key=key)
