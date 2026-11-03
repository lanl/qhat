"""
Sizing tapered QPE (see `qhat.common.tapered_window_state`) when the initial state only partly
overlaps the target eigenstate: how many precision and confidence qubits to use, how many times N
to run the circuit, and which estimate to report (the K-th smallest of the N). This implements the
multi-run analysis by Yigit Subasi and Alan Bu (unpublished).

One run returns one phase estimate. With per-run QPE failure u and overlap p between the initial
state and the target eigenstate, the estimate is

- accepted (within 2^-(num_precision_qubits + 1) turns of the target phase) with probability at
  least p (1 - u);
- an underestimate with probability at most u. Nothing lies below the target, so only a failed
  QPE gives an underestimate;
- otherwise an overestimate.

So one run fails with probability at most 1 - p (1 - u). Reporting the K-th smallest of N runs
fails if that estimate is an under- or overestimate; `multi_run_failure` bounds this over every
split of the failure probability allowed above. The window-independent part is in
`qhat.common.multi_run_qpe_sizing`.

Example: phase error 0.2/1.7 turns, probability of failure 1e-3, overlap 0.5, DPSS window: 3
precision and 4 confidence qubits; one run fails with probability up to 0.50, and the 2nd-smallest
of 14 runs fails with probability up to 9.2e-4.

The per-run failure u is `failure_bound_factor * dpss_average_failure_bound`, an average over the
offset of the true phase from the QPE grid. The paper's random phase shift (arXiv:2403.18927,
Sec. VI) would make it a worst-case bound; qhat does not apply that shift.
"""

from typing import Optional, Tuple, Type

import attrs

# COST_CEILING is re-exported for callers of this module.
from qhat.common.multi_run_qpe_sizing import (
    check_targets,
    COST_CEILING,
    fewest_runs_for_failure,
    worst_failure,
)
# The cost is looked up through the module so that it is defined in one place
# (see `multi_run_qpe_sizing.controlled_unitary_cost`).
from qhat.common import multi_run_qpe_sizing
from qhat.common.qpe_window_state import precision_bits
from qhat.common.tapered_window_state import TaperedWindowStateBase, dpss_average_failure_bound

# Confidence qubits tried by `size_tapered_qpe`.
MIN_CONFIDENCE_QUBITS = 2
MAX_CONFIDENCE_QUBITS = 12


@attrs.frozen
class TaperedQPESizing:
    """Qubit counts and run plan. E.g. 3 precision and 4 confidence qubits, run 14 times and report
    the 2nd-smallest estimate (`quantile` 2): `single_run_failure` 0.50, `multi_run_failure`
    9.2e-4 at `overlap` 0.5.
    """

    num_precision_qubits: int
    num_confidence_qubits: int
    num_runs: int
    quantile: int
    single_run_failure: float
    multi_run_failure: float
    overlap: float
    force_lowest_estimate: bool

    @property
    def num_phase_qubits(self) -> int:
        return self.num_precision_qubits + self.num_confidence_qubits

    @property
    def controlled_unitary_cost(self) -> int:
        """`controlled_unitary_cost` of the plan, e.g. 14 runs of 7 phase qubits: 1778."""
        return multi_run_qpe_sizing.controlled_unitary_cost(self.num_phase_qubits,
                                                            self.num_runs)


def split_phase_qubits(phase_error: float, num_phase_qubits: int) -> Tuple[int, int]:
    """(precision, confidence) qubits for a given total: the fewest precision qubits that meet
    `phase_error` (turns), and the rest for confidence. E.g. 0.2/1.7 turns and 9 phase qubits:
    (3, 6).

    At a fixed total each run costs the same, and every precision qubit beyond the fewest takes a
    confidence qubit away, which raises the per-run failure bound and so cannot lower the runs
    needed. Raises a ValueError if fewer than one confidence qubit is left.
    """
    num_precision = precision_bits(phase_error) - 1
    num_confidence = num_phase_qubits - num_precision
    if num_confidence < 1:
        raise ValueError(
            f"{num_phase_qubits} phase qubits cannot meet a phase error of {phase_error:.3g} "
            f"turns with a tapered window: that needs {num_precision} precision qubits plus at "
            "least 1 confidence qubit.")
    return num_precision, num_confidence


def per_run_qpe_failure(window_cls: Type[TaperedWindowStateBase], num_precision_qubits: int,
                        num_confidence_qubits: int) -> float:
    """u = `failure_bound_factor * dpss_average_failure_bound(...)`, at most 1. E.g. the DPSS with
    3 precision and 4 confidence qubits: 2.7e-4."""
    return min(1.0, window_cls.failure_bound_factor
               * dpss_average_failure_bound(num_precision_qubits, num_confidence_qubits))


def multi_run_failure(window_cls: Type[TaperedWindowStateBase], num_precision_qubits: int,
                      num_confidence_qubits: int, num_runs: int, quantile: int,
                      overlap: float) -> float:
    """Bound on the probability that the `quantile`-th smallest of `num_runs` estimates misses the
    target phase by more than 2^-(num_precision_qubits + 1) turns. E.g. the DPSS with 3 + 4
    qubits, 14 runs, quantile 2: 9.2e-4 at overlap 0.5, 6.8e-6 at overlap 1.
    """
    u = per_run_qpe_failure(window_cls, num_precision_qubits, num_confidence_qubits)
    return worst_failure(num_runs, quantile, u, overlap * (1 - u))


def fewest_runs(window_cls: Type[TaperedWindowStateBase], num_precision_qubits: int,
                num_confidence_qubits: int, probability_of_failure: float, overlap: float,
                force_lowest_estimate: bool = False) -> Optional[TaperedQPESizing]:
    """The fewest runs N (and the best K for that N; K = 1 if `force_lowest_estimate`) that meet
    `probability_of_failure` with these qubit counts, or None if none does within
    `COST_CEILING`.
    """
    check_targets(probability_of_failure, overlap)
    return _plan(window_cls, num_precision_qubits, num_confidence_qubits, probability_of_failure,
                 overlap, force_lowest_estimate, COST_CEILING)


def size_tapered_qpe(window_cls: Type[TaperedWindowStateBase], phase_error: float,
                     probability_of_failure: float, overlap: float,
                     force_lowest_estimate: bool = False) -> Optional[TaperedQPESizing]:
    """The cheapest plan meeting `phase_error` (turns) and `probability_of_failure` (for the
    reported estimate of N runs), or None if none costs at most `COST_CEILING`.

    Precision qubits: the fewest with 2^-(num_precision_qubits + 1) <= phase_error. Confidence
    qubits (MIN_CONFIDENCE_QUBITS to MAX_CONFIDENCE_QUBITS), N and K minimize
    `controlled_unitary_cost`; ties go to fewer qubits. K = 1 if `force_lowest_estimate`.
    """
    check_targets(probability_of_failure, overlap)
    num_precision = precision_bits(phase_error) - 1
    best, limit = None, COST_CEILING
    for num_confidence in range(MIN_CONFIDENCE_QUBITS, MAX_CONFIDENCE_QUBITS + 1):
        if multi_run_qpe_sizing.controlled_unitary_cost(num_precision + num_confidence) > limit:
            break
        plan = _plan(window_cls, num_precision, num_confidence, probability_of_failure, overlap,
                     force_lowest_estimate, limit)
        if plan is None:
            continue
        cost = plan.controlled_unitary_cost
        if best is None or cost < limit:
            best, limit = plan, cost
    return best


def _plan(window_cls, num_precision, num_confidence, target, overlap, force_lowest_estimate,
          cost_limit) -> Optional[TaperedQPESizing]:
    u = per_run_qpe_failure(window_cls, num_precision, num_confidence)
    hit = overlap * (1 - u)  # smallest probability that a run is accepted
    max_runs = multi_run_qpe_sizing.max_runs_within(cost_limit, num_precision + num_confidence)
    found = fewest_runs_for_failure(u, hit, target, force_lowest_estimate, max_runs)
    if found is None:
        return None
    runs, quantile = found
    return TaperedQPESizing(num_precision, num_confidence, runs, quantile, 1 - hit,
                            worst_failure(runs, quantile, u, hit), overlap,
                            force_lowest_estimate)
