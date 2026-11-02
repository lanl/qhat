"""
Sizing textbook QPE (the rectangular window, `qhat.common.qpe_window_state.RectangularWindowState`)
when the initial state only partly overlaps the target eigenstate: how many phase qubits m to use,
how many times N to run the circuit, and which estimate to report (the K-th smallest of the N).
This applies the multi-run analysis by Yigit Subasi and Alan Bu (unpublished; see
`qhat.common.multi_run_qpe_sizing`) with a per-run bound for the rectangular window.

Per-run bound: with phase error W (turns) and h = W 2^m, one run of textbook QPE misses the true
phase by more than W with probability less than (2/pi^2) psi'(h), psi' being the trigamma
function. This holds for every true phase, so no random phase shift is needed. E.g. W = 0.2/1.7
and m = 5 (h = 3.76): 6.2e-2. Nielsen & Chuang Eq. 5.35 is 2.5-6.6x larger.

Why: the outcome d grid steps from the phase has probability
sin^2(pi f) / pi^2 * sum_k 1/(d - k 2^m - f)^2 (f = offset of the phase from the grid). A
rejected outcome has every term more than h from f, and each tail of 1/(n - f)^2 over integers n
sums to less than psi'(h).

Example: W = 0.2/1.7 turns, probability of failure 1e-3, overlap 0.5: 5 phase qubits; one run fails
with probability up to 0.53, and the 9th-smallest of 38 runs fails with probability up to 8.9e-4.
"""

import math
from typing import Optional

import attrs
from scipy.special import polygamma

from qhat.common.multi_run_qpe_sizing import (
    check_targets,
    COST_CEILING,
    fewest_runs_for_failure,
    worst_failure,
)

# Relative margin on the computed bound, covering floating-point error in polygamma.
_ROUNDOFF_MARGIN = 1e-12


@attrs.frozen
class TextbookQPESizing:
    """Phase qubits and run plan. E.g. 5 phase qubits, run 38 times and report the 9th-smallest
    estimate (`quantile` 9): `single_run_failure` 0.53, `multi_run_failure` 8.9e-4 at
    `overlap` 0.5.
    """

    num_phase_qubits: int
    num_runs: int
    quantile: int
    single_run_failure: float
    multi_run_failure: float
    overlap: float
    force_lowest_estimate: bool


def per_run_qpe_failure(phase_error: float, num_phase_qubits: int) -> float:
    """u = (2/pi^2) psi'(h), h = phase_error * 2^num_phase_qubits, at most 1: a bound on the
    probability that one run misses the true phase by more than `phase_error` turns, for every
    true phase. E.g. 0.2/1.7 turns with 5 phase qubits: 6.2e-2."""
    _check_phase_error(phase_error)
    h = phase_error * 2.0**num_phase_qubits
    return min(1.0, 2 / math.pi**2 * float(polygamma(1, h)) * (1 + _ROUNDOFF_MARGIN))


def multi_run_failure(phase_error: float, num_phase_qubits: int, num_runs: int, quantile: int,
                      overlap: float) -> float:
    """Bound on the probability that the `quantile`-th smallest of `num_runs` estimates misses the
    target phase by more than `phase_error` turns. E.g. 0.2/1.7 turns, 5 phase qubits, 38 runs,
    quantile 9: 8.9e-4 at overlap 0.5.
    """
    u = per_run_qpe_failure(phase_error, num_phase_qubits)
    return worst_failure(num_runs, quantile, u, overlap * (1 - u))


def fewest_runs(phase_error: float, num_phase_qubits: int, probability_of_failure: float,
                overlap: float, force_lowest_estimate: bool = False) -> Optional[TextbookQPESizing]:
    """The fewest runs N (and the best K for that N; K = 1 if `force_lowest_estimate`) that meet
    `probability_of_failure` with this many phase qubits, or None if none does within
    `COST_CEILING`.
    """
    check_targets(probability_of_failure, overlap)
    return _plan(phase_error, num_phase_qubits, probability_of_failure, overlap,
                 force_lowest_estimate, COST_CEILING)


def size_textbook_qpe(phase_error: float, probability_of_failure: float, overlap: float,
                      force_lowest_estimate: bool = False) -> Optional[TextbookQPESizing]:
    """The cheapest plan meeting `phase_error` (turns) and `probability_of_failure` (for the
    reported estimate of N runs), or None if none costs at most `COST_CEILING`.

    The phase qubits m, N and K minimize the cost N * 2^m; ties go to fewer qubits. K = 1 if
    `force_lowest_estimate`.
    """
    _check_phase_error(phase_error)
    check_targets(probability_of_failure, overlap)
    best, limit = None, COST_CEILING
    # Most qubits first: those plans need few runs and are quick to find, so the cost limit is
    # tight by the time the search reaches qubit counts that would need very many runs.
    for num_phase_qubits in range(int(math.log2(COST_CEILING)), 0, -1):
        if 2.0**num_phase_qubits > limit:
            continue
        u = per_run_qpe_failure(phase_error, num_phase_qubits)
        if u >= overlap * (1 - u) and probability_of_failure < 0.5:
            break  # no plan here (see `fewest_runs_for_failure`), and fewer qubits only raise u
        plan = _plan(phase_error, num_phase_qubits, probability_of_failure, overlap,
                     force_lowest_estimate, limit)
        if plan is None:
            continue
        cost = plan.num_runs * 2.0**plan.num_phase_qubits
        if best is None or cost <= limit:
            best, limit = plan, cost
    return best


def _check_phase_error(phase_error: float):
    if not 0 < phase_error < 0.5:
        raise ValueError(f"phase error must be in (0, 0.5) turns, got {phase_error}.")


def _plan(phase_error, num_phase_qubits, target, overlap, force_lowest_estimate,
          cost_limit) -> Optional[TextbookQPESizing]:
    u = per_run_qpe_failure(phase_error, num_phase_qubits)
    hit = overlap * (1 - u)  # smallest probability that a run is accepted
    max_runs = int(cost_limit // 2.0**num_phase_qubits)
    found = fewest_runs_for_failure(u, hit, target, force_lowest_estimate, max_runs)
    if found is None:
        return None
    runs, quantile = found
    return TextbookQPESizing(num_phase_qubits, runs, quantile, 1 - hit,
                             worst_failure(runs, quantile, u, hit), overlap, force_lowest_estimate)
