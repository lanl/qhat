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
split of the failure probability allowed above.

Example: phase error 0.2/1.7 turns, probability of failure 1e-3, overlap 0.5, DPSS window: 3
precision and 4 confidence qubits; one run fails with probability up to 0.50, and the 2nd-smallest
of 14 runs fails with probability up to 9.2e-4.

The per-run failure u is `failure_bound_factor * dpss_average_failure_bound`, an average over the
offset of the true phase from the QPE grid. The paper's random phase shift (arXiv:2403.18927,
Sec. VI) would make it a worst-case bound; qhat does not apply that shift.
"""

import math
from typing import Callable, Optional, Tuple, Type

import attrs
from scipy.special import betainc, betaincc

from qhat.common.qpe_window_state import precision_bits
from qhat.common.tapered_window_state import TaperedWindowStateBase, dpss_average_failure_bound

# Confidence qubits tried by `size_tapered_qpe`.
MIN_CONFIDENCE_QUBITS = 2
MAX_CONFIDENCE_QUBITS = 12
# Largest cost N * 2^(number of phase qubits) considered.
COST_CEILING = 1e18
# With no answer yet, stop trying larger K after this many in a row fail to get closer.
_QUANTILE_PATIENCE = 24
# Larger K tried after the first one that works, in case one of them needs fewer runs.
_QUANTILE_GUARD = 2


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
    return _worst_failure(num_runs, quantile, u, overlap * (1 - u))


def fewest_runs(window_cls: Type[TaperedWindowStateBase], num_precision_qubits: int,
                num_confidence_qubits: int, probability_of_failure: float, overlap: float,
                force_lowest_estimate: bool = False) -> Optional[TaperedQPESizing]:
    """The fewest runs N (and the best K for that N; K = 1 if `force_lowest_estimate`) that meet
    `probability_of_failure` with these qubit counts, or None if none does within
    `COST_CEILING`.
    """
    _check_targets(probability_of_failure, overlap)
    return _plan(window_cls, num_precision_qubits, num_confidence_qubits, probability_of_failure,
                 overlap, force_lowest_estimate, COST_CEILING)


def size_tapered_qpe(window_cls: Type[TaperedWindowStateBase], phase_error: float,
                     probability_of_failure: float, overlap: float,
                     force_lowest_estimate: bool = False) -> Optional[TaperedQPESizing]:
    """The cheapest plan meeting `phase_error` (turns) and `probability_of_failure` (for the
    reported estimate of N runs), or None if none costs at most `COST_CEILING`.

    Precision qubits: the fewest with 2^-(num_precision_qubits + 1) <= phase_error. Confidence
    qubits (MIN_CONFIDENCE_QUBITS to MAX_CONFIDENCE_QUBITS), N and K minimize the cost
    N * 2^(number of phase qubits); ties go to fewer qubits. K = 1 if `force_lowest_estimate`.
    """
    _check_targets(probability_of_failure, overlap)
    num_precision = precision_bits(phase_error) - 1
    best, limit = None, COST_CEILING
    for num_confidence in range(MIN_CONFIDENCE_QUBITS, MAX_CONFIDENCE_QUBITS + 1):
        if 2.0**(num_precision + num_confidence) > limit:
            break
        plan = _plan(window_cls, num_precision, num_confidence, probability_of_failure, overlap,
                     force_lowest_estimate, limit)
        if plan is None:
            continue
        cost = plan.num_runs * 2.0**plan.num_phase_qubits
        if best is None or cost < limit:
            best, limit = plan, cost
    return best


def _check_targets(probability_of_failure: float, overlap: float):
    if not 0 < probability_of_failure < 1:
        raise ValueError(
            f"probability of failure must be in (0, 1), got {probability_of_failure}.")
    if not 0 < overlap <= 1:
        raise ValueError(f"overlap must be in (0, 1], got {overlap}.")


def _plan(window_cls, num_precision, num_confidence, target, overlap, force_lowest_estimate,
          cost_limit) -> Optional[TaperedQPESizing]:
    u = per_run_qpe_failure(window_cls, num_precision, num_confidence)
    hit = overlap * (1 - u)  # smallest probability that a run is accepted
    if hit <= 0:
        return None
    cap = int(cost_limit // 2.0**(num_precision + num_confidence))
    # Even K = 1 fails if every run overestimates, which has probability up to (1 - hit)^N.
    min_runs = math.log(target) / math.log1p(-hit)
    if min_runs > cap:
        return None
    found = _fewest_runs(u, hit, target, force_lowest_estimate,
                         max(1, math.ceil(min_runs - 1e-9)), cap)
    if found is None:
        return None
    runs, quantile = found
    return TaperedQPESizing(num_precision, num_confidence, runs, quantile, 1 - hit,
                            _worst_failure(runs, quantile, u, hit), overlap,
                            force_lowest_estimate)


def _failure(n: int, k: int, under: float, not_over: float) -> float:
    """P(the k-th smallest of n estimates is an under- or overestimate) when each run independently
    underestimates with probability `under` and overestimates with probability 1 - `not_over`.

    Underestimate: at least k runs underestimate. Overestimate: at most k - 1 runs do not
    overestimate. Both terms are computed directly, so a small failure keeps its precision.
    """
    low = betainc(k, n - k + 1, under)
    high = betaincc(k, n - k + 1, not_over)
    return min(1.0, float(low + high))


def _worst_failure(n: int, k: int, u: float, hit: float) -> float:
    """`_failure` at the worse of the two extreme splits of the single-run failure 1 - hit: all
    overestimates, or underestimates u and the rest overestimates. The worst split is one of the
    two."""
    return max(_failure(n, k, 0.0, hit), _failure(n, k, u, min(1.0, u + hit)))


def _first(pred: Callable[[int], bool], lo: int, hi: int) -> int:
    """Smallest x in [lo, hi] with pred(x), for pred false then true, and pred(hi) true."""
    while lo < hi:
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def _last(pred: Callable[[int], bool], lo: int, hi: int) -> int:
    """Largest x in [lo, hi] with pred(x), for pred true then false, and pred(lo) true."""
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if pred(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def _argmin(f: Callable[[int], float], lo: int, hi: int) -> int:
    """x in [lo, hi] minimizing f, for f falling then rising."""
    while hi - lo > 2:
        m1, m2 = lo + (hi - lo) // 3, hi - (hi - lo) // 3
        if f(m1) <= f(m2):
            hi = m2
        else:
            lo = m1
    return min(range(lo, hi + 1), key=f)


def _runs_for_quantile(k: int, u: float, hit: float, target: float, n_min: int,
                       cap: int) -> Tuple[Optional[int], float]:
    """Fewest runs n in [max(k, n_min), cap] with `_worst_failure(n, k, ...) <= target`, or None;
    also the smallest worst failure seen, to tell when larger k stop helping.

    The all-overestimate split only improves with more runs, which allows n >= n_over. In the other
    split, more runs mean more chances to underestimate: its underestimate term alone allows
    n <= n_under, and between n_over and n_under its failure falls, then rises.
    """
    lo = max(k, n_min)
    if lo > cap:
        return None, 1.0
    if _failure(cap, k, 0.0, hit) > target:
        return None, _failure(cap, k, 0.0, hit)
    n_over = _first(lambda n: _failure(n, k, 0.0, hit) <= target, lo, cap)
    split = lambda n: _failure(n, k, u, min(1.0, u + hit))
    if _failure(n_over, k, u, 1.0) > target:
        return None, _worst_failure(n_over, k, u, hit)
    n_under = _last(lambda n: _failure(n, k, u, 1.0) <= target, n_over, cap)
    best = _argmin(split, n_over, n_under)
    if split(best) > target:
        return None, split(best)
    return _first(lambda n: split(n) <= target, n_over, best), split(best)


def _fewest_runs(u: float, hit: float, target: float, force_lowest_estimate: bool, n_min: int,
                 cap: int) -> Optional[Tuple[int, int]]:
    """(N, K) with the fewest runs N in [n_min, cap] meeting `target`, and the K with the smallest
    failure at that N (K = 1 if `force_lowest_estimate`); None if no N works."""
    if force_lowest_estimate:
        runs, _ = _runs_for_quantile(1, u, hit, target, n_min, cap)
        return None if runs is None else (runs, 1)

    best = None
    nearest, stale, extra = math.inf, 0, 0
    k = 1
    # Runs for quantile k number at least k, so once k reaches the best N nothing can beat it.
    while k <= cap and (best is None or k < best[0]):
        # The all-overestimate failure rises with k, so if it misses at cap, so do larger k.
        if _failure(cap, k, 0.0, hit) > target:
            break
        runs, seen = _runs_for_quantile(k, u, hit, target, n_min, cap)
        if runs is not None and (best is None or runs < best[0]):
            best = (runs, k)
        if best is None:
            nearest, stale = (seen, 0) if seen < nearest else (nearest, stale + 1)
            if stale >= _QUANTILE_PATIENCE:
                break
        else:
            extra += 1
            if extra > _QUANTILE_GUARD:
                break
        k += 1
    if best is None:
        return None
    runs, k = best
    return runs, _best_quantile(runs, k, u, hit, target)


def _best_quantile(n: int, k_found: int, u: float, hit: float, target: float) -> int:
    """K with the smallest worst failure at n runs. Searched where each split's dominant term
    meets the target (k_found does), so the failure falls, then rises."""
    k_lo = _first(lambda k: _failure(n, k, u, 1.0) <= target, 1, k_found)
    k_hi = _last(lambda k: _failure(n, k, 0.0, hit) <= target, k_found, n)
    worst = lambda k: _worst_failure(n, k, u, hit)
    k_best = _argmin(worst, k_lo, k_hi)
    return k_best if worst(k_best) < worst(k_found) else k_found
