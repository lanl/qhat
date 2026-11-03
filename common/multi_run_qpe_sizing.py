"""
Window-independent part of sizing QPE for N runs when the initial state only partly overlaps the
target eigenstate: the failure of reporting the K-th smallest of N estimates, and the fewest N (and
best K) meeting a target. This implements the multi-run analysis by Yigit Subasi and Alan Bu
(unpublished). The window-specific sizing is in `qhat.common.tapered_qpe_sizing` and
`qhat.common.textbook_qpe_sizing`.

One run returns one phase estimate. With per-run QPE failure u and overlap p between the initial
state and the target eigenstate (the lowest phase), the estimate is

- accepted (within the phase error of the target phase) with probability at least
  hit = p (1 - u);
- an underestimate with probability at most u. Nothing lies below the target, so only a failed
  QPE gives an underestimate;
- otherwise an overestimate.

So one run fails with probability at most 1 - hit. Reporting the K-th smallest of N runs fails if
that estimate is an under- or overestimate; `worst_failure` bounds this over every split of the
failure probability allowed above.

Example: u = 2.7e-4 and overlap 0.5 (hit = 0.49986): the 2nd-smallest of 14 runs fails with
probability up to 9.2e-4.

This assumes the phase window leaves at least the phase error of margin above the highest
eigenphase; otherwise a correct estimate of a high eigenphase can wrap past phase 1 and decode as
an underestimate.
"""

import math
from typing import Callable, Optional, Tuple

from scipy.special import betainc, betaincc

# Largest cost N * 2^(number of phase qubits) considered.
COST_CEILING = 1e18
# With no answer yet, stop trying larger K after this many in a row fail to get closer.
_QUANTILE_PATIENCE = 24
# Larger K tried after the first one that works, in case one of them needs fewer runs.
_QUANTILE_GUARD = 2
# The K loop starts this far below the first K that looks useful (see
# `_smallest_useful_quantile`): integer rounding makes usefulness jitter near the switch.
_QUANTILE_BACKOFF = 64


def check_targets(probability_of_failure: float, overlap: float):
    """Raise a ValueError unless 0 < `probability_of_failure` < 1 and 0 < `overlap` <= 1."""
    if not 0 < probability_of_failure < 1:
        raise ValueError(
            f"probability of failure must be in (0, 1), got {probability_of_failure}.")
    if not 0 < overlap <= 1:
        raise ValueError(f"overlap must be in (0, 1], got {overlap}.")


def worst_failure(num_runs: int, quantile: int, u: float, hit: float) -> float:
    """Bound on the probability that the `quantile`-th smallest of `num_runs` estimates is an
    under- or overestimate, when each run underestimates with probability at most `u` and is
    accepted with probability at least `hit`. E.g. 14 runs, quantile 2, u = 2.7e-4,
    hit = 0.49986: 9.2e-4.

    The worst split of the single-run failure 1 - hit is one of two extremes: all overestimates,
    or underestimates u and the rest overestimates.
    """
    return max(_failure(num_runs, quantile, 0.0, hit),
               _failure(num_runs, quantile, u, min(1.0, u + hit)))


def fewest_runs_for_failure(u: float, hit: float, probability_of_failure: float,
                            force_lowest_estimate: bool, max_runs: int) -> Optional[Tuple[int, int]]:
    """(N, K): the fewest runs N <= `max_runs` with `worst_failure(N, K, u, hit)` at most
    `probability_of_failure`, and the K with the smallest failure at that N (K = 1 if
    `force_lowest_estimate`). None if no N works. E.g. u = 2.7e-4, hit = 0.49986, target 1e-3:
    (14, 2).
    """
    if hit <= 0 or max_runs < 1:
        return None
    # If underestimates are at least as likely as accepted runs, one of the two extreme splits
    # fails with probability at least 1/2, whatever N and K.
    if u >= hit and probability_of_failure < 0.5:
        return None
    # Even K = 1 fails if every run overestimates, which has probability up to (1 - hit)^N.
    min_runs = 1.0 if hit >= 1 else math.log(probability_of_failure) / math.log1p(-hit)
    if min_runs > max_runs:
        return None
    return _fewest_runs(u, hit, probability_of_failure, force_lowest_estimate,
                        max(1, math.ceil(min_runs - 1e-9)), max_runs)


def _failure(n: int, k: int, under: float, not_over: float) -> float:
    """P(the k-th smallest of n estimates is an under- or overestimate) when each run independently
    underestimates with probability `under` and overestimates with probability 1 - `not_over`.

    Underestimate: at least k runs underestimate. Overestimate: at most k - 1 runs do not
    overestimate. Both terms are computed directly, so a small failure keeps its precision.
    """
    low = betainc(k, n - k + 1, under)
    high = betaincc(k, n - k + 1, not_over)
    return min(1.0, float(low + high))


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
    """Fewest runs n in [max(k, n_min), cap] with `worst_failure(n, k, ...) <= target`, or None;
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
        return None, worst_failure(n_over, k, u, hit)
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
    k_useful = _smallest_useful_quantile(u, hit, target, n_min, cap)
    k = max(1, k_useful - _QUANTILE_BACKOFF)
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
            if stale >= _QUANTILE_PATIENCE and k >= k_useful:
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


def _smallest_useful_quantile(u: float, hit: float, target: float, n_min: int, cap: int) -> int:
    """Roughly the smallest K for which some N can work, so the K loop need not start at 1.

    K is useless if, at the fewest runs n_over(K) that its all-overestimate split allows, its
    underestimate term alone already misses the target: more runs only add underestimates. As K
    grows this switches from useless to useful, apart from jitter of a few K from integer
    rounding of n_over(K), so a galloping search finds the switch approximately; the K loop
    starts `_QUANTILE_BACKOFF` below it. E.g. u = 9.8e-4, hit = 1e-3, target 1e-9: K up to
    about 4e5 are useless.
    """
    def useful(k):
        lo = max(k, n_min)
        if lo > cap or _failure(cap, k, 0.0, hit) > target:
            return True  # let the K loop's own stopping rules handle it
        n_over = _first(lambda n: _failure(n, k, 0.0, hit) <= target, lo, cap)
        return _failure(n_over, k, u, 1.0) <= target

    if useful(1):
        return 1
    k_useless, k = 1, 2
    while not useful(k):
        k_useless, k = k, 2 * k
    return _first(useful, k_useless + 1, k)


def _best_quantile(n: int, k_found: int, u: float, hit: float, target: float) -> int:
    """K with the smallest worst failure at n runs. Searched where each split's dominant term
    meets the target (k_found does), so the failure falls, then rises."""
    k_lo = _first(lambda k: _failure(n, k, u, 1.0) <= target, 1, k_found)
    k_hi = _last(lambda k: _failure(n, k, 0.0, hit) <= target, k_found, n)
    worst = lambda k: worst_failure(n, k, u, hit)
    k_best = _argmin(worst, k_lo, k_hi)
    return k_best if worst(k_best) < worst(k_found) else k_found
