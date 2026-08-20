from __future__ import annotations

import math

import numpy as np
import pandas as pd

try:
    from qhat.analysis.analyze_signed_order_cancellation_structure import (
        load_cancellation_measurements,
        order_interaction_metrics,
        random_order_local_pair_probability,
        random_order_mean_span_fraction,
    )
except ImportError:
    from analysis.analyze_signed_order_cancellation_structure import (
        load_cancellation_measurements,
        order_interaction_metrics,
        random_order_local_pair_probability,
        random_order_mean_span_fraction,
    )


def test_random_order_pair_baselines() -> None:
    # Four positions have six unordered pairs; three are adjacent.
    assert math.isclose(random_order_local_pair_probability(4, 1), 0.5)
    # E[distance] is (N+1)/3, normalized here by N-1.
    assert math.isclose(random_order_mean_span_fraction(4), 5.0 / 9.0)


def test_order_metrics_use_coefficient_product_and_shared_support_weight() -> None:
    coefficients = np.array([2.0, -3.0, 4.0])
    shared = np.array(
        [
            [0, 2, 1],
            [2, 0, 1],
            [1, 1, 0],
        ],
        dtype=np.uint8,
    )
    metrics = order_interaction_metrics(
        coefficients,
        shared,
        order=[0, 1, 2],
        prefix="test",
    )

    # Weights: w01=12 (opposite), w02=8 (same), w12=12 (opposite).
    assert math.isclose(metrics["test_total_interaction_weight"], 32.0)
    assert math.isclose(
        metrics["test_opposite_sign_interaction_weight_fraction"],
        0.75,
    )
    # Weighted separations are 12*1 + 8*2 + 12*1; normalize by 32*(N-1).
    assert math.isclose(metrics["test_weighted_mean_span_fraction"], 40.0 / 64.0)
    # For N=3, the 10% cutoff rounds up to one position.
    assert math.isclose(metrics["test_weighted_local_10pct_fraction"], 0.75)
    assert math.isclose(
        metrics["test_opposite_sign_local_10pct_interaction_fraction"],
        0.75,
    )
    assert math.isclose(
        metrics["test_opposite_sign_local_10pct_enrichment"],
        4.0 / 3.0,
    )


def test_order_metrics_change_locality_but_not_global_sign_share() -> None:
    coefficients = np.array([-4.0, -1.0, 1.0, 4.0])
    shared = np.zeros((4, 4), dtype=np.uint8)
    shared[0, 2] = shared[2, 0] = 1
    shared[1, 3] = shared[3, 1] = 1

    signed = order_interaction_metrics(
        coefficients,
        shared,
        order=[0, 1, 2, 3],
        prefix="signed",
    )
    permuted = order_interaction_metrics(
        coefficients,
        shared,
        order=[0, 3, 1, 2],
        prefix="permuted",
    )

    assert signed["signed_opposite_sign_interaction_weight_fraction"] == 1.0
    assert permuted["permuted_opposite_sign_interaction_weight_fraction"] == 1.0
    assert signed["signed_weighted_local_10pct_fraction"] == 0.0
    assert permuted["permuted_weighted_local_10pct_fraction"] > 0.0


def test_load_cancellation_measurements_builds_signed_strength(tmp_path) -> None:
    path = tmp_path / "ablation.csv"
    pd.DataFrame(
        [
            {
                "case_id": "case",
                "schedule": "fermionic_signed_reference",
                "status": "success",
                "bch_cancellation_ratio": 0.01,
                "bch2_hf_state_norm": 1.0,
                "bch_pair_abs_sum": 100.0,
            },
            {
                "case_id": "case",
                "schedule": "fermionic_magnitude_reference",
                "status": "success",
                "bch_cancellation_ratio": 0.1,
                "bch2_hf_state_norm": 10.0,
                "bch_pair_abs_sum": 100.0,
            },
            {
                "case_id": "case",
                "schedule": "jw_signed_baseline",
                "status": "success",
                "bch_cancellation_ratio": 0.2,
                "bch2_hf_state_norm": 20.0,
                "bch_pair_abs_sum": 100.0,
            },
        ]
    ).to_csv(path, index=False)

    result = load_cancellation_measurements([path]).iloc[0]
    assert math.isclose(result["signed_bch_cancellation_strength"], 2.0)
    assert math.isclose(
        result["signed_bch_cancellation_advantage_vs_fermionic_magnitude"],
        10.0,
    )
    assert math.isclose(
        result["signed_bch_cancellation_advantage_vs_jw_signed"],
        20.0,
    )
