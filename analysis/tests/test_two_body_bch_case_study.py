from __future__ import annotations

import math

import numpy as np

try:
    from qhat.analysis.analyze_two_body_bch_case_study import (
        BCHComponents,
        exact_paired_sign_permutation_p_value,
        schedule_component_metrics,
    )
    from qhat.analysis.benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
    )
except ImportError:
    from analysis.analyze_two_body_bch_case_study import (
        BCHComponents,
        exact_paired_sign_permutation_p_value,
        schedule_component_metrics,
    )
    from analysis.benchmark_b2_coloring_robustness import (
        HFCommutatorEvaluator,
    )


def test_component_metrics_reconstruct_full_bch_vector() -> None:
    evaluator = HFCommutatorEvaluator(
        number_of_pauli_terms=3,
        left_indices=np.asarray([0, 0, 1], dtype=np.int32),
        right_indices=np.asarray([1, 2, 2], dtype=np.int32),
        target_bins=np.asarray([0, 0, 1], dtype=np.int32),
        pair_amplitudes=np.asarray([1.0, 2.0j, -0.5], dtype=np.complex128),
        number_of_target_bins=2,
    )
    one_one = np.asarray([0.4, 0.5j, -0.2], dtype=np.complex128)
    one_two = np.asarray([0.3, 1.0j, -0.1], dtype=np.complex128)
    two_two = evaluator.pair_amplitudes - one_one - one_two
    components = BCHComponents(one_one, one_two, two_two)

    metrics = schedule_component_metrics(evaluator, components, [0, 1, 2])

    assert math.isclose(metrics["full_bch_norm"], evaluator.evaluate([0, 1, 2]))
    assert 0.0 <= metrics["full_cancellation_ratio"] <= 1.0
    assert 0.0 <= metrics["two_body_cancellation_ratio"] <= 1.0


def test_reversing_order_preserves_component_pair_mass() -> None:
    evaluator = HFCommutatorEvaluator(
        number_of_pauli_terms=3,
        left_indices=np.asarray([0, 0], dtype=np.int32),
        right_indices=np.asarray([1, 2], dtype=np.int32),
        target_bins=np.asarray([0, 0], dtype=np.int32),
        pair_amplitudes=np.asarray([1.0, 1.0], dtype=np.complex128),
        number_of_target_bins=1,
    )
    components = BCHComponents(
        one_one=np.asarray([0.2, 0.2], dtype=np.complex128),
        one_two=np.asarray([0.5, 0.1], dtype=np.complex128),
        two_two=np.asarray([0.3, 0.7], dtype=np.complex128),
    )

    forward = schedule_component_metrics(evaluator, components, [0, 1, 2])
    reverse = schedule_component_metrics(evaluator, components, [2, 1, 0])

    assert math.isclose(
        forward["two_body_cancellation_ratio"],
        reverse["two_body_cancellation_ratio"],
    )
    assert math.isclose(forward["full_bch_norm"], reverse["full_bch_norm"])


def test_exact_paired_sign_test() -> None:
    assert exact_paired_sign_permutation_p_value([1.0, 1.0, 1.0]) == 0.25
    assert math.isnan(exact_paired_sign_permutation_p_value([]))
