from __future__ import annotations

import math

import pandas as pd

try:
    from qhat.analysis.analyze_cancellation_hypothesis_validation import (
        REFERENCE,
        aggregate_random_replicates,
        build_case_summary,
        build_matched_pair_summary,
        build_within_case_deltas,
        collapse_conditions_for_inference,
    )
    from qhat.analysis.run_cancellation_hypothesis_validation import (
        DISCOVERY_CASE_IDS,
        selected_panel,
    )
except ImportError:
    from analysis.analyze_cancellation_hypothesis_validation import (
        REFERENCE,
        aggregate_random_replicates,
        build_case_summary,
        build_matched_pair_summary,
        build_within_case_deltas,
        collapse_conditions_for_inference,
    )
    from analysis.run_cancellation_hypothesis_validation import (
        DISCOVERY_CASE_IDS,
        selected_panel,
    )


def manifest() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "case_id": "case-a",
                "molecule": "A",
                "matched_pair": "pair",
                "expected_outcome": "favorable",
                "fermionic_advantage_factor": 2.0,
                "n_qubits": 12,
                "historical_fermionic_one_minus_overlap": 0.01,
            }
        ]
    )


def result_row(
    schedule: str,
    cancellation_ratio: float,
    error: float,
    sample_index: int = -1,
) -> dict[str, object]:
    return {
        "case_id": "case-a",
        "molecule": "A",
        "n_qubits": 12,
        "schedule": schedule,
        "sample_index": sample_index,
        "trotter_steps": 100,
        "evolution_time": 1.0,
        "bch_cancellation_ratio": cancellation_ratio,
        "one_minus_overlap": error,
    }


def test_pilot_is_held_out_and_active_space_matched() -> None:
    panel = selected_panel("pilot")
    assert len(panel) == 4
    assert DISCOVERY_CASE_IDS.isdisjoint(case.case_id for case in panel)

    by_pair: dict[str, list[object]] = {}
    for case in panel:
        by_pair.setdefault(case.matched_pair, []).append(case)
    for cases in by_pair.values():
        assert {case.expected_outcome for case in cases} == {
            "favorable",
            "negative_control",
        }


def test_random_replicates_are_reduced_to_one_log_space_median() -> None:
    results = pd.DataFrame(
        [
            result_row(REFERENCE, 0.1, 0.01),
            result_row("signed_parent_blocks_randomized", 0.2, 0.02, 0),
            result_row("signed_parent_blocks_randomized", 0.8, 0.32, 1),
        ]
    )
    deltas = build_within_case_deltas(results, manifest(), error_floor=1.0e-13)
    aggregated = aggregate_random_replicates(deltas)

    assert len(deltas) == 2
    assert len(aggregated) == 1
    row = aggregated.iloc[0]
    assert row["replicates"] == 2
    assert math.isclose(row["delta_log10_cancellation"], math.log10(4.0))
    assert math.isclose(row["cancellation_ratio_to_signed"], 4.0)
    assert math.isclose(row["delta_log10_error"], math.log10(8.0))
    assert math.isclose(row["error_ratio_to_signed"], 8.0)


def test_case_summary_recomputes_advantage_from_current_run() -> None:
    results = pd.DataFrame(
        [
            result_row(REFERENCE, 0.1, 0.01),
            result_row("jw_signed_baseline", 0.2, 0.03),
            result_row("jw_magnitude_baseline", 0.15, 0.02),
            result_row("fermionic_magnitude_reference", 0.3, 0.04),
        ]
    )
    deltas = build_within_case_deltas(results, manifest(), error_floor=1.0e-13)
    aggregated = aggregate_random_replicates(deltas)
    summary = build_case_summary(aggregated, results, manifest()).iloc[0]

    assert math.isclose(summary["fresh_jw_to_signed_advantage"], 2.0)
    assert summary["fresh_best_jw_schedule"] == "jw_magnitude_baseline"
    assert math.isclose(
        summary["fresh_best_jw_to_signed_bch_cancellation_ratio"],
        1.5,
    )
    assert summary["fresh_observed_outcome"] == "favorable"
    assert bool(summary["performance_label_reproduced"])


def test_repeated_step_counts_do_not_inflate_primary_units() -> None:
    aggregated = pd.DataFrame(
        [
            {
                "case_id": "case-a",
                "schedule": "jw_signed_baseline",
                "trotter_steps": 50,
                "delta_log10_cancellation": math.log10(2.0),
                "delta_log10_error": math.log10(4.0),
            },
            {
                "case_id": "case-a",
                "schedule": "jw_signed_baseline",
                "trotter_steps": 200,
                "delta_log10_cancellation": math.log10(8.0),
                "delta_log10_error": math.log10(16.0),
            },
        ]
    )
    collapsed = collapse_conditions_for_inference(aggregated)

    assert len(collapsed) == 1
    assert math.isclose(collapsed.iloc[0]["cancellation_ratio_to_signed"], 4.0)
    assert math.isclose(collapsed.iloc[0]["error_ratio_to_signed"], 8.0)


def test_matched_pair_summary_uses_favorable_minus_control_deltas() -> None:
    cases = pd.DataFrame(
        [
            {
                "case_id": "fav",
                "matched_pair": "pair",
                "expected_outcome": "favorable",
                "n_qubits": 12,
                "fresh_best_jw_to_signed_bch_cancellation_ratio": 4.0,
                "fresh_jw_to_signed_advantage": 8.0,
            },
            {
                "case_id": "control",
                "matched_pair": "pair",
                "expected_outcome": "negative_control",
                "n_qubits": 12,
                "fresh_best_jw_to_signed_bch_cancellation_ratio": 2.0,
                "fresh_jw_to_signed_advantage": 2.0,
            },
        ]
    )
    pair = build_matched_pair_summary(cases).iloc[0]

    assert math.isclose(pair["delta_log10_relative_cancellation"], math.log10(2.0))
    assert math.isclose(pair["delta_log10_fresh_advantage"], math.log10(4.0))
    assert bool(pair["direction_concordant"])
