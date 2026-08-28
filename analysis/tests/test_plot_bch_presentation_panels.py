from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest
from PIL import Image

try:
    from qhat.analysis.plot_bch_presentation_panels import (
        EXPECTED_HELDOUT_DIRECTION,
        EXPECTED_MATCHED_DIRECTION,
        PRIMARY_BASELINE,
        PROPOSED_ORDERING,
        RANDOM_PARENT_BLOCKS,
        WITHIN_PARENT_SHUFFLE,
        generate_all,
        prepare_heldout_cases,
        prepare_matched_pairs,
        prepare_parent_ablation,
        validate_expected_statistics,
    )
except ImportError:
    from analysis.plot_bch_presentation_panels import (
        EXPECTED_HELDOUT_DIRECTION,
        EXPECTED_MATCHED_DIRECTION,
        PRIMARY_BASELINE,
        PROPOSED_ORDERING,
        RANDOM_PARENT_BLOCKS,
        WITHIN_PARENT_SHUFFLE,
        generate_all,
        prepare_heldout_cases,
        prepare_matched_pairs,
        prepare_parent_ablation,
        validate_expected_statistics,
    )


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CASE_SUMMARY = REPOSITORY_ROOT / (
    "analysis/cancellation_hypothesis_validation/full_analysis/case_summary.csv"
)
MATCHED_PAIRS = REPOSITORY_ROOT / (
    "analysis/cancellation_hypothesis_validation/full_analysis/"
    "matched_pair_summary.csv"
)
PARENT_ABLATION = REPOSITORY_ROOT / (
    "analysis/cancellation_hypothesis_validation/pilot_analysis/"
    "primary_aggregated_deltas.csv"
)


def load_prepared_data():
    heldout, heldout_stats = prepare_heldout_cases(pd.read_csv(CASE_SUMMARY))
    matched, matched_stats = prepare_matched_pairs(
        pd.read_csv(MATCHED_PAIRS),
        heldout,
    )
    parent, medians = prepare_parent_ablation(pd.read_csv(PARENT_ABLATION))
    return heldout, heldout_stats, matched, matched_stats, parent, medians


def test_presentation_statistics_are_recomputed_from_primary_comparison() -> None:
    heldout, heldout_stats, matched, matched_stats, _, medians = load_prepared_data()

    assert set(heldout["baseline_schedule"]) == {PRIMARY_BASELINE}
    assert set(heldout["proposed_schedule"]) == {PROPOSED_ORDERING}
    assert set(matched["baseline_schedule"]) == {PRIMARY_BASELINE}
    assert set(matched["proposed_schedule"]) == {PROPOSED_ORDERING}
    assert (
        heldout_stats.direction_count,
        heldout_stats.direction_total,
    ) == EXPECTED_HELDOUT_DIRECTION
    assert (
        matched_stats.direction_count,
        matched_stats.direction_total,
    ) == EXPECTED_MATCHED_DIRECTION
    validate_expected_statistics(heldout_stats, matched_stats, medians)


def test_parent_ablation_uses_only_requested_schedules_and_collapses_steps() -> None:
    _, _, _, _, parent, medians = load_prepared_data()

    assert set(parent["schedule"]) == {
        RANDOM_PARENT_BLOCKS,
        WITHIN_PARENT_SHUFFLE,
    }
    assert set(parent["reference_schedule"]) == {PROPOSED_ORDERING}
    assert set(parent["steps_collapsed"]) == {"50/100/200"}
    assert len(parent) == 8
    assert parent.groupby("schedule")["case_id"].nunique().to_dict() == {
        RANDOM_PARENT_BLOCKS: 4,
        WITHIN_PARENT_SHUFFLE: 4,
    }
    assert math.isclose(medians[RANDOM_PARENT_BLOCKS][0], 5.947491667, rel_tol=1e-9)
    assert math.isclose(medians[RANDOM_PARENT_BLOCKS][1], 99.118072202, rel_tol=1e-9)
    assert math.isclose(medians[WITHIN_PARENT_SHUFFLE][0], 1.0, abs_tol=1e-12)
    assert math.isclose(
        medians[WITHIN_PARENT_SHUFFLE][1],
        0.999999730,
        abs_tol=1e-9,
    )


def test_wrong_comparator_and_missing_schedule_are_rejected() -> None:
    cases = pd.read_csv(CASE_SUMMARY)
    cases.loc[0, "fresh_jw_magnitude_to_signed_advantage"] *= 2.0
    with pytest.raises(AssertionError, match="Primary error comparator"):
        prepare_heldout_cases(cases)

    ablation = pd.read_csv(PARENT_ABLATION)
    ablation = ablation[ablation["schedule"] != WITHIN_PARENT_SHUFFLE]
    with pytest.raises(AssertionError, match="missing a required schedule"):
        prepare_parent_ablation(ablation)


def test_figure_generation_writes_standalone_vector_and_raster_outputs(
    tmp_path: Path,
) -> None:
    results = generate_all(
        case_summary_path=CASE_SUMMARY,
        matched_pair_path=MATCHED_PAIRS,
        parent_ablation_path=PARENT_ABLATION,
        output_dir=tmp_path,
    )

    expected = {
        "slide2_heldout_cases.png",
        "slide2_heldout_cases.pdf",
        "slide2_heldout_cases_data.csv",
        "slide2_matched_pairs.png",
        "slide2_matched_pairs.pdf",
        "slide2_matched_pairs_data.csv",
        "slide3_parent_block_ablation.png",
        "slide3_parent_block_ablation.pdf",
        "slide3_parent_block_ablation_data.csv",
        "slide3_parent_order_schematic.svg",
        "slide3_parent_order_schematic_data.csv",
        "recomputed_statistics.csv",
    }
    assert {path.name for path in results.output_files} == expected
    for filename in expected:
        path = tmp_path / filename
        assert path.is_file()
        assert path.stat().st_size > 100

    for filename in (
        "slide2_heldout_cases.png",
        "slide2_matched_pairs.png",
        "slide3_parent_block_ablation.png",
    ):
        with Image.open(tmp_path / filename) as image:
            dpi = image.info.get("dpi")
            assert dpi is not None
            assert math.isclose(dpi[0], 300.0, rel_tol=0.01)
            assert math.isclose(dpi[1], 300.0, rel_tol=0.01)
            assert image.convert("RGB").getpixel((0, 0)) == (255, 255, 255)

    for filename in (
        "slide2_heldout_cases.pdf",
        "slide2_matched_pairs.pdf",
        "slide3_parent_block_ablation.pdf",
    ):
        assert (tmp_path / filename).read_bytes().startswith(b"%PDF")

    svg = (tmp_path / "slide3_parent_order_schematic.svg").read_text(
        encoding="utf-8"
    )
    assert svg.startswith("<svg")
    assert "<text" in svg
    assert "Parent-block intervention" in svg
    assert "P1 ↻" in svg
