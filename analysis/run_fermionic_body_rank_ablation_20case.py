#!/usr/bin/env python3
"""
Run the selective 1B-vs-2B fermionic-parent ordering ablation on the exact
20-case BCH validation panel.

COPY THIS ONE FILE TO:
    qhat/analysis/run_fermionic_body_rank_ablation_20case.py

It depends only on the EXISTING repository file:
    qhat/analysis/benchmark_fermionic_structure_ablation.py

It does NOT require benchmark_fermionic_body_rank_ablation.py.

Experiment
----------
For every case in:
    analysis/cancellation_hypothesis_validation/full_panel_manifest.csv

compare:
    1. fermionic_signed_reference
    2. signed_parent_1b_blocks_randomized
    3. signed_parent_2b_blocks_randomized

The selective shuffle:
    * keeps each parent-owned Pauli block intact;
    * keeps the 1B/2B slot pattern fixed;
    * shuffles only parent identities within the selected body rank;
    * leaves the other body rank fixed;
    * reuses the repository's existing BCH and Trotter-error machinery.

Typical use from qhat/:
    python analysis/run_fermionic_body_rank_ablation_20case.py --dry-run
    python analysis/run_fermionic_body_rank_ablation_20case.py
"""

from __future__ import annotations

import argparse
import copy
import csv
import sys
import zlib
from pathlib import Path
from typing import Any, Sequence

import numpy as np

# Because this file is placed in qhat/analysis/, Python can import the existing
# analysis script directly from the same directory.
import benchmark_fermionic_structure_ablation as base


REFERENCE = base.REFERENCE
SHUFFLE_1B = "signed_parent_1b_blocks_randomized"
SHUFFLE_2B = "signed_parent_2b_blocks_randomized"

SELECTIVE_RANDOM_SCHEDULES = (
    SHUFFLE_1B,
    SHUFFLE_2B,
)

_ORIGINAL_BUILD_SCHEDULE_SPECS = base.build_schedule_specs
_ORIGINAL_ORDERING_DEFINITION = base.ordering_definition


def parse_runner_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "analysis/cancellation_hypothesis_validation/full_panel_manifest.csv"
        ),
        help="Frozen 20-case BCH validation manifest.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/fermionic_body_rank_ablation_20case.csv"),
        help="Append/resume result CSV.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        default=[100],
        help="First-order Trotter step counts. Default: 100.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        dest="evolution_time",
        help="Total evolution time. Default: 1.0.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=20,
        help="Random samples for each of the 1B-only and 2B-only shuffles.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260807,
        help="Base random seed.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=base.DEFAULT_TOLERANCE,
        help="Coefficient/commutator tolerance.",
    )
    parser.add_argument(
        "--parallel-threshold",
        type=int,
        default=2**16,
        help="Parallel Pauli-kernel threshold.",
    )
    parser.add_argument(
        "--no-spin-sector",
        action="store_true",
        help="Pass through to the existing ablation benchmark.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve and print the 20 tensor paths, then stop.",
    )
    return parser.parse_args()


def repository_root() -> Path:
    # qhat/analysis/this_file.py -> qhat/
    return Path(__file__).resolve().parents[1]


def resolve_path_from_root(root: Path, path: Path) -> Path:
    if path.is_absolute():
        return path.resolve()
    return (root / path).resolve()


def resolve_tensor(
    *,
    root: Path,
    case_id: str,
    manifest_tensor_path: str | None,
) -> Path:
    """
    Resolve the exact tensor for a manifest row.

    Prefer manifest tensor_path. If that path is absent/stale, search by the
    exact case filename under hamiltonian_generator and require a unique match.
    """
    if manifest_tensor_path:
        candidate = Path(manifest_tensor_path)
        if not candidate.is_absolute():
            candidate = root / candidate
        if candidate.is_file():
            return candidate.resolve()

    expected_name = f"{case_id}.tensors.npz"
    search_root = root / "hamiltonian_generator"
    matches = sorted(search_root.rglob(expected_name))

    if not matches:
        raise FileNotFoundError(
            f"No tensor found for case {case_id!r}.\n"
            f"Expected {expected_name!r} somewhere under:\n"
            f"  {search_root}"
        )

    if len(matches) > 1:
        listed = "\n".join(f"  - {path}" for path in matches)
        raise RuntimeError(
            f"Case {case_id!r} matched more than one tensor:\n"
            f"{listed}\n"
            "Refusing to guess which tensor belongs to the frozen panel."
        )

    return matches[0].resolve()


def load_panel(manifest: Path, root: Path) -> list[tuple[str, Path]]:
    with manifest.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    if not rows:
        raise RuntimeError(f"Manifest is empty: {manifest}")

    if "case_id" not in rows[0]:
        raise RuntimeError(
            f"Manifest has no 'case_id' column: {manifest}"
        )

    if len(rows) != 20:
        raise RuntimeError(
            f"Expected exactly 20 rows in the frozen panel, found {len(rows)} "
            f"in {manifest}."
        )

    panel: list[tuple[str, Path]] = []
    seen: set[str] = set()

    for row_number, row in enumerate(rows, start=2):
        case_id = str(row.get("case_id", "")).strip()
        if not case_id:
            raise RuntimeError(
                f"Manifest row {row_number} has an empty case_id."
            )
        if case_id in seen:
            raise RuntimeError(
                f"Duplicate case_id in manifest: {case_id}"
            )
        seen.add(case_id)

        manifest_tensor_path = str(row.get("tensor_path", "")).strip() or None
        tensor = resolve_tensor(
            root=root,
            case_id=case_id,
            manifest_tensor_path=manifest_tensor_path,
        )
        panel.append((case_id, tensor))

    return panel


def display_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def fermionic_body_rank(term: Any) -> int:
    """
    Classify a complete Hermitian fermionic parent as 1B or 2B.

    The repository's Hermitian parent object stores its normal-ordered
    fermionic monomial keys in term.component_keys.

        1B -> 2 ladder operators
        2B -> 4 ladder operators
    """
    ranks: set[int] = set()

    for component_key in term.component_keys:
        n_operators = len(component_key)

        if n_operators == 0:
            continue
        if n_operators == 2:
            ranks.add(1)
        elif n_operators == 4:
            ranks.add(2)
        else:
            raise ValueError(
                "Unexpected fermionic monomial length while classifying "
                f"parent {term.index}: {n_operators} operators."
            )

    if len(ranks) != 1:
        raise ValueError(
            "Could not assign a unique 1B/2B rank to parent "
            f"{term.index}; observed ranks={sorted(ranks)}."
        )

    return next(iter(ranks))


def selective_schedule_seed(
    base_seed: int,
    case_id: str,
    schedule: str,
    sample_index: int,
) -> int:
    schedule_code = {
        SHUFFLE_1B: 101,
        SHUFFLE_2B: 102,
    }[schedule]

    case_code = zlib.crc32(case_id.encode("utf-8")) & 0xFFFFFFFF
    sequence = np.random.SeedSequence(
        [base_seed, case_code, schedule_code, sample_index]
    )
    return int(sequence.generate_state(1, dtype=np.uint32)[0])


def randomized_selected_body_block_order(
    *,
    buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
    body_rank_by_parent: Sequence[int],
    selected_body_rank: int,
    rng: np.random.Generator,
) -> tuple[list[int], list[int]]:
    """
    Shuffle only 1B or only 2B parent-owned Pauli blocks.

    The body-rank slot pattern is fixed. For example:

        reference:   1B-A  2B-A  2B-B  1B-B  2B-C
        shuffle 2B:  1B-A  2B-C  2B-A  1B-B  2B-B

    Thus 1B blocks never move during the 2B ablation, and vice versa.
    """
    permuted_buckets = [
        (int(parent), list(map(int, bucket)))
        for parent, bucket in buckets
    ]

    selected_slots = [
        slot
        for slot, (parent, _) in enumerate(permuted_buckets)
        if body_rank_by_parent[int(parent)] == selected_body_rank
    ]

    if len(selected_slots) < 2:
        # Nothing meaningful to permute in this rank.
        return (
            base.flatten_reference_buckets(permuted_buckets, fallback),
            [int(parent) for parent, _ in permuted_buckets],
        )

    selected_blocks = [
        permuted_buckets[slot]
        for slot in selected_slots
    ]

    # Avoid counting the identity permutation as a randomized ablation sample.
    identity = np.arange(len(selected_blocks), dtype=np.int64)
    while True:
        permutation = rng.permutation(len(selected_blocks))
        if not np.array_equal(permutation, identity):
            break

    shuffled_blocks = [
        selected_blocks[int(index)]
        for index in permutation
    ]

    for slot, block in zip(selected_slots, shuffled_blocks):
        permuted_buckets[slot] = block

    original_rank_pattern = [
        body_rank_by_parent[int(parent)]
        for parent, _ in buckets
    ]
    shuffled_rank_pattern = [
        body_rank_by_parent[int(parent)]
        for parent, _ in permuted_buckets
    ]

    if original_rank_pattern != shuffled_rank_pattern:
        raise RuntimeError(
            "Selective shuffle changed the 1B/2B slot pattern."
        )

    # Extra check: every unselected parent must remain in the same slot.
    for slot, ((old_parent, _), (new_parent, _)) in enumerate(
        zip(buckets, permuted_buckets)
    ):
        old_parent = int(old_parent)
        new_parent = int(new_parent)
        if (
            body_rank_by_parent[old_parent] != selected_body_rank
            and old_parent != new_parent
        ):
            raise RuntimeError(
                "Selective shuffle moved an unselected parent: "
                f"slot={slot}, old={old_parent}, new={new_parent}."
            )

    indices = base.flatten_reference_buckets(
        permuted_buckets,
        fallback,
    )
    parent_order = [
        int(parent)
        for parent, _ in permuted_buckets
    ]

    return indices, parent_order


def build_schedule_specs(
    case_id: str,
    args: Any,
    raw_pauli_keys: Sequence[Any],
    final_coefficients: dict[Any, complex],
    n_qubits: int,
    hermitian_terms: Sequence[Any],
    fermion_to_pauli_indices: Sequence[Sequence[int]],
    signed_parent_order: Sequence[int],
    magnitude_parent_order: Sequence[int],
    reference_buckets: Sequence[tuple[int, Sequence[int]]],
    fallback: Sequence[int],
) -> list[Any]:
    """
    Reuse the repository implementation for the signed reference, then add
    only the two selective randomized body-rank schedules.
    """
    reference_args = copy.copy(args)
    reference_args.schedules = [REFERENCE]

    specs = _ORIGINAL_BUILD_SCHEDULE_SPECS(
        case_id=case_id,
        args=reference_args,
        raw_pauli_keys=raw_pauli_keys,
        final_coefficients=final_coefficients,
        n_qubits=n_qubits,
        hermitian_terms=hermitian_terms,
        fermion_to_pauli_indices=fermion_to_pauli_indices,
        signed_parent_order=signed_parent_order,
        magnitude_parent_order=magnitude_parent_order,
        reference_buckets=reference_buckets,
        fallback=fallback,
    )

    body_rank_by_parent = [
        fermionic_body_rank(term)
        for term in hermitian_terms
    ]

    n_1b_terms = sum(rank == 1 for rank in body_rank_by_parent)
    n_2b_terms = sum(rank == 2 for rank in body_rank_by_parent)
    n_1b_owned_blocks = sum(
        body_rank_by_parent[int(parent)] == 1
        for parent, _ in reference_buckets
    )
    n_2b_owned_blocks = sum(
        body_rank_by_parent[int(parent)] == 2
        for parent, _ in reference_buckets
    )

    print(
        "  Body-rank parents: "
        f"1B={n_1b_terms} ({n_1b_owned_blocks} owned blocks), "
        f"2B={n_2b_terms} ({n_2b_owned_blocks} owned blocks)"
    )

    for schedule, selected_rank in (
        (SHUFFLE_1B, 1),
        (SHUFFLE_2B, 2),
    ):
        if schedule not in args.schedules:
            continue

        for sample_index in range(args.samples):
            seed = selective_schedule_seed(
                args.seed,
                case_id,
                schedule,
                sample_index,
            )
            rng = np.random.default_rng(seed)

            indices, parent_order = randomized_selected_body_block_order(
                buckets=reference_buckets,
                fallback=fallback,
                body_rank_by_parent=body_rank_by_parent,
                selected_body_rank=selected_rank,
                rng=rng,
            )

            # Same final JW Hamiltonian: only the Pauli schedule is permuted.
            base.validate_pauli_order(
                schedule,
                [raw_pauli_keys[index] for index in indices],
                raw_pauli_keys,
            )

            specs.append(
                base.ScheduleSpec(
                    name=schedule,
                    sample_index=sample_index,
                    random_seed=seed,
                    pauli_order_indices=tuple(indices),
                    parent_order=tuple(parent_order),
                )
            )

    return specs


def ordering_definition(name: str) -> str:
    if name == SHUFFLE_1B:
        return (
            "signed-reference parent-owned Pauli blocks retained internally; "
            "only 1B fermionic parent blocks are randomly permuted among the "
            "original 1B block slots; all 2B blocks remain fixed"
        )

    if name == SHUFFLE_2B:
        return (
            "signed-reference parent-owned Pauli blocks retained internally; "
            "only 2B fermionic parent blocks are randomly permuted among the "
            "original 2B block slots; all 1B blocks remain fixed"
        )

    return _ORIGINAL_ORDERING_DEFINITION(name)


def configure_base_module() -> None:
    """
    Add exactly the schedules required for this ablation to the existing
    repository benchmark.
    """
    base.DETERMINISTIC_SCHEDULES = (REFERENCE,)
    base.RANDOM_SCHEDULES = SELECTIVE_RANDOM_SCHEDULES
    base.ALL_SCHEDULES = (
        base.DETERMINISTIC_SCHEDULES + base.RANDOM_SCHEDULES
    )
    base.build_schedule_specs = build_schedule_specs
    base.ordering_definition = ordering_definition


def run_base_benchmark(
    *,
    root: Path,
    panel: Sequence[tuple[str, Path]],
    args: argparse.Namespace,
) -> None:
    """
    Call the repository's existing base.main() ONCE with all 20 tensors.

    Its parser already supports repeated --tensor and its main loop already
    processes all tensor paths while preserving the existing resume logic.
    """
    configure_base_module()

    command_args: list[str] = [
        "benchmark_fermionic_structure_ablation.py",
    ]

    for _, tensor in panel:
        command_args.extend(
            ["--tensor", display_path(tensor, root)]
        )

    command_args.extend(["--steps"])
    command_args.extend(str(value) for value in args.steps)

    command_args.extend(
        [
            "--time",
            str(args.evolution_time),
            "--samples",
            str(args.samples),
            "--seed",
            str(args.seed),
            "--schedules",
            REFERENCE,
            SHUFFLE_1B,
            SHUFFLE_2B,
            "--output",
            display_path(args.output, root),
            "--tolerance",
            str(args.tolerance),
            "--parallel-threshold",
            str(args.parallel_threshold),
        ]
    )

    if args.no_spin_sector:
        command_args.append("--no-spin-sector")

    old_argv = sys.argv[:]
    try:
        sys.argv = command_args
        base.main()
    finally:
        sys.argv = old_argv


def main() -> None:
    args = parse_runner_args()
    root = repository_root()

    manifest = resolve_path_from_root(root, args.manifest)
    args.output = resolve_path_from_root(root, args.output)

    base_script = root / "analysis" / "benchmark_fermionic_structure_ablation.py"
    if not base_script.is_file():
        raise FileNotFoundError(
            "This one-file runner reuses the existing repository benchmark, "
            "but it was not found:\n"
            f"  {base_script}"
        )

    if not manifest.is_file():
        raise FileNotFoundError(
            f"Missing frozen 20-case manifest:\n  {manifest}"
        )

    panel = load_panel(manifest, root)

    print()
    print("=" * 100)
    print("Resolved exact 20-case BCH validation panel")
    print("=" * 100)
    for index, (case_id, tensor) in enumerate(panel, start=1):
        print(
            f"{index:2d}. {case_id}\n"
            f"    {display_path(tensor, root)}"
        )

    print()
    print(f"Cases:   {len(panel)}")
    print(f"Steps:   {args.steps}")
    print(f"Time:    {args.evolution_time}")
    print(f"Samples: {args.samples} per selective shuffle")
    print(f"Output:  {display_path(args.output, root)}")
    print()

    if args.dry_run:
        print("Dry run passed. No benchmark was executed.")
        return

    # Run from qhat/ so all repository-relative tensor/output paths behave
    # exactly like normal analysis commands.
    import os

    old_cwd = Path.cwd()
    try:
        os.chdir(root)
        run_base_benchmark(
            root=root,
            panel=panel,
            args=args,
        )
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    main()
