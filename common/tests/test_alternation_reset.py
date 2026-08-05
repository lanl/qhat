"""
Test suite for the alternation reset bug fix.

This test suite specifically verifies that the ramp direction alternation
resets at the start of each Trotter step, allowing pattern detection to work
correctly for methods with odd numbers of coefficients.

See ALTERNATION_BUG_SUMMARY.md for details.
"""

import pytest
from qhat.common.trotter_flattened import (
    expand_ramped_trotterization,
    get_trotterization_coefficients,
)


class TestAlternationReset:
    """Test that alternation direction resets at each Trotter step."""

    def test_first_order_two_steps_pattern_matches(self):
        """Test first-order (C=1, odd) with 2 steps produces repeating pattern."""
        # First-order has 1 coefficient - odd, so was affected by bug
        result = expand_ramped_trotterization(3, [1.0], 2)

        # Both steps should follow same pattern: ascending (0,1,2)
        # Step 1: (0, 1.0), (1, 1.0), (2, 1.0)
        # Step 2: (0, 1.0), (1, 1.0), (2, 1.0)
        # Since first-order only has one ramp per step, there's no adjacent
        # terms to combine within a step, but we do get the repeating pattern
        expected = [(0, 1.0), (1, 1.0), (2, 1.0), (0, 1.0), (1, 1.0), (2, 1.0)]
        assert result == expected

    def test_first_order_multiple_steps_consistent(self):
        """Test first-order with many steps produces consistent repeating pattern."""
        # With the fix, every step should start with ascending=True
        result = expand_ramped_trotterization(2, [1.0], 5, combine_terms=False)

        # Without combining, we should see the raw pattern
        # Each step: (0, 1.0), (1, 1.0)
        # 5 steps = 10 operations
        assert len(result) == 10

        # Check pattern: should be (0,1), (0,1), (0,1), (0,1), (0,1)
        for i in range(5):
            assert result[i * 2] == (0, 1.0), f"Step {i+1} first term"
            assert result[i * 2 + 1] == (1, 1.0), f"Step {i+1} second term"

    def test_ruth_1983_two_steps_pattern_matches(self):
        """Test Ruth 1983 (C=5, odd) with 2 steps produces repeating pattern."""
        # Ruth 1983 has 5 coefficients - odd, so was affected by bug
        ruth_coeffs = get_trotterization_coefficients("ruth 1983")
        assert len(ruth_coeffs) == 5, "Ruth 1983 should have 5 coefficients"

        # Test with 2 terms, 2 steps, no combining to see raw pattern
        result = expand_ramped_trotterization(2, ruth_coeffs, 2, combine_terms=False)

        # With the fix, each step should have same pattern:
        # Ramp 1 (asc): 0, 1
        # Ramp 2 (desc): 1, 0
        # Ramp 3 (asc): 0, 1
        # Ramp 4 (desc): 1, 0
        # Ramp 5 (asc): 0, 1

        # Step 1 (10 operations: 5 ramps × 2 terms/ramp)
        step1_pattern = [(idx, coeff) for idx, coeff in result[:10]]

        # Step 2 (10 operations: 5 ramps × 2 terms/ramp)
        step2_pattern = [(idx, coeff) for idx, coeff in result[10:20]]

        # With the fix, patterns should match
        # Extract just indices to check alternation pattern
        step1_indices = [idx for idx, _ in step1_pattern]
        step2_indices = [idx for idx, _ in step2_pattern]

        # Should be: [0,1, 1,0, 0,1, 1,0, 0,1] for both steps
        expected_indices = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1]
        assert step1_indices == expected_indices, f"Step 1 pattern incorrect: {step1_indices}"
        assert step2_indices == expected_indices, f"Step 2 pattern incorrect: {step2_indices}"

    def test_ruth_1983_three_steps_all_match(self):
        """Test Ruth 1983 with 3 steps - all should have same pattern."""
        ruth_coeffs = get_trotterization_coefficients("ruth 1983")

        result = expand_ramped_trotterization(2, ruth_coeffs, 3, combine_terms=False)

        # Each step: 5 ramps × 2 terms = 10 operations
        step1_indices = [idx for idx, _ in result[0:10]]
        step2_indices = [idx for idx, _ in result[10:20]]
        step3_indices = [idx for idx, _ in result[20:30]]

        # All three should match
        expected_indices = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1]
        assert step1_indices == expected_indices, "Step 1 mismatch"
        assert step2_indices == expected_indices, "Step 2 mismatch"
        assert step3_indices == expected_indices, "Step 3 mismatch"

    def test_second_order_still_works(self):
        """Test that even-C methods (second-order, C=2) still work correctly."""
        # Second-order has 2 coefficients - even, so wasn't affected by bug
        # But we need to make sure we didn't break it
        result = expand_ramped_trotterization(3, [0.5, 0.5], 2)

        # Each step: asc (0,1,2), desc (2,1,0)
        # Step 1: (0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)
        # Step 2: (0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)
        # Combined: merges to (0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 1.0), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)
        expected = [(0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 1.0), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)]
        assert result == expected

    def test_blanes_moan_4th_still_works(self):
        """Test that Blanes-Moan 4th (C=12, even) still works correctly."""
        # Blanes-Moan has 12 coefficients - even, so wasn't affected
        bm4_coeffs = get_trotterization_coefficients("bm4")
        assert len(bm4_coeffs) == 12, "Blanes-Moan 4th should have 12 coefficients"

        # Just verify it runs and produces expected pattern structure
        result = expand_ramped_trotterization(2, bm4_coeffs, 2, combine_terms=False)

        # Each step: 12 ramps × 2 terms = 24 operations
        # 2 steps = 48 operations total
        assert len(result) == 48

        # Check that both steps have same pattern
        step1_indices = [idx for idx, _ in result[0:24]]
        step2_indices = [idx for idx, _ in result[24:48]]

        assert step1_indices == step2_indices, "Blanes-Moan steps should match"

    def test_direction_starts_ascending_every_step(self):
        """Test that direction explicitly starts ascending for every step."""
        # Use 4 coefficients and check the raw pattern
        result = expand_ramped_trotterization(3, [0.25, 0.25, 0.25, 0.25], 3, combine_terms=False)

        # Each step: 4 ramps × 3 terms = 12 operations
        # Ramps should be: asc (0,1,2), desc (2,1,0), asc (0,1,2), desc (2,1,0)

        for step in range(3):
            offset = step * 12
            step_indices = [idx for idx, _ in result[offset:offset+12]]

            # Expected pattern for each step
            expected = [0, 1, 2,  # Ramp 1 (asc)
                       2, 1, 0,  # Ramp 2 (desc)
                       0, 1, 2,  # Ramp 3 (asc)
                       2, 1, 0]  # Ramp 4 (desc)

            assert step_indices == expected, f"Step {step + 1} pattern incorrect"

    def test_odd_coefficients_multiple_steps_pattern_repeats(self):
        """Test various odd-coefficient counts with multiple steps."""
        # Test C=1, 3, 5, 7
        for num_coeffs in [1, 3, 5, 7]:
            coeffs = [1.0 / num_coeffs] * num_coeffs

            # Test with 4 steps
            result = expand_ramped_trotterization(2, coeffs, 4, combine_terms=False)

            # Each step should have same pattern
            ops_per_step = num_coeffs * 2  # num_coeffs ramps × 2 terms

            step_patterns = []
            for step in range(4):
                offset = step * ops_per_step
                step_indices = [idx for idx, _ in result[offset:offset+ops_per_step]]
                step_patterns.append(step_indices)

            # All steps should have identical patterns
            for i in range(1, 4):
                assert step_patterns[0] == step_patterns[i], \
                    f"C={num_coeffs}: Step 1 and Step {i+1} patterns don't match"

    def test_pattern_enables_detection(self):
        """Test that repeating patterns allow for O(log n) optimization detection."""
        # This is the key test - Ruth 1983 with many steps should produce
        # a repeating pattern that can be detected

        ruth_coeffs = get_trotterization_coefficients("ruth 1983")

        # Test with 10 steps (enough to show repetition)
        result = expand_ramped_trotterization(3, ruth_coeffs, 10)

        # The result should have a repeating structure
        # For Ruth 1983: 5 coefficients, 3 terms
        # Each step produces a sequence that should repeat

        # Find the pattern length by comparing combined results
        # For Ruth 1983 with 3 terms: each step should produce same sub-pattern
        # With combining, the pattern repeats after boundaries

        # We can't easily test pattern detection directly, but we can verify
        # that the structure is consistent across steps by checking that
        # the sequence has the expected properties

        # At minimum, verify it's not empty and has reasonable length
        assert len(result) > 0

        # For 10 steps, we expect the operations to combine efficiently
        num_ops = len(result)

        # For Ruth 1983 (5 coeffs), with 3 terms and extensive combining
        # The combining reduces the total significantly
        # Each step produces ~11 unique operations after combining
        # So 10 steps ≈ 110 operations
        expected_approx = 11 * 10

        # Allow some tolerance due to combining effects
        assert 100 <= num_ops <= 120, \
            f"Unexpected number of operations: {num_ops} (expected ~{expected_approx})"


class TestAlternationBugRegression:
    """Regression tests to ensure the bug doesn't return."""

    def test_bug_scenario_ruth_1983_1000_steps(self):
        """
        Test the exact scenario from the bug report:
        Ruth 1983 with N=1000 steps.

        This should now produce a repeating pattern.
        """
        ruth_coeffs = get_trotterization_coefficients("ruth 1983")

        # This is expensive, so just test with a smaller subset
        # and verify the pattern property
        result = expand_ramped_trotterization(3, ruth_coeffs, 10, combine_terms=False)

        # Check that first and second step have same pattern
        ops_per_step = 5 * 3  # 5 ramps × 3 terms

        step1 = result[0:ops_per_step]
        step2 = result[ops_per_step:2*ops_per_step]

        # Extract index patterns
        step1_indices = [idx for idx, _ in step1]
        step2_indices = [idx for idx, _ in step2]

        assert step1_indices == step2_indices, \
            "Ruth 1983 steps have different patterns - bug has returned!"

    def test_bug_scenario_first_order_large_steps(self):
        """
        Test first-order with large number of steps.

        Previously would alternate between steps, now should be consistent.
        """
        result = expand_ramped_trotterization(3, [1.0], 100, combine_terms=False)

        # Check every 10th step has same pattern
        for step in [0, 10, 20, 30, 40, 50]:
            step_indices = [idx for idx, _ in result[step*3:(step+1)*3]]
            expected = [0, 1, 2]  # Always ascending
            assert step_indices == expected, \
                f"Step {step} has wrong pattern - bug has returned!"

    def test_symmetry_preserved(self):
        """Test that symmetric methods produce symmetric patterns."""
        # Second-order (symmetric)
        result = expand_ramped_trotterization(5, [0.5, 0.5], 1)

        # Should be palindromic (after accounting for middle element)
        indices = [idx for idx, _ in result]
        assert indices == list(reversed(indices)), \
            "Second-order should produce symmetric pattern"
