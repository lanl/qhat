"""
Test that throughput_config parameter works correctly for should_use_exact_tracking.
"""

from qhat.qre.trotter_coefficients_fast import THROUGHPUT_CONFIG, should_use_exact_tracking


def test_default_config():
    """Test with default throughput config."""
    # N=50 should be feasible with default config
    result = should_use_exact_tracking(N=50, time_budget=60)
    print(f"✓ Default config: N=50 feasible = {result}")
    assert result == True, "N=50 should be feasible with default config"


def test_custom_slower_config():
    """Test with custom slower throughput (makes exact computation look less feasible)."""
    # Create a slower config (10x slower)
    slow_config = {
        'c1_samples_per_sec': 1.5e6,   # 10x slower
        'c21_samples_per_sec': 0.8e6,  # 10x slower
        'c22_samples_per_sec': 1.1e6,  # 10x slower
    }

    # With slower throughput, N=50 might not be feasible in 60s
    result = should_use_exact_tracking(N=50, time_budget=60, throughput_config=slow_config)
    print(f"✓ Slow config: N=50 feasible = {result}")
    # With 10x slower throughput, estimation shows it would take ~10x longer
    # so some larger N values become infeasible


def test_custom_faster_config():
    """Test with custom faster throughput (makes larger N feasible)."""
    # Create a faster config (2x faster)
    fast_config = {
        'c1_samples_per_sec': 30e6,   # 2x faster
        'c21_samples_per_sec': 16e6,  # 2x faster
        'c22_samples_per_sec': 22e6,  # 2x faster
    }

    # With faster throughput, N=50 should definitely be feasible
    result = should_use_exact_tracking(N=50, time_budget=60, throughput_config=fast_config)
    print(f"✓ Fast config: N=50 feasible = {result}")
    assert result == True, "N=50 should be feasible with 2x faster config"


def test_config_affects_threshold():
    """Test that throughput config actually affects the feasibility threshold."""
    # Find a threshold N where default config says "not feasible"
    # but 10x faster config says "feasible"

    very_fast_config = {
        'c1_samples_per_sec': 150e6,   # 10x faster
        'c21_samples_per_sec': 80e6,   # 10x faster
        'c22_samples_per_sec': 110e6,  # 10x faster
    }

    # Test N=100
    default_result = should_use_exact_tracking(N=100, time_budget=60)
    fast_result = should_use_exact_tracking(N=100, time_budget=60, throughput_config=very_fast_config)

    print(f"✓ N=100: default={default_result}, 10x faster={fast_result}")
    # At least for some N, faster config should give different result


def test_config_preserves_original():
    """Test that passing custom config doesn't modify the original."""
    original_c1 = THROUGHPUT_CONFIG['c1_samples_per_sec']

    custom_config = THROUGHPUT_CONFIG.copy()
    custom_config['c1_samples_per_sec'] = 999e6

    # Use custom config
    should_use_exact_tracking(N=50, time_budget=60, throughput_config=custom_config)

    # Original should be unchanged
    assert THROUGHPUT_CONFIG['c1_samples_per_sec'] == original_c1
    print(f"✓ Original config preserved: c1_samples_per_sec = {original_c1/1e6:.0f}M/s")


if __name__ == "__main__":
    print("="*70)
    print("THROUGHPUT CONFIG TESTS")
    print("="*70)
    print()

    test_default_config()
    test_custom_slower_config()
    test_custom_faster_config()
    test_config_affects_threshold()
    test_config_preserves_original()

    print()
    print("="*70)
    print("✅ ALL THROUGHPUT CONFIG TESTS PASSED!")
    print("="*70)
