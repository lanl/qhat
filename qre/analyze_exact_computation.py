"""
Analyze the feasibility of exact computation with early termination.

For small systems, we can compute exact values by tracking which
pairs/triples we've sampled. Once we've seen all of them, we have
the exact answer and can stop early.
"""

import math
import numpy as np

print("="*70)
print("EXACT COMPUTATION FEASIBILITY ANALYSIS")
print("="*70)

print("\nFor a Hamiltonian with N terms, the total number of combinations is:")
print("  C1:  N*(N-1)/2 pairs")
print("  C21: Σₖ C(N-k-1, 2) triples")
print("  C22: N*(N-1)/2 pairs")

print("\nMemory required to track which combinations we've seen:")
print("  Using a set of tuples: ~8 bytes per entry (optimistic)")
print("  Using a numpy boolean array: 1 bit per entry")

print("\n" + "="*70)
print("MEMORY REQUIREMENTS")
print("="*70)

print(f"\n{'N':>5} | {'C1 pairs':>10} | {'C21 triples':>12} | {'C22 pairs':>10} | "
      f"{'Memory (set)':>15} | {'Memory (bits)':>15}")
print("-"*95)

for N in [10, 20, 50, 100, 200, 500, 1000]:
    c1_count = N * (N - 1) // 2
    c21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    c22_count = N * (N - 1) // 2

    total_combinations = c1_count + c21_count + c22_count

    # Memory for sets (approximate: 8 bytes per entry for pair, 12 for triple)
    mem_set = (c1_count * 8 + c21_count * 12 + c22_count * 8) / (1024**2)

    # Memory for bit arrays
    mem_bits = (c1_count + c21_count + c22_count) / (8 * 1024**2)

    print(f"{N:5} | {c1_count:10,} | {c21_count:12,} | {c22_count:10,} | "
          f"{mem_set:12.2f} MB | {mem_bits:12.2f} MB")

print("\n" + "="*70)
print("PERFORMANCE IMPACT ANALYSIS")
print("="*70)

print("""
Checking if a sample was already seen:
  - Python set: O(1) average case, ~100 ns per check
  - NumPy boolean array: O(1), ~50 ns per check (if we can index directly)

For a batch of 10,000 samples:
  - Set checks: ~1 ms overhead
  - Array checks: ~0.5 ms overhead

Compared to batch processing time (~10-50 ms), this is <5% overhead.
""")

print("="*70)
print("WHEN IS EXACT COMPUTATION FASTER?")
print("="*70)

print("\nWith current performance (fast version):")
print("  - C1:  ~15M samples/second")
print("  - C21: ~8M samples/second")
print("  - C22: ~11M samples/second")

print("\nTime to sample all combinations (if we get lucky and no duplicates):")
print(f"\n{'N':>5} | {'C1 time':>10} | {'C21 time':>10} | {'C22 time':>10} | {'Total':>10}")
print("-"*60)

for N in [10, 20, 50, 100, 200]:
    c1_count = N * (N - 1) // 2
    c21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    c22_count = N * (N - 1) // 2

    t_c1 = c1_count / 15e6  # 15M samples/sec
    t_c21 = c21_count / 8e6   # 8M samples/sec
    t_c22 = c22_count / 11e6  # 11M samples/sec

    total = t_c1 + t_c21 + t_c22

    print(f"{N:5} | {t_c1:8.3f}s | {t_c21:8.3f}s | {t_c22:8.3f}s | {total:8.3f}s")

print("\nNote: This assumes no duplicate samples, which is optimistic.")
print("In practice, with random sampling, we need ~N*ln(N) samples to")
print("cover all N items (coupon collector problem).")

print("\n" + "="*70)
print("REALISTIC EXACT COMPUTATION TIME")
print("="*70)

print("\nAccounting for duplicate sampling (coupon collector):")
print("Expected samples needed ≈ N * ln(N) for pairs")
print("Expected samples needed ≈ N * (ln(N))^1.5 for triples")

print(f"\n{'N':>5} | {'C1 time':>10} | {'C21 time':>10} | {'C22 time':>10} | {'Total':>10}")
print("-"*60)

for N in [10, 20, 50, 100, 200]:
    c1_count = N * (N - 1) // 2
    c21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
    c22_count = N * (N - 1) // 2

    # Coupon collector: need N*ln(N) samples on average
    c1_samples = c1_count * np.log(c1_count) if c1_count > 0 else 0
    c21_samples = c21_count * np.log(c21_count) if c21_count > 0 else 0
    c22_samples = c22_count * np.log(c22_count) if c22_count > 0 else 0

    t_c1 = c1_samples / 15e6
    t_c21 = c21_samples / 8e6
    t_c22 = c22_samples / 11e6

    total = t_c1 + t_c21 + t_c22

    print(f"{N:5} | {t_c1:8.3f}s | {t_c21:8.3f}s | {t_c22:8.3f}s | {total:8.3f}s")

print("\n" + "="*70)
print("RECOMMENDATION")
print("="*70)

print("""
Early termination is HIGHLY BENEFICIAL for:
✓ N ≤ 50: Exact computation completes in < 5 seconds
✓ N ≤ 100: Exact computation completes in < 30 seconds
✓ Memory usage is negligible (< 10 MB for N=100)

Implementation strategy:
1. Track seen pairs/triples in Python sets (simple, fast enough)
2. After each batch, check if len(seen) == total_possible
3. If complete, compute exact sum from accumulated values and return

Overhead:
- Memory: < 1 MB for N < 50, < 10 MB for N < 100
- Performance: < 5% overhead (set membership checks)
- Benefit: Exact results + potential early termination

For N > 200:
- C21 memory becomes large (>50 MB)
- Exact computation takes longer than 60s anyway
- Can still implement but less beneficial

VERDICT: Implement for all sizes, with optional disable for very large N.
""")

print("\n" + "="*70)
print("CONCRETE EXAMPLE")
print("="*70)

N = 50
c1_count = N * (N - 1) // 2
c21_count = sum(math.comb(N - k - 1, 2) for k in range(N - 1))
c22_count = N * (N - 1) // 2

print(f"\nFor N={N} Pauli terms:")
print(f"  C1:  {c1_count:,} pairs")
print(f"  C21: {c21_count:,} triples")
print(f"  C22: {c22_count:,} pairs")

c1_samples = c1_count * np.log(c1_count)
c21_samples = c21_count * np.log(c21_count)
c22_samples = c22_count * np.log(c22_count)

t_c1 = c1_samples / 15e6
t_c21 = c21_samples / 8e6
t_c22 = c22_samples / 11e6

print(f"\nExpected samples to complete:")
print(f"  C1:  {c1_samples:,.0f} samples (~{t_c1:.2f}s)")
print(f"  C21: {c21_samples:,.0f} samples (~{t_c21:.2f}s)")
print(f"  C22: {c22_samples:,.0f} samples (~{t_c22:.2f}s)")
print(f"  Total: ~{t_c1 + t_c21 + t_c22:.2f}s for EXACT result")

print(f"\nWith 60s time budget and Monte Carlo:")
print(f"  Gets ~{15e6 * 20:.0e} C1 samples (Monte Carlo)")
print(f"  Gets ~{8e6 * 20:.0e} C21 samples (Monte Carlo)")
print(f"  Gets ~{11e6 * 20:.0e} C22 samples (Monte Carlo)")

print(f"\nWith early termination:")
print(f"  Gets EXACT result in ~{t_c1 + t_c21 + t_c22:.1f}s")
print(f"  Saves ~{60 - (t_c1 + t_c21 + t_c22):.1f}s")
print(f"  ✅ Significantly better!")

print("\n" + "="*70)
