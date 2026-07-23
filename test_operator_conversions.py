"""
Comprehensive test of OperatorRepresentation with realistic data.
Tests all four combinations (H/U × shifted/unshifted) and diagnoses Frobenius norm issue.
"""
import numpy as np
import scipy.linalg
from qhat.analysis.operators import OperatorRepresentation

# Load realistic data
H_prime = np.load('analysis/Be-H_exact_hamiltonian.npz')['matrix']
U_prime_approx = np.load('analysis/Be-H_unitary_matrix.npz')['matrix']

# Parameters
E = 3.845617142543042  # Energy shift
t = 1.4176481405134915  # Timestep
dim = 64

print("=" * 80)
print("COMPREHENSIVE OPERATOR CONVERSION TEST")
print("=" * 80)
print(f"\nParameters:")
print(f"  Energy shift E = {E:.6f}")
print(f"  Timestep t = {t:.6f}")
print(f"  Matrix dimension = {dim}")

# ============================================================================
# PART 1: Verify what's in the files
# ============================================================================
print("\n" + "=" * 80)
print("PART 1: What's in the files?")
print("=" * 80)

H_prime_eigs = np.linalg.eigvalsh(H_prime)
print(f"\nH_prime (from file) eigenvalues:")
print(f"  Range: [{H_prime_eigs.min():.6f}, {H_prime_eigs.max():.6f}]")
print(f"  This is H' = H + E·I (shifted)")

U_prime_eigs = np.linalg.eigvals(U_prime_approx)
phases = np.angle(U_prime_eigs)
energies = -phases / t
print(f"\nU_prime_approx encoded energies:")
print(f"  Range: [{energies.min():.6f}, {energies.max():.6f}]")
print(f"  Should match H_prime: [{H_prime_eigs.min():.6f}, {H_prime_eigs.max():.6f}]")
match = np.allclose(np.sort(energies), np.sort(H_prime_eigs), atol=0.01)
print(f"  Match (±0.01): {match} ✓" if match else f"  Match: {match} ✗")

# ============================================================================
# PART 2: Test OperatorRepresentation - Create wrapper objects
# ============================================================================
print("\n" + "=" * 80)
print("PART 2: Create OperatorRepresentation wrappers")
print("=" * 80)

# Create exact operator (from H_prime)
exact_op = OperatorRepresentation(
    data=H_prime,
    operator_type='hamiltonian',
    energy_shifted=True,  # Input is H' (shifted)
    representation='dense_matrix',
    timestep=t,
    energy_shift=E
)
print("\n✓ Created exact_op (H', shifted)")

# Create approx operator (from U_prime)
approx_op = OperatorRepresentation(
    data=U_prime_approx,
    operator_type='time_evolution',
    energy_shifted=True,  # Input is U' (from H')
    representation='dense_matrix',
    timestep=t,
    energy_shift=E
)
print("✓ Created approx_op (U', shifted)")

# ============================================================================
# PART 3: Test all four forms from exact_op
# ============================================================================
print("\n" + "=" * 80)
print("PART 3: Extract all four forms from exact_op")
print("=" * 80)

# 3a. H (unshifted, physical)
H_unshifted = exact_op.get(operator_type='hamiltonian', energy_shifted=False)
H_unshifted_eigs = np.linalg.eigvalsh(H_unshifted)
print(f"\n3a. H (unshifted, physical):")
print(f"  Eigenvalues: [{H_unshifted_eigs.min():.6f}, {H_unshifted_eigs.max():.6f}]")
print(f"  Expected: H_prime - E = [{(H_prime_eigs - E).min():.6f}, {(H_prime_eigs - E).max():.6f}]")
match = np.allclose(H_unshifted_eigs, H_prime_eigs - E, atol=1e-10)
print(f"  Match: {match} {'✓' if match else '✗'}")

# 3b. H' (shifted) - should match input
H_shifted = exact_op.get(operator_type='hamiltonian', energy_shifted=True)
H_shifted_eigs = np.linalg.eigvalsh(H_shifted)
print(f"\n3b. H' (shifted) - should match input:")
print(f"  Eigenvalues: [{H_shifted_eigs.min():.6f}, {H_shifted_eigs.max():.6f}]")
print(f"  Expected: [{H_prime_eigs.min():.6f}, {H_prime_eigs.max():.6f}]")
match = np.allclose(H_shifted_eigs, H_prime_eigs, atol=1e-10)
print(f"  Match: {match} {'✓' if match else '✗'}")

# 3c. U (unshifted, physical) - from H
U_unshifted_from_H = exact_op.get(operator_type='time_evolution', energy_shifted=False)
U_unshifted_eigs = np.linalg.eigvals(U_unshifted_from_H)
phases = np.angle(U_unshifted_eigs)
energies = -phases / t
print(f"\n3c. U (unshifted, physical) - from H:")
print(f"  Encoded energies: [{energies.min():.6f}, {energies.max():.6f}]")
print(f"  Should match H_unshifted: [{H_unshifted_eigs.min():.6f}, {H_unshifted_eigs.max():.6f}]")
# For negative energies, phase wrapping is complex, so just check a few eigenvalues
print(f"  Sample energies (sorted): {np.sort(energies)[:5]}")
print(f"  Sample H eigenvalues: {np.sort(H_unshifted_eigs)[:5]}")

# 3d. U' (shifted) - from H'
U_shifted_from_H = exact_op.get(operator_type='time_evolution', energy_shifted=True)
U_shifted_eigs = np.linalg.eigvals(U_shifted_from_H)
phases = np.angle(U_shifted_eigs)
energies = -phases / t
print(f"\n3d. U' (shifted) - from H':")
print(f"  Encoded energies: [{energies.min():.6f}, {energies.max():.6f}]")
print(f"  Should match H_prime: [{H_prime_eigs.min():.6f}, {H_prime_eigs.max():.6f}]")

# ============================================================================
# PART 4: Test all four forms from approx_op
# ============================================================================
print("\n" + "=" * 80)
print("PART 4: Extract all four forms from approx_op")
print("=" * 80)

# 4a. H' (shifted) - extracted from U'
H_shifted_from_U = approx_op.get(operator_type='hamiltonian', energy_shifted=True)
H_shifted_from_U_eigs = np.linalg.eigvalsh(H_shifted_from_U)
print(f"\n4a. H' (shifted) - extracted from U':")
print(f"  Eigenvalues: [{H_shifted_from_U_eigs.min():.6f}, {H_shifted_from_U_eigs.max():.6f}]")
print(f"  Expected (from U_prime): ~[{H_prime_eigs.min():.6f}, {H_prime_eigs.max():.6f}]")

# 4b. H (unshifted) - extracted from U'
H_unshifted_from_U = approx_op.get(operator_type='hamiltonian', energy_shifted=False)
H_unshifted_from_U_eigs = np.linalg.eigvalsh(H_unshifted_from_U)
print(f"\n4b. H (unshifted) - extracted from U':")
print(f"  Eigenvalues: [{H_unshifted_from_U_eigs.min():.6f}, {H_unshifted_from_U_eigs.max():.6f}]")
print(f"  Expected: ~[{(H_prime_eigs - E).min():.6f}, {(H_prime_eigs - E).max():.6f}]")

# 4c. U' (shifted) - should match input
U_shifted_approx = approx_op.get(operator_type='time_evolution', energy_shifted=True)
print(f"\n4c. U' (shifted) - should match input:")
norm_diff = np.linalg.norm(U_shifted_approx - U_prime_approx, 'fro')
print(f"  ||retrieved - input||_F = {norm_diff:.6e}")
match = norm_diff < 1e-10
print(f"  Match: {match} {'✓' if match else '✗'}")

# 4d. U (unshifted) - this is what we compare
U_unshifted_approx = approx_op.get(operator_type='time_evolution', energy_shifted=False)
U_unshifted_approx_eigs = np.linalg.eigvals(U_unshifted_approx)
phases = np.angle(U_unshifted_approx_eigs)
energies = -phases / t
print(f"\n4d. U (unshifted) - this is what we compare:")
print(f"  Encoded energies: [{energies.min():.6f}, {energies.max():.6f}]")
print(f"  Should match H_unshifted: [{H_unshifted_eigs.min():.6f}, {H_unshifted_eigs.max():.6f}]")

# ============================================================================
# PART 5: Direct comparison - what error_analysis does
# ============================================================================
print("\n" + "=" * 80)
print("PART 5: Direct comparison (what error_analysis does)")
print("=" * 80)

frobenius = np.linalg.norm(U_unshifted_from_H - U_unshifted_approx, 'fro')
print(f"\n||U_exact(unshifted) - U_approx(unshifted)||_F = {frobenius:.6f}")
print(f"  (This is what error_analysis reports)")

# ============================================================================
# PART 6: Alternative comparisons to diagnose the issue
# ============================================================================
print("\n" + "=" * 80)
print("PART 6: Alternative comparisons")
print("=" * 80)

# 6a. Compare on shifted scale
frobenius_shifted = np.linalg.norm(U_shifted_from_H - U_shifted_approx, 'fro')
print(f"\n6a. Shifted scale:")
print(f"  ||U'_exact - U'_approx||_F = {frobenius_shifted:.6f}")

# 6b. Direct computation without OperatorRepresentation
print(f"\n6b. Direct computation (no OperatorRepresentation):")

# Compute U_exact directly from H_physical
H_physical = H_prime - E * np.eye(dim)
U_exact_direct = scipy.linalg.expm(-1j * H_physical * t)

# Compute U_approx by phase unshifting U_prime
phase = np.exp(1j * E * t)
U_approx_direct = phase * U_prime_approx

frobenius_direct = np.linalg.norm(U_exact_direct - U_approx_direct, 'fro')
print(f"  ||U_exact(direct) - U_approx(direct)||_F = {frobenius_direct:.6f}")

# 6c. Compare OperatorRepresentation result vs direct
print(f"\n6c. Consistency check:")
diff_exact = np.linalg.norm(U_unshifted_from_H - U_exact_direct, 'fro')
diff_approx = np.linalg.norm(U_unshifted_approx - U_approx_direct, 'fro')
print(f"  ||U_exact(OpRep) - U_exact(direct)||_F = {diff_exact:.6e}")
print(f"  ||U_approx(OpRep) - U_approx(direct)||_F = {diff_approx:.6e}")

if diff_exact > 1e-10 or diff_approx > 1e-10:
    print(f"  ⚠ WARNING: OperatorRepresentation doesn't match direct computation!")
else:
    print(f"  ✓ OperatorRepresentation matches direct computation")

# ============================================================================
# PART 7: Check unitarity
# ============================================================================
print("\n" + "=" * 80)
print("PART 7: Check unitarity")
print("=" * 80)

identity = np.eye(dim)
unitarity_U_unshifted = np.linalg.norm(U_unshifted_from_H.conj().T @ U_unshifted_from_H - identity, 'fro')
unitarity_U_approx = np.linalg.norm(U_unshifted_approx.conj().T @ U_unshifted_approx - identity, 'fro')

print(f"\nUnitarity check:")
print(f"  ||U_exact† U_exact - I||_F = {unitarity_U_unshifted:.6e}")
print(f"  ||U_approx† U_approx - I||_F = {unitarity_U_approx:.6e}")

if unitarity_U_unshifted > 1e-10:
    print(f"  ⚠ WARNING: U_exact is not unitary!")
if unitarity_U_approx > 1e-10:
    print(f"  ⚠ WARNING: U_approx is not unitary!")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"\nFrobenius norm comparisons:")
print(f"  Via OperatorRepresentation (unshifted): {frobenius:.6f}")
print(f"  Via OperatorRepresentation (shifted):   {frobenius_shifted:.6f}")
print(f"  Direct computation (no OpRep):          {frobenius_direct:.6f}")
print(f"\nExpected (from earlier tests): ~0.038")
print(f"\nConsistency:")
print(f"  OperatorRepresentation matches direct: {diff_exact < 1e-10 and diff_approx < 1e-10}")
print(f"  Both operators are unitary: {unitarity_U_unshifted < 1e-10 and unitarity_U_approx < 1e-10}")

if frobenius > 0.1:
    print(f"\n⚠ WARNING: Frobenius error ({frobenius:.6f}) is unexpectedly high!")
    print(f"  Expected: ~0.038")
    print(f"  Possible causes:")
    print(f"    - Incorrect unshifting")
    print(f"    - Phase wrapping issues for negative energies")
    print(f"    - Eigenvalue extraction issues")
else:
    print(f"\n✓ Frobenius error looks reasonable")
