"""Quick benchmark to compare sample throughput."""
import numpy as np
from openfermion import QubitOperator
from trotter_coefficients import trotter_error_estimator
from trotter_coefficients_fast import trotter_error_estimator_fast

# Generate test data
np.random.seed(42)
n_terms = 100
n_qubits = 15
pauli_ops = ['X', 'Y', 'Z']
terms = []

for _ in range(n_terms):
    weight = np.random.randint(1, 4)
    qubits = np.random.choice(n_qubits, size=weight, replace=False)
    ops = np.random.choice(pauli_ops, size=weight)
    pauli_string = ' '.join([f'{op}{q}' for q, op in zip(qubits, ops)])
    coeff = np.random.randn() + 1j * np.random.randn()
    terms.append(QubitOperator(pauli_string, coeff))

print("="*70)
print("ORIGINAL IMPLEMENTATION")
print("="*70)
C1_orig, C2_orig = trotter_error_estimator(terms, 3.0, batch_size=100)
print(f"C1 = {C1_orig:.6f}, C2 = {C2_orig:.6f}")

print("\n" + "="*70)
print("FAST IMPLEMENTATION")
print("="*70)
C1_fast, C2_fast = trotter_error_estimator_fast(terms, 3.0, batch_size=10000)
print(f"C1 = {C1_fast:.6f}, C2 = {C2_fast:.6f}")
