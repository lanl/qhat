# QHAT Common Software Modules

Modules available in this directory are used throughout QHAT. They are generally not intended to
be used directly by users.

## Pauli String and Time Evolution

### Representation-Independent Pauli Strings

`pauli_string.py` provides an immutable, hashable `PauliString` value type. Construct a value from
either a dense string or a sparse mapping, then request the representation needed by the caller:

```python
from qhat.common.pauli_string import PauliString

dense = PauliString.from_dense("IIXYIZ")
sparse = PauliString.from_sparse({2: "X", 3: "Y", 5: "Z"}, num_qubits=6)

assert dense == sparse
assert dense.to_dense() == "IIXYIZ"
assert dense.to_sparse() == ((2, "X"), (3, "Y"), (5, "Z"))
assert dense.to_sparse_dict() == {2: "X", 3: "Y", 5: "Z"}
```

Coefficients remain separate so `PauliString` values can be used as Hamiltonian dictionary keys.

### Core Evolution Modules

- **`pauli_string_evolution.py`**: Implements time evolution under a single Pauli string Hamiltonian:
  `U = exp(-i * coefficient * PauliString * t / ℏ)`. Used as a building block for Trotterization.
  - Provides both unitary matrix generation and resource estimation
  - Compatible with pyLIQTR and Qualtran
  - Comprehensively tested

- **`commuting_pauli_string_evolution.py`**: Implements exact time evolution under a sum of
  commuting Pauli string Hamiltonians. Since the terms commute, the evolution can be exactly
  decomposed into a product of individual evolutions:
  `exp(-i * (c₁P₁ + c₂P₂ + ...) * t / ℏ) = exp(-i * c₁P₁ * t / ℏ) * exp(-i * c₂P₂ * t / ℏ) * ...`
  - Validates that all terms pairwise commute
  - Comprehensively tested

- **`dense_pauli_exp.py`**: Lower-level implementation using DensePauliString and
  DensePauliExponential. Provides helper functions for working with Cirq gates and Qualtran bloqs.
  - Used by the original Trotter implementation
  - Contains utility functions for register manipulation
  - ⚠️ **Minimal test coverage** (only example scripts, not proper unit tests)

## Trotterization

### Modern Implementation (Recommended)

- **`trotter_flattened.py`**: Flattened QHAT implementation of Trotterization with ramped coefficients.
  Supports first through eighth-order methods. Features:
  - Flat expansion of all Trotter steps and ramps into a single sequence
  - Optional combining of adjacent identical terms (reduces operation count, benefit varies by Hamiltonian)
  - Evolves each Hamiltonian term with the bloq from its `generate_evolution_operator` (see
    `terms/`), so the terms can be of different types
  - **Comprehensively tested**
  - **Recommended for all use cases**

### Original Implementation (Deprecated)

- **`trotter_original.py`**: ⚠️ **Deprecated; superseded by `trotter_flattened.py` and will be
  removed in a future release.** Original QHAT implementation with nested bloq structure
  (RampedTrotterizedUnitary → RampedTrotterStep → TrotterRamp). `build_ramped_trotterized_unitary`
  emits a `DeprecationWarning`. To migrate, import `build_ramped_trotterized_unitary` from
  `trotter_flattened` instead (same positional arguments; pass `combine_terms=False` to reproduce
  the original gate counts).

## Hamiltonian Representations

- **`MixedFermionBosonOperator.py`**: Represents Hamiltonians with both fermionic and bosonic
  degrees of freedom.
  - Stores fermionic terms (one-body and two-body tensors)
  - Stores bosonic terms (number coefficient)
  - Stores fermion-boson interaction tensor
  - Tracks number of fermionic and bosonic states
  - Maps to qubit representations via transformation functions

## Boson Encoding

- **`bosons_binary.py`**: Binary encoding of bosonic creation and annihilation operators.
  - Maps bosons to qubits using binary representation
  - Provides `creation()` and `annihilation()` operators
  - Calculates required qubits per bosonic state based on maximum boson occupancy
  - Includes utility functions for Pauli string to matrix conversion

## Quantum Phase Estimation

- **`QPE_Kitaev.py`**: Implementation of Kitaev's "standard" quantum phase estimation algorithm.
  - Based on Kitaev's original algorithm from arXiv:quant-ph/9511026
  - Builds QPE algorithm using Hadamard gates, controlled unitaries, and inverse QFT
  - ⚠️ **Note:** May be superseded by Qualtran's TextbookQPE and other QPE implementations
  - Currently used for educational/exploratory purposes

## Utilities

- **`pauli_utils.py`**: Shared Pauli string utility functions used by both source and test code:
  - `validate_pauli_string()`: Validates Pauli string format (only I, X, Y, Z characters)
  - `get_pauli_matrix()`: Returns matrix representation of Pauli operators
  - `pauli_string_to_matrix()`: Converts Pauli string to full matrix
  - `analytical_evolution()`: Computes analytical time evolution for validation
  - Used across multiple modules to ensure consistent validation and testing

## Example/Development Files

These files are examples or development test scripts, not production code:

- **`dps_test.py`**: Example/test script for DensePauliString and DensePauliExponential
- **`trotter_test.py`**: Example/test script for the flattened Trotter implementation

⚠️ **Note:** These are not proper unit tests and are not integrated into the test suite.

## Testing

Modules with comprehensive test coverage have corresponding test files:

```bash
# Run all tests
pytest common/tests/ -v

# Individual test suites
pytest common/tests/test_pauli_string_evolution.py -v
pytest common/tests/test_commuting_pauli_string_evolution.py -v
pytest common/tests/test_trotter_flattened.py -v
```

## Module Status Summary

| Module | Test Coverage | Status |
|--------|---------------|--------|
| `pauli_string_evolution.py` | ✅ Comprehensive | Recommended |
| `commuting_pauli_string_evolution.py` | ✅ Comprehensive | Recommended |
| `trotter_flattened.py` | ✅ Comprehensive | Recommended |
| `trotter_original.py` | ⚠️ Minimal | Deprecated |
| `dense_pauli_exp.py` | ⚠️ Minimal | Used by legacy code |
| `MixedFermionBosonOperator.py` | ⚠️ None | Production use |
| `bosons_binary.py` | ⚠️ None | Production use |
| `QPE_Kitaev.py` | ⚠️ None | Exploratory |
| `pauli_utils.py` | N/A | Utility |

## Documentation

For detailed documentation on Trotterization:
- Quick guide: [`../TROTTER_IMPLEMENTATION.md`](../TROTTER_IMPLEMENTATION.md)
- Comprehensive guide: [`TROTTERIZATION.md`](TROTTERIZATION.md)
