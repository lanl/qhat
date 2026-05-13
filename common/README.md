# QHAT Common Software Modules

Modules available in this directory are used throughout QHAT. They are generally not intended to
be used directly by users.

## Pauli String and Time Evolution

- **`pauli_string_evolution.py`**: Implements time evolution under a single Pauli string Hamiltonian:
  `U = exp(-i * coefficient * PauliString * t / ℏ)`. Used as a building block for Trotterization.
  - Provides both unitary matrix generation and resource estimation
  - Compatible with pyLIQTR and Qualtran

- **`commuting_pauli_string_evolution.py`**: Implements exact time evolution under a sum of
  commuting Pauli string Hamiltonians. Since the terms commute, the evolution can be exactly
  decomposed into a product of individual evolutions:
  `exp(-i * (c₁P₁ + c₂P₂ + ...) * t / ℏ) = exp(-i * c₁P₁ * t / ℏ) * exp(-i * c₂P₂ * t / ℏ) * ...`
  - Validates that all terms pairwise commute

- **`dense_pauli_exp.py`**: Lower-level implementation using DensePauliString and
  DensePauliExponential. Provides helper functions for working with Cirq gates and Qualtran bloqs.
  - Used by the Trotter implementation
  - Contains utility functions for register manipulation
  - ⚠️ **Minimal test coverage** (only example scripts, not proper unit tests)

## Trotterization

- **`trotter.py`**: QHAT implementation of Trotterization with nested bloq structure
  (RampedTrotterizedUnitary → RampedTrotterStep → TrotterRamp).
  - Supports ramped coefficients for higher-order methods
  - ⚠️ **Little to no test coverage**

## Hamiltonian Representations

- **`LCPSHamiltonian.py`**: Linear Combination of Pauli Strings (LCPS) Hamiltonian representation.
  Organizes Pauli strings into commuting groups for efficient processing.
  - Groups are containers of terms where all terms within a group commute
  - Each term consists of a Pauli string (dense format like "IXZZYI") and a coefficient
  - Validates structure and commutativity
  - Compatible with pyLIQTR ProblemInstance interface

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

- **`QPE_Kitaev.py`**: Implementation of Kitaev's "standard" quantum phase estimation circuit.
  - Based on Kitaev's original algorithm from arXiv:quant-ph/9511026
  - Builds QPE circuit using Hadamard gates, controlled unitaries, and inverse QFT
  - ⚠️ **Note:** May be superseded by Qualtran's TextbookQPE and other QPE implementations
  - Currently used for educational/exploratory purposes

## Testing Utilities

- **`test_utils.py`**: Shared testing utilities including:
  - `validate_pauli_string()`: Validates Pauli string format (only I, X, Y, Z characters)
  - Used across multiple modules to ensure consistent validation

## Example/Development Files

These files are examples or development test scripts, not production code:

- **`dps_test.py`**: Example/test script for DensePauliString and DensePauliExponential
- **`trotter_test.py`**: Example/test script for the Trotter implementation

⚠️ **Note:** These are not proper unit tests and are not integrated into the test suite.
