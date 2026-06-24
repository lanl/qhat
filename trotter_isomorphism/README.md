# Trotter Isomorphism

This package compares Jordan-Wigner (JW) and Bravyi-Kitaev (BK) Pauli
Hamiltonians that come from the same fermionic Hamiltonian.

It answers two questions:

- Are the JW and BK Pauli commutation graphs isomorphic?
- When the JW and BK product-formula terms are matched by that isomorphism, are
  the ideal Trotter errors the same?

The intended QHAT workflow is:

1. Use `hamiltonian_generator` to build an active-space fermionic Hamiltonian.
2. Convert that `InteractionOperator` to an OpenFermion `FermionOperator`.
3. Build matched JW/BK data with `build_jw_bk_analysis`.
4. Compare graph isomorphism, operator-norm Trotter error, and state-overlap
   Trotter error.

## Importing From A Checkout

From the repository checkout, put the parent directory of `qhat` on
`PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/Repos
```

For this local checkout:

```bash
export PYTHONPATH=/Users/albertlee0125/Repos
```

Then import with:

```python
from qhat.trotter_isomorphism import build_jw_bk_analysis
```

## Method Numbers

The `method` argument selects the product formula:

| Method | Meaning |
| --- | --- |
| `1` or `"first order"` | First-order Lie-Trotter formula |
| `2` or `"second order"` | Second-order symmetric/Strang formula |
| `4` or `"fourth order"` | Fourth-order Suzuki composition |
| `8` or `"eighth order"` | Eighth-order Morales-style coefficient set |

For a fixed matched term order, increasing `num_steps` decreases the time step
`dt = time / num_steps`.

## Quick Example With A FermionOperator

```python
import numpy as np
from scipy.linalg import norm
from openfermion import FermionOperator

from qhat.trotter_isomorphism import (
    build_clifford_from_z_map,
    build_jw_bk_analysis,
    compare_jw_bk_operator_norm_errors,
)

fermion_op = (
    FermionOperator((), 0.125)
    + FermionOperator(((0, 1), (1, 0)), 1.0)
    + FermionOperator(((1, 1), (0, 0)), 1.0)
    + FermionOperator(((0, 1), (0, 0)), 0.25)
    + FermionOperator(((1, 1), (1, 0)), -0.75)
)

analysis = build_jw_bk_analysis(
    fermion_op,
    n_qubits=2,
    atol=1e-10,
    remove_identity=False,
    build_matrices=True,
)

C, B, A = build_clifford_from_z_map(analysis.jw_to_bk, analysis.n_qubits)
conjugation_error = norm(C @ analysis.H_jw @ C.conj().T - analysis.H_bk, ord=2)

result = compare_jw_bk_operator_norm_errors(
    H_jw=analysis.H_jw,
    H_bk=analysis.H_bk,
    jw_terms=analysis.jw_terms,
    bk_terms=analysis.bk_terms,
    time=0.7,
    num_steps=2,
    n_qubits=analysis.n_qubits,
    method=2,
)

print(conjugation_error)
print(result["error_jw"], result["error_bk"], result["error_difference"])
```

Expected result: `conjugation_error` and `error_difference` should be close to
machine precision.

## Using QHAT Hamiltonian Generator

This example builds H2/STO-3G through QHAT's Hamiltonian Generator functions,
then compares JW and BK Trotter errors.

OpenFermion may try to cache molecule HDF5 files in its installed package
directory. When running from a restricted environment, redirect its data
directory to a writable location before calling the generator.

```python
import logging
import numpy as np
from scipy.linalg import norm
from openfermion.transforms import get_fermion_operator
import openfermion.config as of_config
import openfermion.chem.molecular_data as molecular_data

from qhat.common.logging_utils import configure_logging
from qhat.hamiltonian_generator.hamgen import compute_Hartree_Fock, apply_active_space
from qhat.hamiltonian_generator.hamgen_types import (
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    State,
)
from qhat.trotter_isomorphism import (
    build_clifford_from_z_map,
    build_jw_bk_analysis,
    compare_jw_bk_operator_norm_errors,
)

of_config.DATA_DIRECTORY = "/private/tmp/qhat_hamgen_compare"
molecular_data.DATA_DIRECTORY = "/private/tmp/qhat_hamgen_compare"

configure_logging(level="warning")
logging.getLogger("pyscf").setLevel(logging.ERROR)

general = GeneralConfigurationUser()
general.file_stub = "/private/tmp/qhat_hamgen_compare/h2_sto3g"
general.logfile = "/private/tmp/qhat_hamgen_compare/hamgen.log"

hamiltonian = HamiltonianConfiguration()
hamiltonian.add_atom("H", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.7414, 0.0, 0.0)
hamiltonian.basis = "sto-3g"
hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 2

state = State("hamgen H2 comparison", general, hamiltonian)

hf = compute_Hartree_Fock(state)
active_ham = apply_active_space(state, hf)
fermion_op = get_fermion_operator(active_ham)

analysis = build_jw_bk_analysis(
    fermion_op,
    n_qubits=active_ham.n_qubits,
    atol=1e-10,
    remove_identity=False,
    build_matrices=True,
)

C, B, A = build_clifford_from_z_map(analysis.jw_to_bk, analysis.n_qubits)
print(norm(C @ analysis.H_jw @ C.conj().T - analysis.H_bk, ord=2))

for method in (1, 2, 4):
    for num_steps in (1, 2, 5):
        result = compare_jw_bk_operator_norm_errors(
            H_jw=analysis.H_jw,
            H_bk=analysis.H_bk,
            jw_terms=analysis.jw_terms,
            bk_terms=analysis.bk_terms,
            time=0.7,
            num_steps=num_steps,
            n_qubits=analysis.n_qubits,
            method=method,
        )
        print(method, num_steps, result)
```

For H2/STO-3G with all four spin orbitals active, this produced:

```text
JW/BK Pauli vertices: 15 / 15
JW/BK noncommuting edges: 16 / 16
||C H_JW C^dag - H_BK||_2: 0.000e+00
```

and the ideal operator-norm Trotter errors agreed to floating-point precision:

```text
method  steps       JW error              BK error              |diff|
1       1      6.7541123070063658e-02  6.7541123070063672e-02  1.388e-17
1       2      1.2549311368744981e-02  1.2549311368744999e-02  1.735e-17
1       5      2.4270373388112477e-03  2.4270373388112477e-03  0.000e+00
2       1      1.2549311368744981e-02  1.2549311368744999e-02  1.735e-17
2       2      3.0869899813873191e-03  3.0869899813873191e-03  0.000e+00
2       5      4.9171238123811687e-04  4.9171238123811687e-04  0.000e+00
4       1      8.5450393644864267e-05  8.5450393644864267e-05  0.000e+00
4       2      5.2007589982979328e-06  5.2007589982979328e-06  0.000e+00
4       5      1.3217750866465288e-07  1.3217750866465288e-07  0.000e+00
```

## Comparing The Isomorphism

The isomorphism is built by expanding the original fermionic Hamiltonian into
Majorana monomials. Each Majorana monomial has both a JW Pauli image and a BK
Pauli image. The map is:

```text
JW(gamma_S) -> BK(gamma_S)
```

Use:

```python
from qhat.trotter_isomorphism import build_jw_bk_analysis, anticommutes

analysis = build_jw_bk_analysis(fermion_op, n_qubits=n_qubits)

jw_vertices = list(analysis.jw_to_bk.keys())
bk_vertices = [analysis.jw_to_bk[p] for p in jw_vertices]

for i, p in enumerate(jw_vertices):
    for q in jw_vertices[i + 1:]:
        assert anticommutes(p, q) == anticommutes(analysis.jw_to_bk[p], analysis.jw_to_bk[q])
```

`build_jw_bk_analysis` already calls this verification internally through
`verify_graph_isomorphism`.

## Comparing Operator-Norm Trotter Error

Use `compare_jw_bk_operator_norm_errors`:

```python
from qhat.trotter_isomorphism import compare_jw_bk_operator_norm_errors

result = compare_jw_bk_operator_norm_errors(
    H_jw=analysis.H_jw,
    H_bk=analysis.H_bk,
    jw_terms=analysis.jw_terms,
    bk_terms=analysis.bk_terms,
    time=1.0,
    num_steps=4,
    n_qubits=analysis.n_qubits,
    method=2,
)

print(result["error_jw"])
print(result["error_bk"])
print(result["error_difference"])
```

The dense comparison is intended for small active spaces because it constructs
`2^n x 2^n` matrices and matrix exponentials.

## Comparing State-Dependent Error

For a state-dependent diagnostic, transform the JW initial state into BK with
the Clifford/permutation matrix `C`:

```python
from qhat.trotter_isomorphism import (
    bitstring_state,
    build_clifford_from_z_map,
    compare_state_dependent_trotter_error,
)

C, B, A = build_clifford_from_z_map(analysis.jw_to_bk, analysis.n_qubits)
psi0_jw = bitstring_state([1, 1, 0, 0])

result = compare_state_dependent_trotter_error(
    psi0_jw=psi0_jw,
    C_jw_to_bk=C,
    H_jw=analysis.H_jw,
    H_bk=analysis.H_bk,
    jw_terms=analysis.jw_terms,
    bk_terms=analysis.bk_terms,
    time=0.7,
    num_steps=2,
    n_qubits=analysis.n_qubits,
    method=2,
)

print(result["amp_jw"], result["amp_bk"])
print(result["fidelity_difference"])
```

## Why The Ideal Errors Are The Same

JW and BK are different qubit encodings of the same fermionic operator. For the
standard BK transform, there is a Clifford/permutation unitary `C` such that:

```text
H_BK = C H_JW C^dag
```

The Majorana-labeled map gives a matched term ordering:

```text
P_j^JW -> P_j^BK = C P_j^JW C^dag
```

For any product formula built from that matched order,

```text
U_trotter,BK = C U_trotter,JW C^dag
U_exact,BK   = C U_exact,JW C^dag
```

Therefore the BK error matrix is just a unitary conjugation of the JW error
matrix:

```text
E_BK = U_trotter,BK - U_exact,BK
     = C (U_trotter,JW - U_exact,JW) C^dag
     = C E_JW C^dag
```

Unitary conjugation preserves operator norms, so:

```text
||E_BK||_2 = ||E_JW||_2
```

For state-dependent overlap, the BK initial state must also be transformed:

```text
|psi0>_BK = C |psi0>_JW
```

With that state convention, the overlap amplitudes agree:

```text
<psi0_JW| U_trotter,JW^dag U_exact,JW |psi0_JW>
=
<psi0_BK| U_trotter,BK^dag U_exact,BK |psi0_BK>
```

This equality is for ideal mathematical product formulas. Hardware compilation,
finite precision synthesis, measurement noise, or different term orderings can
break the practical equality.

## Running Tests

```bash
PYTHONPATH=/Users/albertlee0125/Repos python -m pytest -q trotter_isomorphism/tests
```

With OpenFermion installed, the current tests pass:

```text
6 passed
```
