# Issue #31 thresholding hotfix

This is a QHAT-only hotfix.  It does not edit or monkey-patch OpenFermion.

## What changed

Add one setting to a Hamiltonian-generator configuration:

```python
# Existing/default behavior
hamiltonian.coefficient_threshold = 1e-8

# Disable magnitude thresholding
hamiltonian.coefficient_threshold = 0.0
```

QHAT now owns two small operations that previously had hidden OpenFermion cutoffs:

1. Spatial integrals are expanded into spin-orbital tensors in
   `hamiltonian_generator/thresholding.py`, using the configured threshold.
2. JW/BK mapping is done one unit-coefficient fermion term at a time.  QHAT restores the actual
   coefficients, sums matching Pauli strings, and applies the configured threshold only once at
   the end.

The second step matters because OpenFermion also removes small coefficients while incrementally
building a `QubitOperator`.  Fixing only the tensor conversion does not disable thresholding in the
final Pauli file.

## Cache handling

The active-space pickle already contains thresholded tensors.  QHAT now records the threshold on
that object and checks it when loading the cache.  If the requested value differs, QHAT reloads the
Hartree-Fock pickle and rebuilds the active-space tensors.  Old pickles are treated as `1e-8`, which
was the only behavior before this hotfix.

The output filenames are unchanged.  Running the same `general.file_stub` with a new threshold
overwrites the active-space and final files, so use another stub if both versions must be kept.

This hotfix is intentionally limited to `hamiltonian_generator/hamgen.py`.  Its generated `.dat`
Pauli file follows the configured threshold.  The separate analysis path that remaps a
`.tensors.npz` file is unchanged and still uses OpenFermion's normal mapping behavior.

## Li2 result

Using the bundled five-spatial-orbital Li2 fixture:

| Threshold | One-body tensor entries | Two-body tensor entries | JW/BK Pauli terms | Ground energy (Ha) |
| ---: | ---: | ---: | ---: | ---: |
| `1e-8` | 14 | 436 | 156 | about `-14.6426835699246` |
| `0.0` | 50 | 2500 | 876 | about `-14.6426835699246` |

This reproduces the Issue #31 observation: `1e-8` removes 36 one-body and 2,064 two-body tensor
entries.  The one-norm of the difference between the final Hamiltonians is only about `1.30e-14`
Ha, so their ground energies are numerically equivalent for this fixture.

## Important zero-threshold caveat

`0.0` keeps floating-point symmetry and cancellation residue as well as intentionally small terms.
The extra Li2 terms are mostly machine-scale noise, so zero does **not** guarantee identical term
counts at every bond length.  A small positive threshold is normally more stable; zero is useful
for auditing and sensitivity checks.

For `m` dense active spatial orbitals, the threshold-free spin tensors have `2*m^2` one-body and
`4*m^4` two-body entries.  For this Li2 case (`m = 5`), that gives 50 and 2,500.

## Files in the hotfix

- `hamiltonian_generator/hamgen.py`: uses the configurable QHAT path and validates cached tensors.
- `hamiltonian_generator/hamgen_types.py`: exposes and validates the configuration setting.
- `hamiltonian_generator/thresholding.py`: contains the local tensor conversion and JW/BK mapper.
- `hamiltonian_generator/config.py`: shows the setting in the example configuration.
- `hamiltonian_generator/README.md`: documents normal usage.
