# Issue #31 tensor-stage threshold hotfix

This is a narrow QHAT-only change.  OpenFermion is not edited or monkey-patched.

## What changed

The Hamiltonian configuration now accepts a finite, non-negative cutoff:

```python
hamiltonian.coefficient_threshold = 1e-8  # default
hamiltonian.coefficient_threshold = 0.0   # no tensor-stage filtering
```

QHAT uses a local `spinorb_from_spatial` function when it expands active-space spatial integrals
into spin-orbital tensors.  The rest of the generator stays on its original path:

- OpenFermion's standard Jordan-Wigner and Bravyi-Kitaev functions are unchanged.
- The active-space pickle behavior is unchanged.
- The `.tensors.npz` format remains `constant`, `one_body`, and `two_body`.
- Metadata records only the configured spin-orbital threshold.

## Limitations

OpenFermion's mapping functions can still remove small coefficients.  Setting this option to
`0.0` therefore disables filtering only during tensor construction, not in the final Pauli sum.
Pauli-level thresholding is a separate follow-up.

An existing `[stub]_[astag].pickle` is reused without checking its threshold.  When changing the
setting, remove that active-space pickle or use a new `general.file_stub`.  Removing only the
`.tensors.npz` file is not enough.

## Li2 check

For the bundled five-spatial-orbital Li2 case:

| Tensor threshold | One-body entries | Two-body entries | Stock JW terms | Stock BK terms |
| ---: | ---: | ---: | ---: | ---: |
| `1e-8` | 14 | 436 | 156 | 156 |
| `0.0` | 50 | 2500 | 156 | 156 |

The tensor counts confirm that the setting works at the intended stage.  The unchanged mapped
term counts confirm that final Pauli thresholding remains outside this hotfix.
