# New-molecule validation and novelty audit

Start with [RESULTS_KO.md](RESULTS_KO.md) and [NOVELTY_KO.md](NOVELTY_KO.md).
The prospective design is [PROTOCOL.md](PROTOCOL.md).

New files only; existing QHAT source and previous campaigns are untouched.
The helper directory `../codex_mechanism_experiment` must be preserved alongside
this directory. A matching Python environment with PySCF and OpenFermion is needed.

Run from this directory:

```sh
/Users/albertlee0125/miniconda3/envs/qhat/bin/python -B -m unittest -v test_external_validation
/Users/albertlee0125/miniconda3/envs/qhat/bin/python -B external_validation.py
```

Fresh runs create unique directories without overwriting outputs. Regenerating
degenerate orbitals on different library/platform versions need not give identical
tensors; use frozen inputs for numerical reproduction:

```sh
/Users/albertlee0125/miniconda3/envs/qhat/bin/python -B external_validation.py --inputs-from runs/transfer-vtrgy3ds
/Users/albertlee0125/miniconda3/envs/qhat/bin/python -B analyze_transfer.py runs/transfer-vtrgy3ds --repeat runs/transfer-hux1k985
```

## Completed artifacts

- Main: `runs/transfer-vtrgy3ds`, 272/272 target rows.
- Independent process repeat: `runs/transfer-hux1k985`, same frozen tensors, 272/272.
- Analysis: `runs/transfer-vtrgy3ds/analysis-x8cr3me0`.
- Predictions before target evolution:
  `runs/transfer-vtrgy3ds/predictions_before_target_evolution.json`.
- Frozen prediction SHA256:
  `b8140c1255cd32392a9625a96141a00f64176e1d943e97f3712ae7eed2c349ea`.
- `transfer_validation.png`: inspected two-panel figure (predicted vs actual I;
  all eight anchor error ratios, including fermionic losses).

Two analyses exist because an initial bitwise repeat assertion failed on tiny
floating-point differences. The final report explicitly uses numerical tolerance;
it does not claim bitwise-equal outcome JSON. Predictions and input/order hashes
are identical. No scientific success criterion was changed.

## Tests run

2026-09-17: `test_external_validation` 4/4; previous `test_mechanism_experiment`
and `test_prediction` 13/13. Each campaign also validates SCF energy conversion,
full-space execution, Hamiltonian action, derivative convergence, state norm,
defect reconstruction, source/input immutability, and complete row counts.
