# Fixed-input ordering mechanism experiment

This directory is additive: QHAT's original physics implementation and historical
results have not been edited. Start with `PROTOCOL.md`; generated reports explain
both findings and limitations. This is a fixed-case mechanism study, not a
prospective benchmark or a claim of universal ordering superiority.

For the completed September 17 run, read `RESULTS_KO.md` first. The final
numerical report/figures are in
`runs/mechanism-t06qcgac/analysis-dl3ym9tx/`; the separately declared structural
and prediction checks are in `runs/mechanism-t06qcgac/exploratory-kdaizs_h/`.

## Run

Use `/Users/albertlee0125/miniconda3/envs/qhat/bin/python` in this directory.

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/albertlee0125/miniconda3/envs/qhat/bin/python -m unittest -v test_mechanism_experiment.py
PYTHONDONTWRITEBYTECODE=1 /Users/albertlee0125/miniconda3/envs/qhat/bin/python mechanism_experiment.py --controls
```

Every execution creates a new `runs/mechanism-*` directory; nothing is resumed or
overwritten. The full protocol is 108 baseline rows plus 21 structural controls.
Physics is single-threaded and the campaign checks a 45-minute budget before
each configuration. Do not launch duplicate campaigns while one is running.

The runner copies the three frozen input tensors from the first reproducibility
campaign. If those inputs are missing, it stops; it never regenerates integrals.
Each successful campaign includes its own input files. The local `.gitignore`
explicitly permits these small tensor archives and ignores computational caches.
No commit or push is performed automatically.

## Reading the output

- `manifest.json`: repository commit, source hashes, settings and environment.
- `inputs/`: exact tensor archives used, with hashes in the setup files.
- `*.setup.json`: Pauli/order hashes and full-space/Hamiltonian validation.
- `*.operators.json`: merged Hamiltonian, schedules and fixed parent ownership.
- `*.row*.json`, `results.json`: physical errors and decomposition measurements.
- `*.trace.json`: individual local defect powers and overlap contributions.
- `completion.json` or `failure.json`: final state of the run.

Generate figures and a Korean report with:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/albertlee0125/miniconda3/envs/qhat/bin/python analyze_mechanism.py /absolute/path/to/completed/campaign --repeat-campaign /absolute/path/to/f2/anchor/pilot
```

The report writer creates an exclusive `analysis-*` subdirectory. Open
`REPORT.md` and the three PNG figures. `analysis.json` includes validation
residuals and, if requested, independent-process repeat comparisons.

The separately declared exploratory follow-up is documented in
`EXPLORATORY_ADDENDUM.md`. After the main campaign completes, run:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/albertlee0125/miniconda3/envs/qhat/bin/python -m unittest -v test_prediction.py
PYTHONDONTWRITEBYTECODE=1 /Users/albertlee0125/miniconda3/envs/qhat/bin/python validate_structure_and_prediction.py /absolute/path/to/completed/campaign
```

This writes symbolic block-commutator checks, unfitted F2 finite-time predictions,
and three extra product-formula results at T=0.75,r=100. Predictions are saved
before the three extra states are computed. This is an exploratory check on an
already studied molecule, not independent-molecule or broadly prospective evidence.

## Interpretation warnings

E is `1-|normalized overlap|`, not standard infidelity. The internal stable
calculation uses the orthogonal residual to avoid subtracting near-unit values.
Both exact and Trotter dynamics use the same scalar energy convention.

The decomposition `I=D+X=D*C` is exact, but is not by itself a predictive theory.
`C<1` is net destructive interference, `C>1` is net constructive accumulation.
A smaller positive C can explain an ordering advantage without net negative
interference. When comparing different step counts, C generally scales with r;
compare C/r or compare orderings at the same r.

Sector-conditional errors diagnose error channels; they do not grant free
postselection or establish lower experimental cost. Random block controls are
not matched perturbation budgets and are not independent molecular samples.
