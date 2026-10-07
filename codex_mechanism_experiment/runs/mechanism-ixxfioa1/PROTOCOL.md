# Ordering mechanism experiment — preregistered local protocol

Prepared before inspecting the new experiment's outcomes.

## Question and scope

On the fixed F2, NH3 and Li2 tensors from first-repro-y23sjf9v, determine
whether initial-HF BCH norm, loss of particle/spin sector weight, and
interference of transported local defects account for ordering-dependent
finite-time state errors. This is a mechanistic case study, not an independent
held-out test, a discovery of error interference, or proof of universal superiority.

The QHAT repository and previous results remain unchanged. All new inputs,
code snapshots, settings, hashes and results live in a new campaign directory.
No old campaign is resumed or overwritten. A failed numerical check stops the
campaign. Inputs are copied from the previously frozen campaign, not regenerated.

## Fixed design

- Baselines: fermionic signed, JW signed, JW magnitude.
- Each of F2, NH3, Li2: T = 0.25, 0.5, 1; r = 25, 50, 100, 200.
- Baseline grid: 108 rows. HF, first-order, coefficient cutoff 1e-12.
- At T=1, r=100 on each molecule: round-robin descendants, within-block
  permutation seed 20260917, and five whole-block permutations with seeds
  1101, 1102, 1103, 1104, 1105. 21 additional rows; total 129.
- Parent ownership is frozen from signed first occurrence before controls.
  All orders use exactly the same merged nonidentity Pauli terms and coefficients.
  Round robin preserves order inside each original block but breaks contiguity;
  it does not isolate a single commutator perturbation. Random controls are
  descriptive, not five independent molecules or a well-powered significance test.
- Exact defect decomposition: all 36 F2 baseline rows, six NH3/Li2 baseline
  anchors (T=1,r=100), and all 21 controls; total 63 rows.
- Other 66 rows measure final-state errors and convergence only.
- One serial process, one numerical-library thread. No hardware runs or new
  molecular-integral generation. Runtime watchdog: 45 minutes per campaign.

## Numerics and identities

The exact Hamiltonian is assembled directly from the same final nonidentity
Pauli list as the product formula. Consequently both use the same scalar
convention; no fermionic/JW identity offset mismatch is present.

The computational subspace is the affine GF(2) span of all Pauli flip masks
acting on HF. It is closed under EACH Pauli operator, not merely under H.
This is exact, retaining particle- and spin-changing intermediate/final states.
Restricted Pauli actions are checked against QHAT's original full-space actions.
The assembled H action and unitaries are checked independently on a small dense
problem; its Hermiticity and symmetry residuals are recorded for every molecule.

Let S be one product-formula step, U=exp(-i H T/r), x_j=U^j x_0,
d_j=(S-U)x_j, v_j=S^(r-1-j)d_j, and phi=x_r/||x_r||.
Then S^r x_0-x_r = sum_j v_j. Set Q=I-|phi><phi|.

Stable normalized infidelity:
  I = ||Q S^r x_0||^2 / ||S^r x_0||^2.
Our primary overlap error is E=I/(1+sqrt(1-I)) = 1-|normalized overlap|.
Do not confuse this E with standard infidelity I=1-|overlap|^2.

Projected local power:
  D = sum_j ||Q v_j||^2 / ||S^r x_0||^2
    = sum_j (||d_j||^2 - |<b_j|d_j>|^2) / ||S^r x_0||^2,
  b_j=(S^dagger)^(r-1-j)phi.
Projected interference is X=I-D and coherence factor C=I/D.
X<0 (C<1) means cancellation; X>0 means net constructive accumulation.
A lower C alone is NOT evidence of cancellation when C remains above one.

This decomposition is an identity and is not by itself a novel explanatory
theory. Compare D and C between fixed-H orders, and test structural controls.
For every decomposed row explicitly propagate the defect-sum recurrence
z_(j+1)=S z_j+d_j and compare z_r with S^r x_0-x_r. Also reconstruct complex
overlap from backward contributions; independently sum the local projected
power. Small-system tests explicitly propagate every v_j and check cross terms.

Split Q S^r x_0 into exact-spin-sector and outside components; their squared
norms sum exactly to I. Separately record particle-number/spin leakage and
conditional in-sector infidelity. Conditional metrics are diagnostics, not a
free postselection performance advantage.

The initial second-order defect coefficient is computed by differentiating
the actual ordered Pauli product at zero step size. Twice its norm matches the
repository's HF commutator norm. It is not assumed to predict finite-T error.

## Checks and analysis decisions

- Normalize metric definitions; record, do not hide, state norm drift.
- Direct/reconstructed state-vector residual <= 1e-10; overlap <= 1e-10.
- Exact final number/spin leakage must be <= 1e-18.
- Full-space one-step vs closed-subspace action residual <= 1e-10.
- Frozen tensor, baseline Hamiltonian/order hashes must equal first experiment.
- Anchor E values agree with old values at rtol=0.005, atol=2e-12; old values
  include small norm-drift/subtractive-cancellation effects. Report differences.
- E <= 1e-12 is marked below analysis floor; no ratio/slope claims there.
- Evaluate rankings, r convergence, sector error shares, and D/C contributions.
  Report failures and null controls; never select only favorable molecules.
- Independent process rerun of three F2 anchors validates new implementation
  repeatability, but does not increase independent physical sample size.

## Deliverables

Self-contained input tensors, manifest and source hashes; numerical results
and traces in JSON; a machine-readable validation summary; scientific PNG
figures and a Korean interpretation report. No automatic commit or push.
