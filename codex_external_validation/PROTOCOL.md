# Prospective molecular transfer check — fixed before target calculations

Date: 2026-09-17. This is a local timestamped protocol, not a public preregistration.
The mechanism and predictor were developed on F2/NH3/Li2, including exploratory
F2 time tests. No target evolution results for the molecules below were inspected
when specifying this protocol. Existing tensor-library names were inspected only
to avoid reusing its 14 molecular families.

## Fixed inputs and scope

- Neutral singlets HCN, CH2O, H2S, SiH4, chosen for distinct composition and shape,
  not for ordering performance. Two uniform coordinate scales: 1.00 and 1.25.
- Idealized geometries in Angstrom, explicitly NOT optimized equilibrium structures:
  HCN: C=(0,0,0), N=(0,0,1.16), H=(0,0,-1.06).
  CH2O: C=(0,0,0), O=(0,0,1.21), H=(+/-0.94,0,-0.59).
  H2S: S=(0,0,0), H=(+/-1.34*sin(46deg),0,1.34*cos(46deg)).
  SiH4: Si at origin; H at (a,a,a),(a,-a,-a),(-a,a,-a),(-a,-a,a), a=1.48/sqrt(3).
- Basis 6-31g; PySCF RHF, symmetry disabled, conv_tol=1e-12, max_cycle=200,
  one numerical thread. Retain canonical orbital order; fix each orbital's sign
  so its largest-magnitude AO coefficient is positive. No optimization of orbitals
  for Trotter performance. Degenerate orbital rotations remain library-dependent.
- CAS(4 electrons,5 spatial orbitals): the two highest occupied and three lowest
  vacant orbitals; all lower occupied frozen, higher vacant excluded. Thus 10
  spin orbitals, interleaved alpha/beta, HF initial state. This is a controlled
  algorithmic test, not a chemically converged prediction of the full molecules.
- Integral conversion uses the same OpenFermion-PySCF compute_integrals and
  MolecularData.get_molecular_hamiltonian convention as QHAT; no hand-written
  two-electron index or factor-of-two conversion. Preserve nuclear/core constant
  in tensor. Remove the common identity only for both compared evolutions.
- Fresh inputs and orbital metadata are saved. Check active HF expectation plus
  frozen-core/nuclear constant against PySCF RHF energy (tolerance 1e-8 Hartree).
  Abort and report SCF, conversion, or numeric failures; do not substitute cases.

## Experiments and prior criteria

1. Three unchanged fixed-Pauli-list orderings: fermionic signed/first occurrence,
   JW signed, JW descending magnitude. Pauli/fermion threshold 1e-12.
2. T in {0.25,0.75,1.0} inverse Hartree; r in {50,100,200}: 216 target rows.
3. Before ANY target evolution, save all 216 unfitted leading-error predictions
   and their hash. The model is eta'=-iH eta+B_pi exp(-iHt)HF, eta(0)=0;
   I_pred=(T/r)^2 ||Q_T eta(T)||^2. B_pi is the dt^2 local-defect coefficient.
   No Trotter outputs enter predictions. Setup one-step checks at dt=.01 and
   local derivative checks at dt=.002,.001 are numerical validations, not target
   outcomes. They are disclosed and not used to fit or select inputs.
4. Primary metric I=1-|overlap|^2; also E=1-|overlap|. Stable normalized orthogonal
   residual evaluation. Rows with E<=1e-12 are reported but excluded from relative
   accuracy and rank judgments. Predeclared prediction accuracy criterion:
   <=5% relative I error for EACH eligible row at r=100 or 200. Report failures.
5. Pairwise rank comparison at r=100, actual relative separation >1%, with both
   rows above the floor. Do not call near ties successes. Report molecule-family
   counts; geometries, times, steps and seeds are NOT independent molecules.
6. Compare the dynamic predictor to a static HF local-BCH proxy
   I_static=(T^4/r^2)||Q_0 B_pi HF||^2, also fixed without fitting. This is a
   deliberately simple local-error proxy, NOT a reimplementation of any paper's
   best bound/algorithm. Accurate dynamic prediction alone is not a new theorem.
7. At r=100 all three times, compute exact projected defect decomposition I=D+X,
   C=I/D. C<1 is net destructive accumulation; a reduction in C above 1 is only
   reduced constructive accumulation. Record inside-sector error and leakage.
8. At T=1,r=100 add seven previously defined controls: round robin, within-bucket
   seed20260917, whole-bucket permutations seeds1101..1105: 56 extra target rows.
   Use frozen first-occurrence buckets and identical merged Pauli coefficients.
9. Check bucket Pauli commutation and coefficient-l1 commutator bounds against
   N_alpha/N_beta/N; preserve tiny coefficients in symbolic subtraction.
   Structural and symmetry claims apply only to inputs satisfying these checks.

## Safeguards and limits

All files are additive. No tracked code, previous experiment, or Overleaf edit.
Freeze code/protocol, input hashes, physical ordering hashes, environment, commit.
Closure under individual Pauli flips retains leakage; never project Trotter
evolution into a conserved-number sector. Validate full-vs-reduced one-step
evolution, Hamiltonian action, norm drift, and defect reconstruction.
These four new families give a prospective numerical transfer test, not broad
chemical universality, independent experimental replication, or priority proof.
No optimized ordering, gate-count advantage, fault-tolerant resource advantage,
or scalable classical predictor is claimed. All counterexamples are retained.
