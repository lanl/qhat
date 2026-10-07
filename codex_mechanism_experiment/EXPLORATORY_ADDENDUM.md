# Exploratory follow-up declared during the fixed campaign

After the F2 baseline grid and structural controls showed a time-dependent
ranking reversal, two bounded follow-up checks were specified. These are not
part of the original preregistered 129 settings, and are not independent-molecule
validation. Their output is stored in a separate exclusive directory.

1. Check symbolic Pauli commutation within every frozen parent bucket and
   compute coefficient-one-norm bounds on commutators of the bucket Hamiltonian
   with N, Nalpha and Nbeta. Commuting descendants imply exact within-bucket
   permutation invariance. If the summed bucket Hamiltonian also commutes with
   number operators, whole-bucket permutations preserve those symmetries. Check
   this sufficient condition on the actual globally merged coefficients; do
   not assume it from a fermionic label.

2. On F2 only, form the dt^2 coefficient D of S(dt)-exp(-iHdt) as a linear map
   from the exact initial spin sector to the full reachable space. No fitting
   to observed finite-time errors. Propagate the leading global error using
     eta'(t) = -i H eta(t) + D psi_exact(t), eta(0)=0.
   For dt=T/r, predict I ~= dt^2 ||Q_T eta(T)||^2, where Q_T removes the
   component parallel to the exact final state. This is the usual first-order
   perturbative construction, not a newly claimed theorem or efficient selector.
   Compare with the original grid and also predict a new time T=0.75,r=100.
   Save the predictions before calculating the three new T=0.75 Trotter states.
   Record relative residuals, including failures, rather than fitting parameters.

Independent checks: D acting on HF matches the previously validated initial
defect; on a different state in the exact sector, finite-difference defects
approach D linearly with step size. The inhomogeneous evolution is solved by an
augmented sparse exponential. Exact evolution and the numerical product formula
still share the same final nonidentity Pauli Hamiltonian and energy offset.
