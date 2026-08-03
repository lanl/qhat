# B₂ Random Color-Group Ordering Summary

**Metric:** $1-|\langle\psi_{\mathrm{exact}}|\psi_{\mathrm{Trotter}}\rangle|$ (lower is better).

**Setup:** B2/STO-6G, first-order Trotter evolution, $t=1$, 100 steps, and 100 random color-group permutations per coloring method.

Only rows with `schedule = random_groups` are included in the random-order statistics.

| Qubits | Signed-coefficient baseline | Raw JW | JW mean +/- SD | JW median | JW min-max | JW beating baseline | Fermionic mean +/- SD | Fermionic median | Fermionic min-max | Fermionic beating baseline |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 3.0744e-07 | 2.4663e-07 | 1.2939e-06 +/- 7.2729e-07 | 1.2188e-06 | 1.9641e-07--3.5180e-06 | 9/100 (9.0%) | 9.2625e-07 +/- 6.7464e-07 | 7.7170e-07 | 1.2339e-07--2.9224e-06 | 18/100 (18.0%) |
| 16 | 1.0295e-06 | 4.2453e-06 | 3.8824e-06 +/- 1.7082e-06 | 3.5376e-06 | 1.2133e-06--8.8309e-06 | 0/100 (0.0%) | 2.9332e-06 +/- 1.7139e-06 | 2.3316e-06 | 7.7133e-07--1.0168e-05 | 4/100 (4.0%) |
| 18 | 8.7480e-07 | 2.6310e-06 | 6.5722e-06 +/- 4.7150e-06 | 5.2655e-06 | 1.0954e-06--2.7818e-05 | 0/100 (0.0%) | 4.4265e-06 +/- 2.1227e-06 | 3.9201e-06 | 1.4179e-06--1.2044e-05 | 0/100 (0.0%) |

## Main findings

- The fermionic random-group distribution has a lower mean and median than the JW-coloring distribution for all three active spaces.
- The corrected signed-coefficient/lexicographic baseline is stronger than the typical random coloring schedule in every case.
- Fermionic schedules beating the corrected baseline occur in 18% of the 12-qubit samples, 4% of the 16-qubit samples, and 0% of the 18-qubit samples.
- The 16-qubit fermionic coloring therefore contains a small number of competitive schedules, but random group ordering is not a reliable improvement over the corrected baseline.

The shaded regions in the figure represent mean plus or minus one sample standard deviation.
