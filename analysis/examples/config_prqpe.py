"""
Partially-Randomized QPE Resource Estimation Example (analytic "prqpe" estimator)

This configuration runs the analytic prqpe resource estimator end-to-end on a real
molecule. Rather than constructing an explicit circuit and counting gates, prqpe
computes Toffoli-based fault-tolerant costs (Toffoli counts, logical qubit counts,
per-circuit Toffoli budgets, and number of circuits) directly from the Pauli
decomposition of the Hamiltonian using closed-form formulas for partially-randomized
QPE (Gunther et al., arXiv:2503.05647).

Because the cost is analytic, the unitary encoding is a no-op (encode_none): the
Hamiltonian is passed through unchanged to the estimator.

This example uses the Be-H second-quantization tensors that ship alongside it
(Be-H_1.30_sto-6g_as-003-003.tensors.npz, a 3-occupied/3-vacant active space).

Usage (run from this examples directory so the relative .npz path resolves):
    cd qhat/analysis/examples
    python -m qhat.analysis.driver config_prqpe.py

It writes a TOML results summary (named by a content hash) and a log file into the
current working directory; look for the [resource_estimates] table.
"""

# =================================================================================================
# GENERAL CONFIGURATION
# =================================================================================================

general.print_verbose()  # Options: print_default(), print_verbose(), print_debug()
general.logfile = "prqpe_analysis.log"

# =================================================================================================
# HAMILTONIAN SPECIFICATION
# =================================================================================================
# Load a second-quantized Hamiltonian from file and map to qubits.
# Supports: .npz (NumPy), .h5/.hdf5 (HDF5). The path is relative to the directory you
# run the driver from (see the Usage note above).

hamiltonian.load_second_quantization(
    "Be-H_1.30_sto-6g_as-003-003.tensors.npz",
    fermion_to_qubit_transform="JW"  # Jordan-Wigner
)

# To generate tensors for your own molecule, use qhat.hamiltonian_generator
# (see qhat/hamiltonian_generator/README.md), or build them with openfermionpyscf as
# shown in the companion script run_prqpe_molecule.py.
#
# Alternative: load a Pauli-string Hamiltonian directly
# hamiltonian.load_pauli_strings("path/to/your_pauli_strings.json")

# =================================================================================================
# UNITARY ENCODING
# =================================================================================================
# The analytic prqpe estimator does not need an explicit circuit encoding, so the
# unitary step is a no-op: the Hamiltonian is passed through unchanged.

unitary.encode_none()

# =================================================================================================
# ALGORITHM SELECTION
# =================================================================================================
# Partially-randomized QPE.

algorithm.method = "QPE: partially randomized"
algorithm.overlap = 1.0  # Overlap of the initial (guiding) state with the target eigenstate

# Optional partially-randomized QPE knobs:
# algorithm.xi = None                  # depth-reduction factor (auto if None)
# algorithm.randomizer = "rte"         # randomizer strategy ("rte" or "qdrift")
# algorithm.commuting_group_size = None

# =================================================================================================
# ANALYSIS CONFIGURATION
# =================================================================================================
# Request the analytic prqpe resource estimate.

analysis.resource_estimator = "prqpe"
analysis.prqpe_target_precision = 0.0016  # Target phase-estimation precision (Hartree)

# C_gs (ground-state Trotter constant):
#   - Omit to auto-estimate on small systems (diagonalization-based, <= 14 qubits).
#   - Set explicitly for large systems where auto-estimation is impractical:
# analysis.prqpe_C_gs = 1.2e-3
#   - Or extrapolate from a fitted (A, b): C_gs = A * lambda**b
# analysis.prqpe_cgs_rule = (A, b)

# Expected output (Toffoli-based costs, written to the [resource_estimates] TOML table):
#   - toffoli_count: total Toffoli gate count
#   - logical_qubits: number of logical qubits required
#   - max_toffoli_per_circuit: maximum Toffoli count per circuit
#   - num_circuits: number of circuits
#   - C_gs / method / metadata
