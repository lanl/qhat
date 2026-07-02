# Benchmark configuration for LiH with SymUCCSD
# From paper Table I: 12 qubits, C2v symmetry, 44→20 parameters (55% reduction)

general.file_stub = "LiH_benchmark"

# LiH geometry from CCCBDB - bond length 1.5949 Angstroms
hamiltonian.add_atom("Li", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.0, 0.0, 1.5949)
hamiltonian.basis = "sto-3g"

# Active space for 12 qubits:
# Li: 3 electrons (1s² 2s¹), H: 1 electron
# Total: 4 electrons = 2 occupied spatial orbitals
# 12 qubits = 6 spatial orbitals → 2 occ + 4 vac
hamiltonian.num_active_occupied = 4  # spin orbitals (2 spatial × 2 spins)
hamiltonian.num_active_vacant = 8    # spin orbitals (4 spatial × 2 spins)

# Enable symmetry detection and reduction
hamiltonian.enable_symmetry = True
hamiltonian.apply_symmetry_reduction = True

hamiltonian.f2q_mapping = "Jordan-Wigner"
