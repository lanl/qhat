# Test configuration for BeH2 with SymUCCSD symmetry reduction
# Expected: ~74% reduction (90 → 23 parameters) according to paper

general.file_stub = "BeH2_symtest"

# BeH2 geometry from CCCBDB (linear molecule)
# Be at origin, H atoms at ±1.3264 Angstroms
hamiltonian.add_atom("Be", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.0, 0.0, -1.3264)
hamiltonian.add_atom("H", 0.0, 0.0, 1.3264)
hamiltonian.basis = "sto-3g"

# Active space from paper: 14 qubits
# Be: 4 electrons, 2 core + 2 valence
# Each H: 1 electron
# Total: 6 electrons = 3 occupied spatial orbitals
# With STO-3G: Be has ~5 basis functions, each H has 1
# Active space: 3 occupied, 4 vacant (doubling gives 6 spin-occ, 8 spin-vac = 14 total)
hamiltonian.num_active_occupied = 6  # spin orbitals
hamiltonian.num_active_vacant = 8    # spin orbitals

# Enable symmetry detection and reduction
hamiltonian.enable_symmetry = True
hamiltonian.apply_symmetry_reduction = True

# Jordan-Wigner mapping
hamiltonian.f2q_mapping = "Jordan-Wigner"
