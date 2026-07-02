# Test configuration for Phase 1: PySCF symmetry detection
# Using Li2 molecule as a simple test case

general.file_stub = "Li2_test_symmetry"

# Li2 geometry (bond length 2.673 Angstroms from CCCBDB)
hamiltonian.add_atom("Li", 0.0, 0.0, -1.3365)
hamiltonian.add_atom("Li", 0.0, 0.0,  1.3365)
hamiltonian.basis = "sto-3g"

# Active space: Li2 has 6 electrons (3 occupied spatial orbitals)
# With STO-3G: 10 spatial orbitals total, 3 occupied, 7 vacant
# Need even number for frozen occupied: freeze 2, keep 1 active
# This gives: 1 active occupied, 2 active vacant
hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 2

# Enable symmetry detection (Phase 1 test)
hamiltonian.enable_symmetry = True

# Use Jordan-Wigner mapping
hamiltonian.f2q_mapping = "Jordan-Wigner"
