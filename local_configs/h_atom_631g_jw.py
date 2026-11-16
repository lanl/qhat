general.print_verbose()
general.file_stub = "/tmp/qhat_h_atom_631g/h_atom_631g"
general.file_format = "default"
general.logfile = "/tmp/qhat_h_atom_631g/h_atom_631g_jw.log"

hamiltonian.add_atom("H", 0.0, 0.0, 0.0)

# Larger than STO-3G; should provide more than one spatial orbital for H.
hamiltonian.basis = "6-31g"

# Hydrogen has one electron.
# For two spatial orbitals, there are four spin orbitals:
# one active occupied spin orbital and three active vacant spin orbitals.
hamiltonian.num_active_occupied = 1
hamiltonian.num_active_vacant = 3

hamiltonian.f2q_mapping = "JW"
