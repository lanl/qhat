general.print_verbose()
general.file_stub = "/tmp/qhat_h2_sto3g/h2_sto3g"
general.file_format = "default"
general.logfile = "/tmp/qhat_h2_sto3g/h2_sto3g_jw.log"

# H2 molecule at equilibrium bond length, in Angstroms.
hamiltonian.add_atom("H", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.7414, 0.0, 0.0)

hamiltonian.basis = "sto-3g"

# H2/STO-3G has 2 occupied spin orbitals and 2 vacant spin orbitals.
# Total active spin orbitals = 4, so n_qubits = 4.
hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 2

hamiltonian.f2q_mapping = "JW"
