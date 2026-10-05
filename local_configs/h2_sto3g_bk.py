general.print_verbose()
general.file_stub = "/tmp/qhat_h2_sto3g/h2_sto3g"
general.file_format = "default"
general.logfile = "/tmp/qhat_h2_sto3g/h2_sto3g_bk.log"

# H2 molecule at equilibrium bond length, in Angstroms.
hamiltonian.add_atom("H", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.7414, 0.0, 0.0)

hamiltonian.basis = "sto-3g"

hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 2

hamiltonian.f2q_mapping = "BK"
