general.print_verbose()
general.file_stub = "/tmp/qhat_h_atom_631g/h_atom_631g"
general.file_format = "default"
general.logfile = "/tmp/qhat_h_atom_631g/h_atom_631g_bk.log"

hamiltonian.add_atom("H", 0.0, 0.0, 0.0)

hamiltonian.basis = "6-31g"

hamiltonian.num_active_occupied = 1
hamiltonian.num_active_vacant = 3

hamiltonian.f2q_mapping = "BK"
