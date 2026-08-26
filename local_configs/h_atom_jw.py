# local_configs/h_atom_jw.py

general.print_verbose()
general.file_stub = "/tmp/qhat_h_atom_smoke/h_atom"
general.file_format = "default"
general.logfile = "/tmp/qhat_h_atom_smoke/h_atom_jw.log"

# Single hydrogen atom at the origin.
hamiltonian.add_atom("H", 0.0, 0.0, 0.0)

hamiltonian.basis = "sto-3g"

# H/STO-3G has 1 electron and 2 spin orbitals.
hamiltonian.num_active_occupied = 1
hamiltonian.num_active_vacant = 1

# Generate the Jordan-Wigner qubit Hamiltonian.
hamiltonian.f2q_mapping = "JW"
