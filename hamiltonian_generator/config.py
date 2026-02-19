general.print_verbose()
general.file_stub = "diatomic_lithium"
general.file_format = "default"

L = 2.0
for i in range(2):
            hamiltonian.add_atom("Li", i * L, 0, 0)

hamiltonian.basis = "sto-3g"

hamiltonian.num_active_occupied = 4
hamiltonian.num_active_vacant = 6

