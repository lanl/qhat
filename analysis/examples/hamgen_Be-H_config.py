# _________________________________________________________________________________________________
# General configuration

general.print_verbose()                         # Additional information printed out
general.file_stub = "Be-H_1.30_sto-6g"
general.file_format = "default"                 # Use default Pauli string style (not HamLib style)

# _________________________________________________________________________________________________
# Describe the Hamiltonian

hamiltonian.add_atom("Be", 0.00, 0, 0)
hamiltonian.add_atom("H" , 1.30, 0, 0)

hamiltonian.basis = "sto-3g"                    # Select the atomic basis functions

hamiltonian.num_active_occupied = 3             # Specify the active space
hamiltonian.num_active_vacant = 3

