# _________________________________________________________________________________________________
# General configuration

general.print_verbose()                         # Additional information printed out
general.file_stub = "diatomic_lithium"          # Base name that all filenames are built from
general.file_format = "default"                 # Use default Pauli string style (not HamLib style)

# _________________________________________________________________________________________________
# Describe the Hamiltonian

L = 2.0
for i in range(2):
    hamiltonian.add_atom("Li", i * L, 0, 0)     # Add two lithium atoms at (0,0,0) and (L,0,0)

hamiltonian.basis = "sto-3g"                    # Select the atomic basis functions

# Absolute cutoff applied while building spin-orbital tensors; the default is 1e-8.
hamiltonian.coefficient_threshold = 1e-8

hamiltonian.num_active_occupied = 4             # Specify the active space
hamiltonian.num_active_vacant = 6
