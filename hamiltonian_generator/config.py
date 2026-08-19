# _________________________________________________________________________________________________
# General configuration

general.print_verbose()                           # Additional information printed out
general.file_stub = "diatomic_lithium"            # Base name that all filenames are built from
general.file_format = "default"                   # Default Pauli string style (not HamLib style)

# Output directory - keeps all generated files organized
# Override with: -p output_directory=custom_dir
# Defaults to current directory (old behavior)
general.output_directory = params.get("output_directory", "")

# Cache directory - where to look for reusable intermediate files (ham1, ham2)
# If not set, falls back to output_directory
# Override with: -p cache_directory=previous_run
general.cache_directory = params.get("cache_directory", "")

# _________________________________________________________________________________________________
# Describe the Hamiltonian

L = params.get("L", 2.0)
for i in range(2):
    hamiltonian.add_atom("Li", i * L, 0, 0)       # Add two lithium atoms at (0,0,0) and (L,0,0)

hamiltonian.basis = params.get("basis", "sto-3g") # Select the atomic basis functions

hamiltonian.num_active_occupied = 4               # Specify the active space
hamiltonian.num_active_vacant = 6

