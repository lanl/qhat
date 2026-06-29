# H2 / STO-3G / full active space

general.print_verbose()
general.logfile = "experiments/h2/hamgen_h2_sto3g.log"
general.file_stub = "experiments/h2/h2_0.7414_sto3g"
general.file_format = "default"

# H2 near equilibrium bond length
L = 0.7414

hamiltonian.add_atom("H", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H", 0.0, 0.0, L)

hamiltonian.basis = "sto-3g"

# H2/STO-3G has 2 electrons and 4 spin orbitals.
# Use full active space: 2 occupied spin orbitals + 2 vacant spin orbitals.
hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 2

# For Option 1, start with JW.
hamiltonian.f2q_mapping = "JW"
