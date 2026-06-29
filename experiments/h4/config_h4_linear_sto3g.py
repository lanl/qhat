# H4 linear chain / STO-3G / full active space

general.print_verbose()
general.logfile = "experiments/h4/hamgen_h4_linear_sto3g.log"
general.file_stub = "experiments/h4/h4_linear_0.7414_sto3g"
general.file_format = "default"

# Linear H4 chain with nearest-neighbor spacing L
L = 0.7414

hamiltonian.add_atom("H", 0.0, 0.0, 0.0 * L)
hamiltonian.add_atom("H", 0.0, 0.0, 1.0 * L)
hamiltonian.add_atom("H", 0.0, 0.0, 2.0 * L)
hamiltonian.add_atom("H", 0.0, 0.0, 3.0 * L)

hamiltonian.basis = "sto-3g"

# H4/STO-3G has 4 electrons and 8 spin orbitals.
# Full active space: 4 occupied spin orbitals + 4 vacant spin orbitals.
hamiltonian.num_active_occupied = 4
hamiltonian.num_active_vacant = 4

# Mapping is only needed if full hamgen is run.
# The ham2-only script ignores this and stops before mapping.
hamiltonian.f2q_mapping = "JW"
