# LiH / STO-3G / reduced active space

general.print_verbose()
general.logfile = "experiments/lih/hamgen_lih_sto3g_as_2_6.log"
general.file_stub = "experiments/lih/lih_1.45_sto3g"
general.file_format = "default"

# LiH bond length, common test value
L = 1.45

hamiltonian.add_atom("Li", 0.0, 0.0, 0.0)
hamiltonian.add_atom("H",  0.0, 0.0, L)

hamiltonian.basis = "sto-3g"

# LiH has 4 electrons total.
# This active space freezes the deepest occupied spin orbitals
# and keeps 2 active occupied + 6 active vacant spin orbitals.
hamiltonian.num_active_occupied = 2
hamiltonian.num_active_vacant = 6

hamiltonian.f2q_mapping = "JW"
