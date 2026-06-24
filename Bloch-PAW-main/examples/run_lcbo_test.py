"""LCBO end-to-end test from QRE_test.ipynb, updated for bloch_paw package."""
from __future__ import annotations
import numpy as np
from gpaw import GPAW
from ase.build import bulk
import math
import h5py

from bloch_paw.extractor import PawExtractor
from bloch_paw.reader import PawReader
from bloch_paw.one_norm import OneNormCalculator
from bloch_paw.resources import ResourceEstimator

k_mesh_size = (3, 3, 3)
cell_size = (1, 1, 1)
Nb = 5
h_cut = 0.30
atoms = bulk('H', 'fcc', a=3.67).repeat(cell_size)
calc = GPAW(mode="lcao",
            basis='dzp',
            h=h_cut,
            kpts={'size': k_mesh_size, 'gamma': True},
            xc='PBE',
            nbands=Nb,
            symmetry={'point_group': False, 'time_reversal': False})
atoms.calc = calc
calc.get_potential_energy(atoms)

extractor = PawExtractor(calc, thr_rho=1.0e-3, thr_D=1.0e-2, thr_C=1.0e-2,
                               thr_h=1e-5, thr_kappa=1e-5, nbands=Nb)

hdf5_file = extractor.export_hdf5(filepath="paw_ingredients.h5", write_two_body=False)

hamiltonian = PawReader(hdf5_file)
with hamiltonian:
    inputs = hamiltonian.to_calculator_inputs(lazy=True)
    one_norm_calc = OneNormCalculator(**inputs, thr_rank=3e-5,
                                            sv_floor=1e-12, scale_floor=1e-12)
    lam = one_norm_calc.lambda_one_norm()
    avg_rank_l, avg_rank_zero = one_norm_calc.compute_average_rank()

est = ResourceEstimator.from_hdf5("paw_ingredients.h5")

Rl  = avg_rank_l
R0  = avg_rank_zero

t_per_query = est.toffoli_count_per_be(Rl=Rl, R0=R0)

eps_chem = 1.6e-3    # chemical accuracy (Hartree)
eps_qpe = eps_chem / 5  # QPE share of error budget
q_total = est.total_qubits(Rl=Rl, R0=R0, lam=lam, eps_qpe=eps_qpe)

print(f"Toffoli count per query: {t_per_query}")
print(f"Total qubit count:       {q_total}")
print(f"One-norm lambda:         {lam}")
print()
print("Reference (notebook):")
print("  Toffoli count per query: 206251")
print("  Total qubit count:       103384")
print("  One-norm lambda:         319.4945806013492")
