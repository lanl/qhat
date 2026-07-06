#!/usr/bin/env python
"""
End-to-end PRQPE resource estimate for a molecule, the whole way through QHAT.

This single, self-contained script:
  1. builds a molecule's second-quantized Hamiltonian with openfermionpyscf,
  2. writes the integral tensors in the .npz layout QHAT's loader expects,
  3. writes a QHAT analysis config that selects the analytic "prqpe" estimator, and
  4. runs the QHAT analysis driver on it -- the exact same entry point as
     `python -m qhat.analysis.driver <config>` -- then prints the resource estimates.

So it exercises the real, config-driven QHAT pipeline end-to-end (load second
quantization -> Jordan-Wigner -> no-op "analytic" encoding -> partially-randomized
QPE -> analytic prqpe resource estimate -> TOML summary).

Run it with the project venv (qhat + openfermionpyscf installed):
    python run_prqpe_molecule.py

Edit MOLECULE below to estimate a different system. Small systems (<= 14 qubits)
auto-estimate the ground-state Trotter constant C_gs by diagonalization; for larger
systems set analysis.prqpe_C_gs in the config (see config_prqpe.py).
"""

import os
import subprocess
import sys
import tempfile

try:
    import tomllib  # Python 3.11+ standard library
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib

import numpy as np
from openfermion import MolecularData
from openfermionpyscf import run_pyscf

# --- Molecule definition (edit me) -------------------------------------------------------------
MOLECULE = {
    "name": "H2_sto-3g",
    "geometry": [("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 0.74))],  # Angstrom
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1,  # spin multiplicity, 2S + 1
}
TARGET_PRECISION = 0.0016  # Hartree (chemical accuracy)


def build_tensors_npz(path):
    """Build the molecule and save (constant, one_body, two_body) as a QHAT-format .npz.

    QHAT's load_second_quantization reads ``one_body`` and ``two_body`` (and an
    optional scalar ``constant``) and reconstructs an openfermion InteractionOperator,
    so we save exactly the tensors of ``get_molecular_hamiltonian()``.
    """
    molecule = MolecularData(
        geometry=MOLECULE["geometry"],
        basis=MOLECULE["basis"],
        multiplicity=MOLECULE["multiplicity"],
        charge=MOLECULE["charge"],
    )
    molecule = run_pyscf(molecule, run_scf=True)
    ham = molecule.get_molecular_hamiltonian()  # openfermion InteractionOperator
    np.savez_compressed(
        path,
        constant=ham.constant,
        one_body=ham.one_body_tensor,
        two_body=ham.two_body_tensor,
    )
    return ham.one_body_tensor.shape[0]  # number of spin-orbitals (= qubits under JW)


CONFIG_TEMPLATE = '''\
general.print_verbose()
general.logfile = "prqpe_{name}.log"

hamiltonian.load_second_quantization("{npz}", fermion_to_qubit_transform="JW")

unitary.encode_none()

algorithm.method = "QPE: partially randomized"
algorithm.overlap = 1.0

analysis.resource_estimator = "prqpe"
analysis.prqpe_target_precision = {eps}
'''


def main():
    name = MOLECULE["name"]
    with tempfile.TemporaryDirectory() as workdir:
        npz = os.path.join(workdir, name + ".tensors.npz")
        n_qubits = build_tensors_npz(npz)
        print(f"Built {name}: {n_qubits} spin-orbitals (qubits under Jordan-Wigner).")

        config_path = os.path.join(workdir, "config.py")
        with open(config_path, "w") as fh:
            fh.write(CONFIG_TEMPLATE.format(name=name, npz=npz, eps=TARGET_PRECISION))

        # Run the real QHAT driver, exactly like: python -m qhat.analysis.driver config.py
        print("Running the QHAT analysis driver (prqpe estimator)...")
        proc = subprocess.run(
            [sys.executable, "-m", "qhat.analysis.driver", "config.py"],
            cwd=workdir, capture_output=True, text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + "\n" + proc.stderr + "\n")
            raise SystemExit(f"QHAT driver failed (exit {proc.returncode})")

        tomls = [f for f in os.listdir(workdir) if f.endswith(".toml")]
        if not tomls:
            raise SystemExit("QHAT driver produced no TOML summary")
        with open(os.path.join(workdir, tomls[0]), "rb") as fh:
            summary = tomllib.load(fh)

    est = summary["resource_estimates"]
    print(f"\n=== PRQPE resource estimate for {name} ===")
    print(f"  estimator              : {est['estimator']} ({est['method']})")
    print(f"  toffoli_count          : {est['toffoli_count']:.6g}")
    print(f"  logical_qubits         : {est['logical_qubits']}")
    print(f"  max_toffoli_per_circuit: {est['max_toffoli_per_circuit']:.6g}")
    print(f"  num_circuits           : {est['num_circuits']}")
    print(f"  C_gs                   : {est['C_gs']:.6g}")


if __name__ == "__main__":
    main()
