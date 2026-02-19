# Disable qsharp's telemetry.  qsharp checks the environment variable
# QSHARP_PYTHON_TELEMETRY on import, so this needs to appear before
# that import occurs.
import os

os.environ["QSHARP_PYTHON_TELEMETRY"] = "none"

import math

from qre_configuration import load_configuration
from qre_hamiltonian import get_physical_hamiltonian
from qre_unitary import encode_as_unitary
from qre_circuit import build_qpe_circuit, compute_initial_phase_qubits
from qre_analysis import analyze_circuit

# =================================================================================================

def run():

    # Configuration _______________________________________________________________________________

    state = load_configuration()

    # Hamiltonian _________________________________________________________________________________

    physical_hamiltonian = get_physical_hamiltonian(
            state.config_general,
            state.config_hamiltonian)

    # Compute Trotterization parameters ___________________________________________________________

    tevol_hbar = None
    P0 = None

    if state.config_unitary.method == "ramped trotter":

        # first-pass computation of energy bounds
        Elo1, Ehi1 = physical_hamiltonian.compute_initial_energy_bounds(
                state.config_general, state.config_hamiltonian)

        # energy-shift Hamiltonian
        physical_hamiltonian.energy_shift(-1 * Elo1)
        Elo2 = Elo1 - Elo1
        Ehi2 = Ehi1 - Elo1
        state.config_general.log_verbose(f"-- shifted bounds = [{Elo2}, {Ehi2})")
        tevol_hbar = 2 * math.pi / (Ehi2 - Elo2)
        state.config_general.log_verbose(f"-- preliminary evolution time = {tevol_hbar} * hbar")

        # preliminiary number of phase qubits, with upper bound correction
        P0, Elo3, Ehi3 = compute_initial_phase_qubits(
                state.config_general,
                state.config_circuit,
                Elo2, Ehi2)
        tevol_hbar = 2 * math.pi / (Ehi3 - Elo3)
        state.config_general.log_verbose(f"-- optimized evolution time = {tevol_hbar} * hbar")

    # Unitary _____________________________________________________________________________________

    unitary_hamiltonian = encode_as_unitary(
            state.config_general,
            state.config_unitary,
            physical_hamiltonian,
            tevol_hbar)

    # QPE Circuit _________________________________________________________________________________

    qpe_circuit = build_qpe_circuit(
            state.config_general,
            state.config_circuit,
            unitary_hamiltonian,
            P0)

    # Analysis ____________________________________________________________________________________

    state.store_results(analyze_circuit(state.config_general, state.config_analysis, qpe_circuit))

    # Save Results ________________________________________________________________________________

    state.show_results()
    state.save_summary()

# =================================================================================================

if __name__ == "__main__":
    run()
