# Disable qsharp's telemetry.  qsharp checks the environment variable
# QSHARP_PYTHON_TELEMETRY on import, so this needs to appear before
# that import occurs.
import os

os.environ["QSHARP_PYTHON_TELEMETRY"] = "none"

import logging
import math

from qhat.common.logging_utils import configure_logging
from qhat.analysis.algorithm import build_algorithm, compute_initial_phase_qubits
from qhat.analysis.analysis import analyze_algorithm, validate_and_autocomplete_analysis_config
from qhat.analysis.configuration import load_configuration
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.unitary import encode_as_unitary

logger = logging.getLogger(__name__)

# =================================================================================================

def run():

    # Configuration _______________________________________________________________________________

    state = load_configuration()

    # Configure logging based on user settings
    configure_logging(
        level=state.config_general.loglevel,
        logfile=state.config_general.logfile
    )

    logger.info("=" * 99)
    logger.info("ANALYSIS DRIVER START")
    logger.info("=" * 99)
    logger.info(f"Logfile: {state.config_general.logfile}")
    logger.info(f"Git hash: {state.config_general.git_hash}")

    # Validate analysis configuration _____________________________________________________________

    logger.info("Validating analysis configuration...")
    validate_and_autocomplete_analysis_config(state.config_analysis)
    logger.info("Configuration validated successfully.")

    # Hamiltonian _________________________________________________________________________________

    physical_hamiltonian = get_physical_hamiltonian(state.config_hamiltonian)

    # Compute Trotterization parameters ___________________________________________________________

    tevol_hbar = None
    P0 = None

    if state.config_unitary.method == "ramped trotter":

        # first-pass computation of energy bounds
        Elo1, Ehi1 = physical_hamiltonian.compute_initial_energy_bounds(state.config_hamiltonian)

        # energy-shift Hamiltonian
        physical_hamiltonian.energy_shift(-1 * Elo1)
        Elo2 = Elo1 - Elo1
        Ehi2 = Ehi1 - Elo1
        logger.verbose(f"-- shifted bounds = [{Elo2}, {Ehi2})")
        tevol_hbar = 2 * math.pi / (Ehi2 - Elo2)
        logger.verbose(f"-- preliminary evolution time = {tevol_hbar} * hbar")

        # preliminiary number of phase qubits, with upper bound correction
        P0, Elo3, Ehi3 = compute_initial_phase_qubits(state.config_algorithm, Elo2, Ehi2)
        tevol_hbar = 2 * math.pi / (Ehi3 - Elo3)
        logger.verbose(f"-- optimized evolution time = {tevol_hbar} * hbar")

    # Unitary _____________________________________________________________________________________

    unitary_hamiltonian = encode_as_unitary(
            state.config_unitary,
            physical_hamiltonian,
            tevol_hbar)

    # Algorithm ___________________________________________________________________________________

    algorithm = build_algorithm(
            state.config_algorithm,
            unitary_hamiltonian,
            P0)

    # Analysis ____________________________________________________________________________________

    state.store_results(analyze_algorithm(state.config_analysis, algorithm, hamiltonian=physical_hamiltonian))

    # Save Results ________________________________________________________________________________

    state.show_results()
    state.save_summary()

# =================================================================================================

if __name__ == "__main__":
    run()
