# Disable qsharp's telemetry.  qsharp checks the environment variable
# QSHARP_PYTHON_TELEMETRY on import, so this needs to appear before
# that import occurs.
import os

os.environ["QSHARP_PYTHON_TELEMETRY"] = "none"

import logging
import math

from qhat.common.logging_utils import configure_logging
from qhat.analysis.algorithm import build_algorithm, compute_initial_phase_qubits
from qhat.analysis.analysis import analyze_algorithm
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

    # Hamiltonian _________________________________________________________________________________

    physical_hamiltonian = get_physical_hamiltonian(state.config_hamiltonian)

    # Compute Trotterization parameters ___________________________________________________________

    tevol_hbar = None
    P0 = None

    if state.config_unitary.method == "ramped trotter":

        # first-pass computation of energy bounds
        Elo1, Ehi1 = physical_hamiltonian.compute_initial_energy_bounds(state.config_hamiltonian)

        # energy-shift Hamiltonian to center at zero
        # This maps eigenvalues from [Elo1, Ehi1] to [-(Ehi1-Elo1)/2, +(Ehi1-Elo1)/2]
        # which enables phase angles in [-π, +π] (matching np.angle output directly)
        E0 = (Elo1 + Ehi1) / 2
        physical_hamiltonian.energy_shift(-E0)
        Elo2 = Elo1 - E0
        Ehi2 = Ehi1 - E0
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

    # Extract timestep for eigendecomposition analysis (if available)
    # For ramped trotter, use tevol_hbar; otherwise check config_unitary.timestep
    timestep = tevol_hbar if tevol_hbar is not None else getattr(state.config_unitary, 'timestep', None)

    # Extract energy shift for correcting eigenvalue comparisons
    energy_shift = physical_hamiltonian.get_energy_shift()

    state.store_results(analyze_algorithm(
        state.config_analysis,
        algorithm,
        state.config_unitary.method,
        approximate_time_evolution=unitary_hamiltonian,
        exact_hamiltonian=physical_hamiltonian,
        timestep=timestep,
        energy_shift=energy_shift
    ))

    # Save Results ________________________________________________________________________________

    state.show_results()
    state.save_summary()

# =================================================================================================

if __name__ == "__main__":
    run()
