"""Compare ground-state energies from two QHAT tensor archives."""

from pathlib import Path

import numpy as np
from openfermion import (
    InteractionOperator,
    get_sparse_operator,
    jw_get_ground_state_at_particle_number,
)


# ---------------------------------------------------------------------------
# Experiment settings
# ---------------------------------------------------------------------------

MOLECULE_NAME = "H2"

DEFAULT_THRESHOLD_FILE = Path(
    "diatomic_hydrogen_as-002-002.tensors.npz"
)

ZERO_THRESHOLD_FILE = Path(
    "diatomic_hydrogen_0_as-002-002.tensors.npz"
)

# Number of active electrons.
N_ELECTRONS = 2


# ---------------------------------------------------------------------------
# Load QHAT tensors
# ---------------------------------------------------------------------------

def load_interaction_operator(npz_path):
    """Load a QHAT .tensors.npz file as an InteractionOperator."""
    npz_path = Path(npz_path)

    if not npz_path.is_file():
        raise FileNotFoundError(
            f"Tensor file not found: {npz_path}"
        )

    with np.load(npz_path) as tensors:
        required_arrays = {
            "constant",
            "one_body",
            "two_body",
        }

        available_arrays = set(tensors.files)
        missing_arrays = required_arrays - available_arrays

        if missing_arrays:
            raise KeyError(
                f"{npz_path} is missing arrays: "
                f"{sorted(missing_arrays)}"
            )

        constant = np.asarray(
            tensors["constant"]
        ).item()

        one_body = np.array(
            tensors["one_body"],
            copy=True,
        )

        two_body = np.array(
            tensors["two_body"],
            copy=True,
        )

    if one_body.ndim != 2:
        raise ValueError(
            f"one_body must be a matrix, but its shape is "
            f"{one_body.shape}"
        )

    if one_body.shape[0] != one_body.shape[1]:
        raise ValueError(
            f"one_body must be square, but its shape is "
            f"{one_body.shape}"
        )

    n_spin_orbitals = one_body.shape[0]

    expected_two_body_shape = (
        n_spin_orbitals,
        n_spin_orbitals,
        n_spin_orbitals,
        n_spin_orbitals,
    )

    if two_body.shape != expected_two_body_shape:
        raise ValueError(
            f"two_body has shape {two_body.shape}; expected "
            f"{expected_two_body_shape}"
        )

    if N_ELECTRONS > n_spin_orbitals:
        raise ValueError(
            f"N_ELECTRONS={N_ELECTRONS} exceeds the number of "
            f"spin orbitals, {n_spin_orbitals}"
        )

    # The saved two_body array is already the final
    # InteractionOperator.two_body_tensor.
    #
    # Do not multiply it by 0.5 again.
    interaction_operator = InteractionOperator(
        constant=constant,
        one_body_tensor=one_body,
        two_body_tensor=two_body,
    )

    return interaction_operator, n_spin_orbitals


# ---------------------------------------------------------------------------
# Ground-state calculation
# ---------------------------------------------------------------------------

def calculate_ground_state(
    interaction_operator,
    n_spin_orbitals,
    n_electrons,
):
    """Calculate the ground state in a fixed-electron-number sector."""
    sparse_hamiltonian = get_sparse_operator(
        interaction_operator,
        n_qubits=n_spin_orbitals,
    )

    energy, state = jw_get_ground_state_at_particle_number(
        sparse_hamiltonian,
        n_electrons,
    )

    energy = float(
        np.real_if_close(energy)
    )

    return energy, state


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def main():
    default_operator, default_n_spin_orbitals = (
        load_interaction_operator(
            DEFAULT_THRESHOLD_FILE
        )
    )

    zero_operator, zero_n_spin_orbitals = (
        load_interaction_operator(
            ZERO_THRESHOLD_FILE
        )
    )

    if default_n_spin_orbitals != zero_n_spin_orbitals:
        raise ValueError(
            "The tensor files have different active-space sizes: "
            f"{default_n_spin_orbitals} and "
            f"{zero_n_spin_orbitals} spin orbitals"
        )

    default_energy, default_state = calculate_ground_state(
        default_operator,
        default_n_spin_orbitals,
        N_ELECTRONS,
    )

    zero_energy, zero_state = calculate_ground_state(
        zero_operator,
        zero_n_spin_orbitals,
        N_ELECTRONS,
    )

    signed_energy_difference = (
        zero_energy - default_energy
    )

    absolute_energy_difference = abs(
        signed_energy_difference
    )

    state_overlap = abs(
        np.vdot(
            default_state,
            zero_state,
        )
    )

    state_fidelity = state_overlap**2

    one_body_difference = np.linalg.norm(
        zero_operator.one_body_tensor
        - default_operator.one_body_tensor
    )

    two_body_difference = np.linalg.norm(
        zero_operator.two_body_tensor
        - default_operator.two_body_tensor
    )

    print("=" * 72)
    print(
        f"{MOLECULE_NAME} active-space ground-state "
        "energy comparison"
    )
    print("=" * 72)

    print(
        f"Number of spin orbitals : "
        f"{default_n_spin_orbitals}"
    )

    print(
        f"Number of electrons     : "
        f"{N_ELECTRONS}"
    )

    print()

    print("Default threshold tensor file:")
    print(f"  {DEFAULT_THRESHOLD_FILE}")
    print(
        f"  Ground-state energy = "
        f"{default_energy:.16f} Hartree"
    )

    print()

    print("Zero-threshold tensor file:")
    print(f"  {ZERO_THRESHOLD_FILE}")
    print(
        f"  Ground-state energy = "
        f"{zero_energy:.16f} Hartree"
    )

    print()

    print("Energy comparison:")
    print(
        "  E(zero threshold) - E(default threshold)"
        f" = {signed_energy_difference:+.16e} Hartree"
    )

    print(
        f"  Absolute difference = "
        f"{absolute_energy_difference:.16e} Hartree"
    )

    print(
        f"  Absolute difference = "
        f"{absolute_energy_difference * 1.0e3:.16e} "
        "mHartree"
    )

    print(
        f"  Absolute difference = "
        f"{absolute_energy_difference * 1.0e6:.16e} "
        "microHartree"
    )

    print()

    print("Ground-state comparison:")
    print(
        f"  State overlap  = "
        f"{state_overlap:.16f}"
    )

    print(
        f"  State fidelity = "
        f"{state_fidelity:.16f}"
    )

    print()

    print("Tensor differences:")
    print(
        f"  ||delta one_body||_F = "
        f"{one_body_difference:.16e}"
    )

    print(
        f"  ||delta two_body||_F = "
        f"{two_body_difference:.16e}"
    )


if __name__ == "__main__":
    main()