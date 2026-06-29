from qhat.common.logging_utils import configure_logging
from qhat.hamiltonian_generator import hamgen


def main():
    state = hamgen.load_configuration()

    configure_logging(
        level=state.config_general.loglevel,
        logfile=state.config_general.logfile,
    )

    ham2 = hamgen.get_ham2(state)

    print("Generated or loaded active-space fermionic Hamiltonian.")
    print(f"ham2 pickle: {state.filename_ham2()}")
    print(f"n_qubits / spin orbitals: {ham2.n_qubits}")

    tensor_file = state.filename_ham2()
    tensor_file = tensor_file[: tensor_file.rfind(".")] + ".tensors.npz"
    print(f"ham2 tensors: {tensor_file}")


if __name__ == "__main__":
    main()
