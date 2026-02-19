from qre_types import GeneralConfigurationUser, \
                      HamiltonianConfiguration, \
                      UnitaryConfiguration, \
                      QPEConfiguration, \
                      AnalysisConfiguration, \
                      State

import argparse

# -------------------------------------------------------------------------------------------------

def load_configuration() -> State:

    # Set up and read command-line arguments
    parser = argparse.ArgumentParser()
    default_config = "config.py"
    parser.add_argument(
            "configuration_file",
            nargs='?',
            default=default_config,
            help=f"Name of the configuration file; defaults to \"{default_config}\"")
    args = parser.parse_args()

    # Read the configuration file
    with open(args.configuration_file, 'r') as fin:
        config_script = fin.read()
    if len(config_script) > 0 and config_script[-1] == '\n':
        config_script = config_script[:-1]

    # Execute the configuration file
    general = GeneralConfigurationUser()
    hamiltonian = HamiltonianConfiguration()
    unitary = UnitaryConfiguration()
    circuit = QPEConfiguration()
    analysis = AnalysisConfiguration()
    def meV_to_Hartree(meV):
        return 3.67493221757e-5 * meV
    exec(config_script)

    # Build the state (does some post-processing of user configuration)
    state = State(config_script, general, hamiltonian, unitary, circuit, analysis)

    state.config_general.log("\n".join([
        f"Contents of configuration file \"{args.configuration_file}\":",
        config_script
        ]))

    return state

