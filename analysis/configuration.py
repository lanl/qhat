import argparse
import logging

from qhat.analysis.config_types import (
    AnalysisConfiguration,
    AlgorithmConfiguration,
    GeneralConfigurationUser,
    HamiltonianConfiguration,
    State,
    UnitaryConfiguration,
)
from qhat.common.config_utils import (
    parse_key_value_params,
    get_standard_exec_namespace,
)

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def load_configuration() -> tuple[State, str]:

    # Set up and read command-line arguments
    parser = argparse.ArgumentParser()
    default_config = "config.py"
    parser.add_argument(
            "configuration_file",
            nargs='?',
            default=default_config,
            help=f"Name of the configuration file; defaults to \"{default_config}\"")

    # Add support for arbitrary key=value arguments
    parser.add_argument(
            '--param', '-p',
            action='append',
            dest='params',
            default=[],
            metavar='KEY=VALUE',
            help='Parameters to pass to the configuration file (e.g., -p distance=1.5)')

    args = parser.parse_args()

    # Parse the key=value parameters
    config_params = parse_key_value_params(args.params)

    # Read the configuration file
    with open(args.configuration_file, 'r') as fin:
        config_script = fin.read()
    if len(config_script) > 0 and config_script[-1] == '\n':
        config_script = config_script[:-1]

    # Execute the configuration file
    general = GeneralConfigurationUser()
    hamiltonian = HamiltonianConfiguration()
    unitary = UnitaryConfiguration()
    algorithm = AlgorithmConfiguration()
    analysis = AnalysisConfiguration()

    # Create namespace with config objects, params, and standard utilities
    exec_namespace = get_standard_exec_namespace()
    exec_namespace.update({
        # configuration objects
        'general': general,
        'hamiltonian': hamiltonian,
        'unitary': unitary,
        'algorithm': algorithm,
        'analysis': analysis,
        # command-line parameters
        'params': config_params,
    })
    exec(config_script, exec_namespace)

    # Build the state (does some post-processing of user configuration)
    state = State(config_script, general, hamiltonian, unitary, algorithm, analysis)

    # Prepare log messages with configuration file contents and any parameters passed
    # These will be logged after logging is configured
    config_file_message = f"Contents of configuration file \"{args.configuration_file}\":\n{config_script}"

    params_message = None
    if config_params:
        params_lines = ["Command-line parameters:"]
        for key, value in config_params.items():
            params_lines.append(f"  {key} = {value!r}")
        params_message = "\n".join(params_lines)

    return state, config_file_message, params_message

