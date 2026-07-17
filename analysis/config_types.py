import ctypes
import logging
import os
import pprint
import subprocess
import sys
import tomlkit
from tomlkit.toml_file import TOMLFile

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

# TODO: Where should this live?
class ConfigurationBase:
    def save_if_present(self, table, name):
        if hasattr(self, name):
            value = getattr(self, name)
            if value is not None:
                table[name] = value

# -------------------------------------------------------------------------------------------------

# TODO: Where should this live?
def value(x, default):
    if x is None:
        return default
    else:
        return x

# -------------------------------------------------------------------------------------------------
# user types
# -------------------------------------------------------------------------------------------------

class GeneralConfigurationUser:
    def __init__(self):
        # Logfile name
        self.logfile = "analysis.log"
        # How much information to print as the script runs
        self._loglevel = "info"
    def print_default(self):
        self._loglevel= "info"
    def print_verbose(self):
        self._loglevel = "verbose"
    def print_debug(self):
        self._loglevel = "debug"

# -------------------------------------------------------------------------------------------------

class HamiltonianConfiguration(ConfigurationBase):
    def __init__(self):
        self.source = None
        self.lower_bound = float('-inf')
        self.upper_bound = float('inf')
        self.exact_energy_lower_bound = False
        self.exact_energy_upper_bound = False
        self.fermion_to_qubit_transform = None
        self.boson_to_qubit_transform = None
        self.max_bosons_per_state = None
    def _only_once(self):
        if self.source is not None:
            print("Already set Hamiltonian source to {self.source}.")
            assert self.source is None
    def load_second_quantization(self,
                                 filename,
                                 fermion_to_qubit_transform="JW",
                                 boson_to_qubit_transform="binary",
                                 max_bosons_per_state=None):
        self._only_once()
        self.filename = filename
        # Auto-detect file format based on extension
        if filename.endswith('.h5') or filename.endswith('.hdf5'):
            self.source = "hdf5"
        elif filename.endswith('.npy') or filename.endswith('.npz'):
            self.source = "numpy"
        else:
            raise ValueError(f"Unable to determine file format from extension: {filename}. "
                           f"Supported extensions: .h5, .hdf5, .npy, .npz")
        self.fermion_to_qubit_transform = fermion_to_qubit_transform
        self.boson_to_qubit_transform = boson_to_qubit_transform
        self.max_bosons_per_state = max_bosons_per_state
    def load_pauli_strings(self, filename):
        self._only_once()
        self.source = "pauli"
        self.filename = filename
    def set_energy_lower_bound(self, value, exact=False):
        self.lower_bound = value
        self.exact_energy_lower_bound = exact
    def set_energy_upper_bound(self, value, exact=False):
        self.upper_bound = value
        self.exact_energy_upper_bound = exact
    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["source"] = self.source
        self.save_if_present(table, "filename")
        self.save_if_present(table, "fermion_to_qubit_transform")
        self.save_if_present(table, "boson_to_qubit_transform")
        self.save_if_present(table, "max_bosons_per_state")
        table["lower_bound"] = self.lower_bound
        table["exact_energy_lower_bound"] = self.exact_energy_lower_bound
        table["upper_bound"] = self.lower_bound
        table["exact_energy_upper_bound"] = self.exact_energy_upper_bound
        self.save_if_present(table, "format")
        self.save_if_present(table, "filter_metadata")
        return table

# -------------------------------------------------------------------------------------------------

class UnitaryConfiguration(ConfigurationBase):
    def __init__(self):
        self.method = None
    def _only_once(self):
        if self.method is not None:
            print("Already set unitary method to {self.method}.")
            assert self.method is None
    def encode_ramped_trotter(self, **kwargs):
        self._only_once()
        self.method = "ramped trotter"
        self.timestep = kwargs.get("timestep", None)
        # TODO: Pass in a single energy error somewhere, split into Trotter and phase errors.
        self.energy_error = kwargs.get("energy_error", None)
        self.error_scale = kwargs.get("error_scale", 1.0)
        self.trotter_implementation = kwargs.get("trotter_implementation", "flattened")
        self.trotter_combine_terms = kwargs.get("trotter_combine_terms", True)
        self.ordering_method = kwargs.get("ordering_method", None)
        self.trotter_order = kwargs.get("trotter_order", None)
    def encode_pauli_lcu(self, **kwargs):
        self._only_once()
        self.method = "pauli lcu"
        self.energy_error = kwargs["energy_error"]
    def encode_double_factorization(self, **kwargs):
        self._only_once()
        self.method = "double factorization"
        self.energy_error = kwargs["energy_error"]
    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["method"] = self.method
        self.save_if_present(table, "timestep")
        self.save_if_present(table, "energy_error")
        self.save_if_present(table, "error_scale")
        self.save_if_present(table, "trotter_implementation")
        self.save_if_present(table, "trotter_combine_terms")
        self.save_if_present(table, "ordering_method")
        self.save_if_present(table, "trotter_order")
        return table

# -------------------------------------------------------------------------------------------------

class AlgorithmConfiguration(ConfigurationBase):
    def __init__(self):
        self.method = "qualtran textbook"
        self.num_phase_qubits = None
        self.probability_of_failure = None
        self.energy_error = None
    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["method"] = self.method
        self.save_if_present(table, "num_phase_qubits")
        self.save_if_present(table, "probability_of_failure")
        self.save_if_present(table, "energy_error")
        return table

# -------------------------------------------------------------------------------------------------

class AnalysisConfiguration(ConfigurationBase):
    def __init__(self):
        self.resource_estimator = None
        self.matrix_output_file = None
        self.numerical_simulation_inputs = None
        self.exact_matrix_output_file = None
        self.matrix_memory_threshold_gb = 16.0
        self.eigendecomposition_matrices = None
        self.enable_eigenvalue_errors = False
        self.error_matrix_norms = None
        self.error_state_inputs = None
        self.exact_simulation_inputs = None

    def _generate_TOML_table(self):
        table = tomlkit.table()
        self.save_if_present(table, "resource_estimator")
        self.save_if_present(table, "matrix_output_file")
        self.save_if_present(table, "numerical_simulation_inputs")
        self.save_if_present(table, "exact_matrix_output_file")
        self.save_if_present(table, "matrix_memory_threshold_gb")
        self.save_if_present(table, "eigendecomposition_matrices")
        self.save_if_present(table, "enable_eigenvalue_errors")
        self.save_if_present(table, "error_matrix_norms")
        self.save_if_present(table, "error_state_inputs")
        self.save_if_present(table, "exact_simulation_inputs")
        return table

# -------------------------------------------------------------------------------------------------
# internal types and support functions
# -------------------------------------------------------------------------------------------------

def _get_git_hash():
    file_path = os.path.realpath(__file__)
    dirpath = os.path.dirname(file_path)
    commands = ";".join([
            f"pushd {dirpath} > /dev/null",
            "if [[ $(git diff --stat) != '' ]]",
            "then echo $(git rev-parse HEAD)-dirty",
            "else git rev-parse HEAD",
            "fi",
            "popd > /dev/null"
        ])
    output = subprocess.run(commands, shell=True, capture_output=True)
    return output.stdout.decode("utf-8")[:-1]

class GeneralConfiguration:
    def __init__(self, user_config: GeneralConfigurationUser):
        self.logfile = user_config.logfile
        self.loglevel = user_config._loglevel
        self.git_hash = _get_git_hash()

    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["logfile"] = self.logfile
        table["loglevel"] = self.loglevel
        table["git_hash"] = self.git_hash
        return table

# -------------------------------------------------------------------------------------------------

class State:
    def __init__(self,
                 config_script: str,
                 general: GeneralConfigurationUser,
                 hamiltonian: HamiltonianConfiguration,
                 unitary: UnitaryConfiguration,
                 algorithm: AlgorithmConfiguration,
                 analysis: AnalysisConfiguration):
        self.config_script = config_script
        self.config_general = GeneralConfiguration(general)
        self.config_hamiltonian = hamiltonian
        self.config_unitary = unitary
        self.config_algorithm = algorithm
        self.config_analysis = analysis
        self.results = {}
        # Generate hashes
        # -- The hash should exclude the config_script (to avoid changes in the hash from
        #    irrelevant changes in the user configuration file) and config_general (these don't
        #    change the result, simply details of how the script operates).
        # TODO: We were hoping to use the hash to save the Hamiltonian data and the final output.
        #       The idea is that if the user sets the filename, they can overwrite their own data.
        #       If we choose to hash the configuration, then the same Hamiltonian will have the
        #       same hash and thus you can just read the file if it exists or create it if not, and
        #       you'll get the right Hamiltonian without overwriting data.
        #       -- Problem #1: Python's built-in hash() function is non-deterministic (apparently
        #          an intentional choice made for security reasons)
        #       -- Problem #2: If we change anything about the configuration (e.g., add a new
        #          option that was previously defaulted) then the hash changes between versions of
        #          our script.
        #       -- I already added the ability for the user to say, "I put my Hamiltonian in a file
        #          named my own way, so use that file."
        #       -- Maybe we should just let the user set the filename for Hamiltonian data, and let
        #          it be their problem if they overwrite their data?
        #          -- We can ameliorate the problem by one or more of the following:
        #             -- Refuse to overwrite an existing file, raise an exception
        #             -- Put in a flag where the user has to explicitly say, "Yes, I want you to
        #                overwrite the data file."  Some users will just always set that flag to
        #                avoid the crashes, but at least this time it's really really their fault.
        #             -- If the file already exists, have some systematic way to append to the
        #                filename and avoid overwriting, then report the new filename.  Think of
        #                how browsers will often download a file and append " (1)" or similar to
        #                the filename before the extension.
        #       -- Alternately, we need to find a stable hash in Python to address problem #1, and
        #          just accept that problem #2 is going to be a thing.
        self.overall_hash = ctypes.c_size_t(
                hash(
                    tuple(
                        self.__dict__[key] for key in self.__dict__.keys()
                        if key not in ("config_script", "config_general", "results",)
                        )
                    )
                ).value
        self.config_hamiltonian.config_hash = ctypes.c_size_t(hash(self.config_hamiltonian)).value
        logger.info(f"Unique hash for Hamiltonian: {self.config_hamiltonian.config_hash}")
        logger.info(f"Unique hash for full analysis: {self.overall_hash}")
    def store_result(self, key, value):
        self.results[key] = value
    def store_results(self, d):
        self.results.update(d)

    def _format_results(self, obj, indent=0, max_line_length=100):
        """
        Custom formatter for nested dictionaries with readable output.

        Rules:
        - Dictionaries: Each key on own line with nested content indented (4 spaces)
        - Lists/sets that fit on one line: Keep inline
        - Long lists: One element per line, indented
        - Scalar values: Same line as key
        """
        import numpy as np

        indent_str = "    " * indent  # Changed from 3 to 4 spaces

        if isinstance(obj, dict):
            if not obj:
                return "{}"
            lines = []
            for key, value in obj.items():
                # Check if value is a numpy array (treat as list/tuple)
                if isinstance(value, np.ndarray):
                    lines.append(f"{indent_str}{key}:")
                    lines.append(self._format_results(value, indent + 1, max_line_length))
                elif isinstance(value, (dict, list, tuple, set)) and value:
                    # Container types get their own indented section
                    lines.append(f"{indent_str}{key}:")
                    lines.append(self._format_results(value, indent + 1, max_line_length))
                else:
                    # Scalars and empty containers stay on same line
                    lines.append(f"{indent_str}{key}: {value}")
            return "\n".join(lines)

        elif isinstance(obj, (list, tuple)) or isinstance(obj, np.ndarray):
            # Convert numpy arrays to list for consistent handling
            if isinstance(obj, np.ndarray):
                obj_list = obj.tolist()
            else:
                obj_list = obj

            if not obj_list:
                return f"{indent_str}[]"

            # Try formatting on one line first (use list representation, not numpy string)
            one_line = str(obj_list)
            if len(one_line) <= max_line_length:
                return f"{indent_str}{one_line}"

            # Too long - one element per line
            lines = []
            for item in obj_list:
                if isinstance(item, (dict, list, tuple, set)) and item:
                    lines.append(self._format_results(item, indent, max_line_length))
                else:
                    lines.append(f"{indent_str}{item}")
            return "\n".join(lines)

        elif isinstance(obj, set):
            if not obj:
                return f"{indent_str}set()"

            # Try formatting on one line first
            one_line = str(obj)
            if len(one_line) <= max_line_length:
                return f"{indent_str}{one_line}"

            # Too long - one element per line
            lines = []
            for item in sorted(obj):  # Sort for consistency
                lines.append(f"{indent_str}{item}")
            return "\n".join(lines)

        else:
            # Scalar value
            return f"{indent_str}{obj}"

    def show_results(self):
        formatted = self._format_results(self.results, indent=0)
        logger.warning(f"results:\n{formatted}")

    def _filter_for_toml(self, obj):
        """Recursively filter out numpy arrays and other non-serializable objects from results."""
        import numpy as np

        if isinstance(obj, np.ndarray):
            # Skip numpy arrays - they're saved separately in .npz files
            return None
        elif isinstance(obj, dict):
            # Recursively filter dictionary values
            filtered = {}
            for k, v in obj.items():
                filtered_v = self._filter_for_toml(v)
                # Only include if not None (unless original was actually None)
                if filtered_v is not None or (v is None):
                    filtered[k] = filtered_v
            return filtered
        elif isinstance(obj, (list, tuple)):
            # Recursively filter list/tuple elements
            filtered = [self._filter_for_toml(item) for item in obj]
            # Remove None values that came from filtered numpy arrays
            return [item for item in filtered if item is not None]
        else:
            # Return as-is for serializable types
            return obj

    def save_summary(self):
        document = tomlkit.document()
        document.add(tomlkit.comment("CONFIGURATION " + 63 * "="))
        document.add(tomlkit.nl())
        document.add("configuration_file", self.config_script)
        document.add("general", self.config_general._generate_TOML_table())
        document.add("hamiltonian", self.config_hamiltonian._generate_TOML_table())
        document.add("unitary", self.config_unitary._generate_TOML_table())
        document.add("algorithm", self.config_algorithm._generate_TOML_table())
        document.add("analysis", self.config_analysis._generate_TOML_table())
        document.add(tomlkit.nl())
        document.add(tomlkit.comment("RESULTS " + 69 * "="))
        for key in self.results:
            # Filter out numpy arrays and other non-serializable objects
            filtered_value = self._filter_for_toml(self.results[key])
            document.add(key, filtered_value)
        filename = ".".join((str(self.overall_hash), "toml"))
        tomlfile = TOMLFile(filename)
        tomlfile.write(document)
        logger.info(f"Summary file saved to \"{filename}\".")

# -------------------------------------------------------------------------------------------------

