import ctypes
import logging
import os
import subprocess
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
# user types
# -------------------------------------------------------------------------------------------------

class GeneralConfigurationUser:
    """
    User-facing configuration for general analysis settings.

    Attributes
    ----------
    logfile : str
        Name of the log file (default: "analysis.log")
    output_directory : str
        Base directory for all output files. If set, all output files (logfile, matrices,
        eigendecompositions, numerical simulation outputs, error analysis) will be written
        to this directory. The directory is created automatically if it doesn't exist.
        If empty or not set, files are written to the current directory.
        Examples:
            "Be-H/" - all outputs go to Be-H/
            "results/run1/" - all outputs go to results/run1/
    """
    def __init__(self):
        # Logfile name
        self.logfile = "analysis.log"
        # How much information to print as the script runs
        self._loglevel = "info"
        # Output directory for all generated files (empty string = current directory)
        self.output_directory = ""
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
            logger.error(f"Already set Hamiltonian source to {self.source}.")
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
            logger.error(f"Already set unitary method to {self.method}.")
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
        self.phase_scale_factor = kwargs.get("phase_scale_factor", 1.01)
        # Validate phase_scale_factor
        if self.phase_scale_factor <= 0:
            raise ValueError(
                f"phase_scale_factor must be positive, got {self.phase_scale_factor}")
        if self.phase_scale_factor < 1.0:
            logger.warning(
                f"phase_scale_factor = {self.phase_scale_factor} < 1.0 will map eigenvalue phases "
                f"outside the range [-π, π], which may cause phase wrapping issues. "
                f"Values >= 1.0 are recommended (default: 1.01).")
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
        self.save_if_present(table, "phase_scale_factor")
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
        # Internal (not user-facing) options ______________________________________________________
        # New flexible output API
        self._operator_output_requests = []
        # External (user-facing) options __________________________________________________________
        # Do resource estimation (e.g., qubit and gate counts)
        # Can be a string ('pyliqtr', 'qualtran', 'cirq') or list (['pyliqtr', 'qualtran'])
        self.resource_estimator = None
        # Write unitary matrix of full algorithm to a file
        self.algorithm_matrix_output_file = None
        # Do numerical simulation with the provided starting state(s)
        # -- numerical_simulation_inputs: using the constructed algorithm
        self.numerical_simulation_inputs = None
        # Memory threshold above which dense matrices are no longer generated
        # - Switch to matrix-free representation, which disables some features such as saving
        #   matrices or eigendecompositions
        self.matrix_memory_threshold_gb = 16.0
        # Compute error based on eigenvalues
        self.enable_eigenvalue_errors = False
        # Compute error based on matrix norm(s)
        self.error_matrix_norms = None
        # Compute error based on reference state(s)
        self.error_state_inputs = None
    # Write an operator to a file
    def save_operator_to_file(self, filename, source, operator_type, energy_shifted, representation):
        """
        Request saving an operator to file.

        All parameters are required to ensure explicit specification.

        Parameters
        ----------
        filename : str
            Output filename (e.g., 'H_exact.npz', 'H_exact_eig.npz')
            Supported extensions: .npz, .h5, .hdf5, .txt (for matrix representation)
        source : {'exact', 'approximate'}
            Which operator to save
            'exact' = true Hamiltonian (no approximations)
            'approximate' = algorithm output (Trotter, LCU, etc.)
        operator_type : {'hamiltonian', 'time_evolution'}
            Which operator form to save
            'hamiltonian' = Hamiltonian matrix H (or its eigenbasis)
            'time_evolution' = Time-evolution operator U = exp(-i*H*t) (or its eigenbasis)
        energy_shifted : bool
            Whether to include energy shift
            False = physical energy scale (can have negative eigenvalues)
            True = QPE energy scale (all eigenvalues positive)
        representation : {'matrix', 'eigendecomposition'}
            How to store the operator
            'matrix' = full matrix representation
            'eigendecomposition' = eigenenergies and eigenvectors

        Output for representation='eigendecomposition' contains:
            eigenenergies : ndarray
                Eigenvalues (sorted ascending)
            eigenvectors : ndarray
                Eigenvectors (columns correspond to eigenvalues)
            metadata : dict
                source, operator_type, energy_shifted, timestamp, dimension

        Examples
        --------
        >>> # Save physical Hamiltonian as a matrix
        >>> analysis.save_operator_to_file(
        ...     filename='H_exact.npz',
        ...     source='exact',
        ...     operator_type='hamiltonian',
        ...     energy_shifted=False,
        ...     representation='matrix'
        ... )

        >>> # Save approximate time-evolution operator (shifted for QPE) as eigendecomposition
        >>> analysis.save_operator_to_file(
        ...     filename='U_approx_eig.npz',
        ...     source='approximate',
        ...     operator_type='time_evolution',
        ...     energy_shifted=True,
        ...     representation='eigendecomposition'
        ... )
        """
        # Validate inputs
        self._validate_output_request(filename, source, operator_type, energy_shifted, representation)

        # Store request
        # TODO: Deduplicate
        self._operator_output_requests.append({
            'filename': filename,
            'source': source,
            'operator_type': operator_type,
            'energy_shifted': energy_shifted,
            'representation': representation
        })

    def _validate_output_request(self, filename, source, operator_type, energy_shifted, representation):
        """Validate output request parameters."""
        valid_sources = ['exact', 'approximate']
        valid_operator_types = ['hamiltonian', 'time_evolution']
        valid_representations = ['matrix', 'eigendecomposition']

        if source not in valid_sources:
            raise ValueError(
                f"source must be one of {valid_sources}, got '{source}'"
            )
        if operator_type not in valid_operator_types:
            raise ValueError(
                f"operator_type must be one of {valid_operator_types}, got '{operator_type}'"
            )
        if not isinstance(energy_shifted, bool):
            raise TypeError(
                f"energy_shifted must be a boolean, got {type(energy_shifted).__name__}"
            )
        if representation not in valid_representations:
            raise ValueError(
                f"representation must be one of {valid_representations}, got '{representation}'"
            )

        # Validate filename
        if not isinstance(filename, str):
            raise TypeError(f"filename must be a string, got {type(filename)}")

        supported_exts = ('.npz', '.h5', '.hdf5', '.txt')
        if not filename.endswith(supported_exts):
            logger.warning(
                f"Filename '{filename}' has unusual extension. "
                f"Supported: {', '.join(supported_exts)}"
            )

    def _generate_TOML_table(self):
        table = tomlkit.table()
        self.save_if_present(table, "resource_estimator")
        self.save_if_present(table, "algorithm_matrix_output_file")
        self.save_if_present(table, "numerical_simulation_inputs")
        self.save_if_present(table, "matrix_memory_threshold_gb")
        self.save_if_present(table, "enable_eigenvalue_errors")
        self.save_if_present(table, "error_matrix_norms")
        self.save_if_present(table, "error_state_inputs")
        if self._operator_output_requests:
            table['operator_output_requests'] = self._operator_output_requests

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
        self.output_directory = user_config.output_directory
        self.git_hash = _get_git_hash()

    def get_output_path(self, filename):
        """
        Get the full output path for a file, respecting output_directory.

        Parameters
        ----------
        filename : str
            The filename or relative path

        Returns
        -------
        str
            Full path with output_directory prepended (if set)

        Notes
        -----
        - If output_directory is empty or None, returns filename unchanged
        - Uses os.path.join() for proper path joining
        - Absolute paths in filename override output_directory
        - Creates parent directories automatically

        Examples
        --------
        >>> config.output_directory = "Be-H/"
        >>> config.get_output_path("analysis.log")
        'Be-H/analysis.log'
        >>> config.get_output_path("logs/debug.log")
        'Be-H/logs/debug.log'
        >>> config.get_output_path("/tmp/file.log")
        '/tmp/file.log'  # absolute path unchanged
        """
        if not self.output_directory:
            output_path = filename
        else:
            output_path = os.path.join(self.output_directory, filename)

        # Create parent directory if it doesn't exist
        parent_dir = os.path.dirname(output_path)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)

        return output_path

    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["logfile"] = self.logfile
        table["loglevel"] = self.loglevel
        if self.output_directory:
            table["output_directory"] = self.output_directory
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
        logger.results(f"results:\n{formatted}")

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
        # Apply output directory if configured
        output_path = self.config_general.get_output_path(filename)
        tomlfile = TOMLFile(output_path)
        logger.info(f"Writing TOML summary file \"{output_path}\".")
        tomlfile.write(document)
        logger.verbose("Summary file saved.")

# -------------------------------------------------------------------------------------------------

