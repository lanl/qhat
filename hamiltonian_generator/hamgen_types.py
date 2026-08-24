import ctypes
import logging
import mendeleev
import os
import pprint
import subprocess
import sys

<<<<<<< HEAD
from qhat.hamiltonian_generator.thresholding import (
    DEFAULT_COEFFICIENT_THRESHOLD,
    validate_coefficient_threshold,
)
=======
from qhat.common.git_utils import get_git_hash
>>>>>>> upstream/main

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------
# user types
# -------------------------------------------------------------------------------------------------

class GeneralConfigurationUser:
    def __init__(self):
        self.file_stub = None
        self.file_format = "default"
        # Logfile name
        self.logfile = "hamgen.log"
        # How much information to print as the script runs
        self._loglevel = "info"
        # Output directory for all generated files (empty string = current directory)
        self.output_directory = ""
        # Cache directory for reusable intermediate files (empty string = use output_directory)
        self.cache_directory = ""
    def print_default(self):
        self._loglevel= "info"
    def print_verbose(self):
        self._loglevel = "verbose"
    def print_debug(self):
        self._loglevel = "debug"

# -------------------------------------------------------------------------------------------------

class HamiltonianConfiguration:
    def __init__(self):
        self._geometry = ()
        self.basis = "sto-3g"
        self.f2q_mapping = "Jordan-Wigner"
        self._coefficient_threshold = DEFAULT_COEFFICIENT_THRESHOLD
        # The number of active occupied and vacant spin (single-occupancy) orbitals in the molecule
        self.num_active_occupied = None
        self.num_active_vacant = None
    @property
    def coefficient_threshold(self):
        return self._coefficient_threshold
    @coefficient_threshold.setter
    def coefficient_threshold(self, value):
        self._coefficient_threshold = validate_coefficient_threshold(value)
    def add_atom(self, element, x, y, z):
        # TODO: Clarify the units involved.  Right now we just pass numbers along, so the user has
        #       to know the right units for the various choices they make.  It would be better to
        #       define the input units and then explicitly convert to the right units for each
        #       package.
        def to_mendeleev_str(element) -> str:
            if isinstance(element, int):
                return mendeleev.element(element).symbol
            elif isinstance(element, str):
                if element.isdigit():
                    return mendeleev.element(int(element)).symbol
                else:
                    return mendeleev.element(element.title()).symbol
        # Tuples are immutable, which makes it harder to work with them.  However, in order to be
        # hashable, we want immutability.  So we do a little extra work to keep the geometry as a
        # tuple through this add_atom() method.
        self._geometry = (*self._geometry,
                          (to_mendeleev_str(element), (float(x), float(y), float(z)))
                         )
    def geometry(self):
        return self._geometry
    # Return the fermion-to-qubit mapping name to a standard name
    def fermion_to_qubit_name(self):
        jw_names = ["Jordan-Wigner", "Jordan Wigner", "JW",
                    "jordan-wigner", "jordan_wigner", "jordan wigner", "jw"]
        bk_names = ["Bravyi-Kitaev", "Bravyi Kitaev", "BK",
                    "bravyi-kitaev", "bravyi_kitaev", "bravyi kitaev", "bk"]
        if self.f2q_mapping in jw_names:
            return "jordan-wigner"
        elif self.f2q_mapping in bk_names:
            return "bravyi-kitaev"
    def as_tag(self):
        # Construct the filename tag for the active space
        # TODO: Make this smarter by padding numbers with zeros based on the largest numbers?
        return f"as-{self.num_active_occupied:03d}-{self.num_active_vacant:03d}"
    def f2q_tag(self):
        if self.fermion_to_qubit_name() == "jordan-wigner":
            return "jw"
        elif self.fermion_to_qubit_name() == "bravyi-kitaev":
            return "bk"

# -------------------------------------------------------------------------------------------------
# internal types and support functions
# -------------------------------------------------------------------------------------------------

class GeneralConfiguration:
    def __init__(self, user_config: GeneralConfigurationUser):
        self.logfile = user_config.logfile
        self.loglevel = user_config._loglevel
        self.output_directory = user_config.output_directory
        self.cache_directory = user_config.cache_directory
        self.git_hash = get_git_hash(reference_file=__file__)
        logger.info(f"Running script with git hash {self.git_hash}")
        assert user_config.file_stub is not None
        self.file_stub = user_config.file_stub
        if user_config.file_format in [ "HamLib", "hamlib" ]:
            self.file_format = "hamlib"
        else:
            self.file_format = "default"
        logger.info(f"Writing to file stub \"{self.file_stub}\" in format \"{self.file_format}\".")

    def get_output_path(self, filename):
        """
        Get the full output path for a file, respecting output_directory.

        Mirrors the implementation in analysis/config_types.py for consistency.

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

    def get_cache_path(self, filename):
        """
        Get the full cache path for a file, respecting cache_directory.

        Used when looking for reusable intermediate files (ham1, ham2).

        Parameters
        ----------
        filename : str
            The filename or relative path

        Returns
        -------
        str
            Full path with cache_directory or output_directory prepended (if set)

        Notes
        -----
        - If cache_directory is set, uses cache_directory
        - Otherwise falls back to output_directory (if set)
        - If both are empty, returns filename unchanged (current directory)
        - Does NOT create directories (cache is read-only)
        """
        if self.cache_directory:
            return os.path.join(self.cache_directory, filename)
        elif self.output_directory:
            return os.path.join(self.output_directory, filename)
        else:
            return filename

    def ham3_ext(self):
        if self.file_format == "hamlib":
            return "dat" # TODO: return "hdf5"
        elif self.file_format == "default":
            return "dat"
        else:
            raise NotImplementedError("ham3_ext should never get here")

# -------------------------------------------------------------------------------------------------

class State:
    def __init__(self,
                 config_script: str,
                 general: GeneralConfigurationUser,
                 hamiltonian: HamiltonianConfiguration):
        self.config_script = config_script
        self.config_general = GeneralConfiguration(general)
        self.config_hamiltonian = hamiltonian
        self.metadata = dict()
    def log(self, *args, **kwargs):
        logger.info(*args, **kwargs)
    def log_verbose(self, *args, **kwargs):
        logger.verbose(*args, **kwargs)
    def log_debug(self, *args, **kwargs):
        logger.debug(*args, **kwargs)
    def filename_ham1(self):
        return "{stub}.pickle".format(stub=self.config_general.file_stub)
    def filename_ham2(self):
        return "{stub}_{astag}.pickle".format(
                stub=self.config_general.file_stub,
                astag=self.config_hamiltonian.as_tag())
    def filename_ham3(self):
        return "{stub}_{astag}_{f2qtag}.{ext}".format(
                stub=self.config_general.file_stub,
                astag=self.config_hamiltonian.as_tag(),
                f2qtag=self.config_hamiltonian.f2q_tag(),
                ext=self.config_general.ham3_ext())

# -------------------------------------------------------------------------------------------------
