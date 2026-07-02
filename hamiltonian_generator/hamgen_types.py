import ctypes
import logging
import mendeleev
import os
import pprint
import subprocess
import sys

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
        # The number of active occupied and vacant spin (single-occupancy) orbitals in the molecule
        self.num_active_occupied = None
        self.num_active_vacant = None
        # SymUCCSD symmetry reduction options
        self.enable_symmetry = False  # Enable PySCF symmetry detection
        self.apply_symmetry_reduction = False  # Filter Hamiltonian by point-group symmetry
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
    def sym_tag(self, point_group: str = None):
        """Return symmetry tag for filenames (e.g., '_symD2h')."""
        if self.apply_symmetry_reduction and point_group:
            return f"_sym{point_group}"
        return ""

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
        logger.info(f"Running script with git hash {self.git_hash}")
        assert user_config.file_stub is not None
        self.file_stub = user_config.file_stub
        if user_config.file_format in [ "HamLib", "hamlib" ]:
            self.file_format = "hamlib"
        else:
            self.file_format = "default"
        logger.info(f"Writing to file stub \"{self.file_stub}\" in format \"{self.file_format}\".")
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
        # Symmetry data (set during Stage 1: Hartree-Fock)
        self.point_group = None  # Point group name (e.g., "D2h", "C2v")
        self.mo_irreps = None  # Array of irrep IDs for each molecular orbital
    def log(self, *args, **kwargs):
        logger.info(*args, **kwargs)
    def log_verbose(self, *args, **kwargs):
        logger.verbose(*args, **kwargs)
    def log_debug(self, *args, **kwargs):
        logger.debug(*args, **kwargs)
    def filename_ham1(self):
        return "{stub}.pickle".format(stub=self.config_general.file_stub)
    def filename_ham2(self):
        symtag = self.config_hamiltonian.sym_tag(self.point_group)
        return "{stub}_{astag}{symtag}.pickle".format(
                stub=self.config_general.file_stub,
                astag=self.config_hamiltonian.as_tag(),
                symtag=symtag)
    def filename_ham3(self):
        symtag = self.config_hamiltonian.sym_tag(self.point_group)
        return "{stub}_{astag}{symtag}_{f2qtag}.{ext}".format(
                stub=self.config_general.file_stub,
                astag=self.config_hamiltonian.as_tag(),
                symtag=symtag,
                f2qtag=self.config_hamiltonian.f2q_tag(),
                ext=self.config_general.ham3_ext())

# -------------------------------------------------------------------------------------------------

