import ctypes
import logging
import mendeleev
import os
import pprint
import subprocess
import sys

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

# from https://stackoverflow.com/a/35804945/1791919
def _addLoggingLevel(levelName, levelNum, methodName=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `levelName` becomes an attribute of the `logging` module with the value
    `levelNum`. `methodName` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `methodName` is not specified, `levelName.lower()` is
    used.

    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present

    Example
    -------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not methodName:
        methodName = levelName.lower()
    if hasattr(logging, levelName):
       raise AttributeError('{} already defined in logging module'.format(levelName))
    if hasattr(logging, methodName):
       raise AttributeError('{} already defined in logging module'.format(methodName))
    if hasattr(logging.getLoggerClass(), methodName):
       raise AttributeError('{} already defined in logger class'.format(methodName))
    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(levelNum):
            self._log(levelNum, message, args, **kwargs)
    def logToRoot(message, *args, **kwargs):
        logging.log(levelNum, message, *args, **kwargs)
    logging.addLevelName(levelNum, levelName)
    setattr(logging, levelName, levelNum)
    setattr(logging.getLoggerClass(), methodName, logForLevel)
    setattr(logging, methodName, logToRoot)

def get_log_level(level_str):
    if level_str.lower() == "info":
        return logging.INFO
    elif level_str.lower() == "verbose":
        return logging.VERBOSE
    elif level_str.lower() == "debug":
        return logging.DEBUG
    else:
        raise ValueError(f"Invalid log level requested: \"{level_str}\"")

def _configure_log(user_logfile, log_level):
    # Build the logger
    logger = logging.getLogger()
    logger.setLevel(log_level)
    # Define log format
    formatter = logging.Formatter('{asctime:23} {levelname:>7s} | {message}', style='{')
    # Configure the logger for stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(log_level)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)
    # Configure the logger for the log file (if the user provides a filename)
    file_handler = logging.FileHandler(user_logfile)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Return the logger
    return logger

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
        _addLoggingLevel("VERBOSE", (logging.INFO + logging.DEBUG) // 2)
        self.logger = _configure_log(user_config.logfile, get_log_level(user_config._loglevel))
        self.git_hash = _get_git_hash()
        self.log(f"Running script with git hash {self.git_hash}")
        assert user_config.file_stub is not None
        self.file_stub = user_config.file_stub
        if user_config.file_format in [ "HamLib", "hamlib" ]:
            self.file_format = "hamlib"
        else:
            self.file_format = "default"
        self.log(f"Writing to file stub \"{self.file_stub}\" in format \"{self.file_format}\".")
    def log(self, *args, **kwargs):
        self.logger.info(*args, **kwargs)
    def log_verbose(self, *args, **kwargs):
        self.logger.verbose(*args, **kwargs)
    def log_debug(self, *args, **kwargs):
        self.logger.debug(*args, **kwargs)
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
        self.config_general.log(*args, **kwargs)
    def log_verbose(self, *args, **kwargs):
        self.config_general.log_verbose(*args, **kwargs)
    def log_debug(self, *args, **kwargs):
        self.config_general.log_debug(*args, **kwargs)
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

