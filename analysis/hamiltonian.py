from functools import cache, reduce
import h5py
import json
import logging
import numpy as np
import scipy.constants as sc

from openfermion import (
    InteractionOperator,
    QubitOperator,
    binary_code_transform,
    bravyi_kitaev,
    count_qubits,
    jordan_wigner,
)
from pyscf import ao2mo, gto, lib, scf

from qhat.analysis.config_types import GeneralConfiguration, HamiltonianConfiguration, value
from qhat.common.bosons_binary import BosonicBinaryEncoding
from qhat.common.MixedFermionBosonOperator import MixedFermionBosonOperator
from qhat.common.pauli_string import PauliString

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def boson_to_qubit_operator(bosonic_operator, Nmax):
    qubit_op = QubitOperator()
    qubits_per_mode = int(np.ceil(np.log2(Nmax + 1)))

fermionic_mapping = {
    "JW" : jordan_wigner,
    "BK" : bravyi_kitaev,
}
bosonic_mapping = {
    "binary" : BosonicBinaryEncoding,
    #"unary"  : BosonicUnaryEncoding,
}

# -------------------------------------------------------------------------------------------------

def sparse_to_dense_pauli(sparse_pauli, num_qubits):
    """Compatibility wrapper for converting a sparse Pauli tuple to a string."""
    return PauliString.from_sparse(sparse_pauli, num_qubits).to_dense()


def dense_to_sparse_pauli(dense_pauli):
    """Compatibility wrapper for converting a Pauli string to a sparse tuple."""
    try:
        return PauliString.from_dense(dense_pauli).to_sparse()
    except ValueError as exc:
        raise ValueError(f"Invalid character in dense pauli string: \"{dense_pauli}\".") from exc

# -------------------------------------------------------------------------------------------------

class LinearCombinationOfPauliStrings:
    def __init__(self, **kwargs):
        formats = [fmt for fmt in ["dense", "sparse"] if fmt in kwargs]
        if not formats:
            raise ValueError("No data provided to LinearCombinationOfPauliStrings.")
        if len(formats) > 1:
            raise ValueError("Too many formats provided to LinearCombinationOfPauliStrings.")

        input_format = formats[0]
        input_data = kwargs[input_format]
        if not isinstance(input_data, dict):
            raise TypeError("Pauli string data must be provided as a dictionary.")

        identity = PauliString.from_sparse((), kwargs["num_qubits"])
        self._nq = identity.num_qubits
        self._data = {}
        for raw_pauli, coefficient in input_data.items():
            if input_format == "dense":
                pauli = PauliString.from_dense(raw_pauli)
                if pauli.num_qubits != self._nq:
                    raise ValueError(
                        f"Dense Pauli string {raw_pauli!r} has length "
                        f"{pauli.num_qubits}, expected {self._nq}."
                    )
            else:
                pauli = PauliString.from_sparse(raw_pauli, self._nq)
            self._data[pauli] = self._data.get(pauli, 0.0) + coefficient

    def num_qubits(self):
        return self._nq

    def get_pauli_strings(self):
        """Return a copy of the canonical PauliString-keyed coefficient dictionary."""
        return dict(self._data)

    def get_dense_pauli_strings(self):
        return {pauli.to_dense() : coef for pauli, coef in self._data.items()}

    def get_sparse_pauli_strings(self):
        return {pauli.to_sparse() : coef for pauli, coef in self._data.items()}

    def energy_shift(self, shift):
        identity = PauliString.from_sparse((), self._nq)
        identity_coefficient = self._data.get(identity, 0.0) + shift
        self._data[identity] = identity_coefficient

# -------------------------------------------------------------------------------------------------

# TODO: The heavy use of isinstance() suggests that perhaps Hamiltonian should be a base class that
#       other things are built on top of?

class Hamiltonian:
    def __init__(self, hamiltonian):
        self._H = hamiltonian
    def get_core_operator(self):
        # TODO: This should probably be replaced by a function that generates an appropriate
        #       pyLIQTR problem instance.
        return self._H
    def set_fermionic_mapping(self, mapping):
        # self._fmap is never used for LinearCombinationOfPauliStrings
        if isinstance(mapping, str):
            self._fmap = fermionic_mapping[mapping]
        else:
            self._fmap = mapping
        if isinstance(self._H, MixedFermionBosonOperator):
            self._H.set_fermionic_encoding(self._fmap)
    def set_bosonic_mapping(self, mapping, max_bosons_per_state):
        # self._bmap is never used for LinearCombinationOfPauliStrings
        if isinstance(mapping, str):
            self._bmap = bosonic_mapping[mapping](max_bosons_per_state)
        else:
            self._bmap = mapping
        if isinstance(self._H, MixedFermionBosonOperator):
            self._H.set_bosonic_encoding(self._bmap)
    def num_qubits(self):
        if isinstance(self._H, InteractionOperator):
            return self._H.n_qubits
        elif isinstance(self._H, MixedFermionBosonOperator):
            return self._H.num_qubits()
        elif isinstance(self._H, LinearCombinationOfPauliStrings):
            return self._H.num_qubits()
        else:
            raise TypeError("Unable to determine the number of qubits.")
    def get_all_pauli_strings(self, return_as="tuples"):
        # Returns all Pauli strings as a flat data structure, specifically a dictionary where the
        # key is the Pauli string and the value is the coefficient.
        # -- If return_as == "tuples": The Pauli string is encoded as a tuple of tuples, where each
        #    inner tuple is (qubit index, Pauli operator), with the qubit index being a
        #    zero-indexed integer and the Pauli operator is a one-character string.  For example,
        #    assuming at least 4 qubits, ((0, 'X'), (3, 'Z')).
        # -- If return_as == "strings": The Pauli string is encoded as a character string, where
        #    each character is a Pauli matrix, explicitly including identity entries.  For example,
        #    assuming 6 qubits, "XIIZII".
        # -- If return_as == "objects": Keys are immutable PauliString instances that can provide
        #    either representation on demand.
        # TODO: I'd prefer that the flag identify not the data structure but the concept: dense vs
        #       sparse, rather than strings vs tuples.
        if return_as == "tuples":
            if isinstance(self._H, InteractionOperator):
                return self._fmap(self._H).terms
            elif isinstance(self._H, MixedFermionBosonOperator):
                # TODO: MixedFermionBosonOperator already has its encodings, so this ignores the
                #       encodings selected by Hamiltonian.  Clean this up.  Probably by deferring
                #       the specification of encodings for MixedFermionBosonOperator?
                return self._H.generate_qubit_operator().terms
            elif isinstance(self._H, LinearCombinationOfPauliStrings):
                return self._H.get_sparse_pauli_strings()
            else:
                raise TypeError(
                    f"Unable to generate Pauli strings from object of type \"{type(self._H)}\".")
        elif return_as == "strings":
            as_objects = self.get_all_pauli_strings(return_as="objects")
            return {pauli.to_dense() : coef for pauli, coef in as_objects.items()}
        elif return_as == "objects":
            if isinstance(self._H, LinearCombinationOfPauliStrings):
                return self._H.get_pauli_strings()
            as_tuples = self.get_all_pauli_strings(return_as="tuples")
            Nq = self.num_qubits()
            return {PauliString.from_sparse(pauli, Nq) : coef
                    for pauli, coef in as_tuples.items()}
        else:
            raise ValueError("  ".join([
                "The value of return_as must be \"tuples\", \"strings\", or \"objects\".",
                f"Unable to return result as \"{return_as}\".",
                ]))
    def get_grouped_terms(self):
        # Returns all Pauli strings as a list of QubitOperator instances..
        # TODO: I think a QubitOperator can hold a sum of terms, so presumably this structure would
        #       still work for grouped terms.
        # TODO: Should this be cached?
        # TODO: If get_all_pauli_strings is used to transform the Hamiltonian into Pauli strings,
        #       then it can't be based on get_grouped terms.  Can get_grouped_terms be based on
        #       get_all_pauli_strings?  My concern is that we don't want to carry multiple copies
        #       of all the Pauli strings if we don't need to, but we also don't want to add lots of
        #       indirection.  Instead of using functools.cache, I may have to manually cache and
        #       check for (a) if the list exists, use it; (b) if the grouped structure exists,
        #       present it in a flattened way; (c) otherwise compute.  And then get_grouped_terms
        #       would have to do (a) if the list doesn't exist, call get_all_pauli_strings; (b) if
        #       the list exists, process it into the grouped data structure, save that, delete the
        #       list.  This is the sort of clean-up / optimization stuff that I should put off
        #       until later, because right now I just need to get it working.
        groups = list()
        for pauli, coef in self.get_all_pauli_strings().items():
            groups.append(QubitOperator(pauli, coef))
        return groups
    def to_matrix(self, memory_threshold_gb, force_dense=None, force_sparse=None):
        """
        Convert the Hamiltonian to its exact matrix representation.

        The choice between dense and sparse representation is based on the memory
        threshold. Dense matrices require (2^num_qubits)^2 * 16 bytes.

        Parameters:
            memory_threshold_gb: float, memory threshold in GB for dense representation
            force_dense: bool or None, force dense matrix representation
            force_sparse: bool or None, force matrix-free operator

        Returns:
            numpy array (dense) or PauliStringOperator (sparse/matrix-free)

        Raises:
            ValueError: If both force_dense and force_sparse are True

        Example:
            >>> hamiltonian = Hamiltonian(...)
            >>> H_exact = hamiltonian.to_matrix(memory_threshold_gb=16.0)
            >>> # For small systems: H_exact is numpy array
            >>> # For large systems: H_exact is PauliStringOperator
        """
        from qhat.analysis.matrix_operations import create_hamiltonian_operator

        # Get Pauli strings as dense format (full string per term)
        pauli_dict = self.get_all_pauli_strings(return_as="strings")
        num_qubits = self.num_qubits()

        logger.verbose(f"Converting Hamiltonian to matrix for {num_qubits} qubits")
        logger.verbose(f"Hamiltonian has {len(pauli_dict)} Pauli terms")

        # Create operator (dense or sparse based on memory threshold)
        return create_hamiltonian_operator(
            pauli_dict, num_qubits, memory_threshold_gb,
            force_dense=force_dense,
            force_sparse=force_sparse
        )
    def energy_shift(self, dE):
        if isinstance(self._H, InteractionOperator):
            t0 = self._H.constant + dE
            t1 = self._H.one_body_tensor
            t2 = self._H.two_body_tensor
            self._H = InteractionOperator(t0, t1, t2)
        elif isinstance(self._H, MixedFermionBosonOperator):
            self._H.energy_shift(dE)
        elif isinstance(self._H, LinearCombinationOfPauliStrings):
            self._H.energy_shift(dE)
        else:
            raise TypeError(
                    f"Unable to shift a fermionic Hamiltonian of type \"{type(self._H)}\".")
    def compute_initial_energy_bounds(
            self,
            config_hamiltonian: HamiltonianConfiguration):
        logger.info("Computing initial energy bounds.")
        pauli_sum = self.get_all_pauli_strings()
        logger.verbose(f"-- number of Pauli strings = {len(pauli_sum)}")
        energy_shift = pauli_sum.get(tuple(), 0.0) # identity term (may not exist in all formats)
        dE = sum(abs(coefficient) for coefficient in pauli_sum.values()) - abs(energy_shift)
        Elo0 = energy_shift - dE
        Ehi0 = energy_shift + dE
        logger.verbose(f"-- energy shift = {energy_shift}")
        logger.verbose(f"-- computed bounds = [{Elo0}, {Ehi0})")
        EloU = value(config_hamiltonian.lower_bound, float('-inf'))
        EhiU = value(config_hamiltonian.upper_bound, float('inf'))
        Elo1 = max(Elo0, EloU)
        Ehi1 = min(Ehi0, EhiU)
        logger.verbose(f"-- limited bounds = [{Elo1}, {Ehi1})")
        if config_hamiltonian.exact_energy_lower_bound:
            assert config_hamiltonian.lower_bound is not None
            Elo1 = config_hamiltonian.lower_bound
        if config_hamiltonian.exact_energy_upper_bound:
            assert config_hamiltonian.upper_bound is not None
            Ehi1 = config_hamiltonian.upper_bound
        logger.verbose(f"-- initial bounds = [{Elo1}, {Ehi1})")
        return (Elo1, Ehi1)

# -------------------------------------------------------------------------------------------------

def _verify_and_construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb):
    assert len(f1.shape) == 2
    Nf = f1.shape[0]
    assert f1.shape[1] == Nf
    assert len(f2.shape) == 4
    assert f2.shape[0] == Nf
    assert f2.shape[1] == Nf
    assert f2.shape[2] == Nf
    assert f2.shape[3] == Nf
    Nb = 0
    if fb is not None:
        assert len(fb.shape) == 3
        assert fb.shape[0] == Nf
        assert fb.shape[1] == Nf
        Nb = fb.shape[2]
    if bs is None and fb is None:
        H = Hamiltonian(InteractionOperator(f0, f1, f2))
        logger.info(f"Fermionic second-quantization Hamiltonian uses {H.num_qubits()} qubits.")
        H.set_fermionic_mapping(config_hamiltonian.fermion_to_qubit_transform)
        logger.info(" ".join(["Mapping fermionic operators to qubit operaturs using",
                                     f"{config_hamiltonian.fermion_to_qubit_transform} method."]))
        return H
    else:
        assert bs is not None
        assert fb is not None
        assert Nb > 0
        H = Hamiltonian(MixedFermionBosonOperator(f0, f1, f2, bs, fb))
        logger.info(" ".join([
            "mixed fermionic-bosonic second-quantization Hamiltonian",
            f"uses {H.num_qubits()} qubits."]))
        H.set_fermionic_mapping(config_hamiltonian.fermion_to_qubit_transform)
        logger.info(" ".join(["Mapping fermionic operators to qubit operaturs using",
                                     f"{config_hamiltonian.fermion_to_qubit_transform} method."]))
        if config_hamiltonian.max_bosons_per_state is None:
            raise ValueError("User did not specify maximum bosons per state.")
        H.set_bosonic_mapping(
                config_hamiltonian.boson_to_qubit_transform,
                config_hamiltonian.max_bosons_per_state)
        logger.info(" ".join(["Mapping bosonic operators to qubit operaturs using",
                                     f"{config_hamiltonian.boson_to_qubit_transform} method",
                                     f"with {config_hamiltonian.max_bosons_per_state} maximum",
                                     "bosons per state."]))
        return H

# -------------------------------------------------------------------------------------------------

def load_hdf5(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading second-quantization Hamiltonian from HDF5 file \"{filename}\".")
    data = h5py.File(filename)
    f0 = 0      # Currently don't support constant terms in HDF5
    f1 = data["1e"]
    f2 = data["2e"]
    bs = None   # Currently don't support bosons in HDF5
    fb = None   # Currently don't support bosons in HDF5
    return _verify_and_construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb)

# -------------------------------------------------------------------------------------------------

def load_numpy(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading second-quantization Hamiltonian from file \"{filename}\".")
    data = np.load(filename)
    def get_optional_scalar(name, default_value):
        x = data.get(name, None)
        if x is None:
            return default_value
        else:
            return x[()] # extract scalar from 0D NumPy array
    f0 = get_optional_scalar("constant", 0)
    f1 = data["one_body"]
    f2 = data["two_body"]
    bs = get_optional_scalar("bosonic_scalar", None)
    fb = data.get("fb_interaction", None)
    return _verify_and_construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb)

# -------------------------------------------------------------------------------------------------

def load_hamlib_hdf5(config_hamiltonian: HamiltonianConfiguration):
    """
    Load Pauli string Hamiltonian from HamLib HDF5 file format.

    HamLib format stores operators as UTF-8 strings in the OpenFermion QubitOperator format:
    (coefficient+0j) [pauli_string] +

    Example: (1.5+0j) [X0 Z3] +\n(-0.5+0j) [Y1 Y2] +
    """
    filename = config_hamiltonian.filename
    logger.info(f"Loading HamLib HDF5 Hamiltonian from file \"{filename}\".")

    # Determine the HDF5 key/path - user can specify it or we'll try to find it
    hdf5_key = getattr(config_hamiltonian, 'hdf5_key', None)

    with h5py.File(filename, 'r') as f:
        # If no key specified, try to auto-detect
        if hdf5_key is None:
            # Get all dataset keys
            all_keys = []
            def collect_keys(name, obj):
                if isinstance(obj, h5py.Dataset):
                    all_keys.append(name)
            f.visititems(collect_keys)

            if len(all_keys) == 0:
                raise ValueError(f"No datasets found in HDF5 file \"{filename}\".")
            elif len(all_keys) == 1:
                hdf5_key = all_keys[0]
                logger.info(f"Auto-detected HDF5 key: \"{hdf5_key}\"")
            else:
                raise ValueError(
                    f"Multiple datasets found in HDF5 file. Please specify hdf5_key. "
                    f"Available keys: {all_keys}")

        # Load the dataset
        dataset = f[hdf5_key]

        # Read metadata if available (HamLib v1.1+)
        metadata = dict(dataset.attrs.items()) if hasattr(dataset, 'attrs') else {}
        if metadata:
            logger.info(f"HamLib metadata: {metadata}")

        # Read the Pauli string data as UTF-8
        hamlib_string = dataset[()].decode("utf-8")

    # Parse the HamLib format
    # Format: (coefficient+0j) [pauli_ops] +\n
    numq = 0
    pauli_dict = {}

    # Split into individual terms (separated by ' +\n')
    terms = hamlib_string.strip().split(' +\n')

    for term_str in terms:
        term_str = term_str.strip()
        if not term_str:
            continue

        # Find the coefficient part (between parentheses)
        paren_end = term_str.find(')')
        if paren_end == -1:
            raise ValueError(f"Invalid HamLib format: missing ')' in term: {term_str}")

        coef_str = term_str[1:paren_end]  # Extract content between ( )
        coefficient = complex(coef_str)

        # Find the Pauli string part (between brackets)
        bracket_start = term_str.find('[', paren_end)
        bracket_end = term_str.find(']', bracket_start)

        if bracket_start == -1 or bracket_end == -1:
            raise ValueError(f"Invalid HamLib format: missing brackets in term: {term_str}")

        pauli_str = term_str[bracket_start+1:bracket_end].strip()

        # Parse the sparse Pauli string (reusing existing logic)
        pauli_tokens = pauli_str.split() if pauli_str else []
        sparse_pauli = tuple()

        for token in pauli_tokens:
            op = token[0]  # Pauli operator: X, Y, or Z
            idx = int(token[1:])  # Qubit index
            numq = max(numq, idx + 1)
            sparse_pauli = (*sparse_pauli, (idx, op))

        pauli_dict[sparse_pauli] = coefficient

    # Validate that all coefficients are real (Hamiltonians must be Hermitian)
    # and convert to float
    for pauli in pauli_dict:
        coef = pauli_dict[pauli]
        # Use relative tolerance (scale-invariant)
        if abs(coef.imag) > abs(coef) * 1e-8:
            imag_ratio_percent = abs(coef.imag) / abs(coef) * 100
            raise ValueError(
                f"Hamiltonian must be Hermitian (real coefficients). "
                f"Found coefficient {coef} where imaginary part is "
                f"{imag_ratio_percent:.4g}% of magnitude (max allowed: 1e-6%).")
        pauli_dict[pauli] = coef.real

    logger.info(f"Loaded {len(pauli_dict)} Pauli terms on {numq} qubits.")

    return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, sparse=pauli_dict))

# -------------------------------------------------------------------------------------------------

def load_pauli(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading Pauli string Hamiltonian from file \"{filename}\".")
    extension = filename[filename.rfind('.')+1:]

    # Check for HamLib HDF5 format
    if extension in ["h5", "hdf5"]:
        return load_hamlib_hdf5(config_hamiltonian)
    elif extension in [ "txt", "dat" ]:
        fmt = None
        numq = 0
        pauli_dict = dict()
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                idx = line.find(' ')
                coef_str = line[:idx].strip()
                pauli = line[idx+1:].strip()
                if pauli[0] == '[':
                    if fmt is not None and fmt != "sparse":
                        raise ValueError("Inconsistent Pauli string file format.")
                    fmt = "sparse"
                    coefficient = complex(coef_str[1:-1])
                    if pauli[-1] == '+':
                        pauli = pauli[:pauli.rfind(']')+1]
                    pauli = pauli[1:-1]
                    pauli_tokens = pauli.split()
                    sparse_pauli = tuple()
                    for token in pauli_tokens:
                        op = token[0]
                        idx = int(token[1:])
                        numq = max(numq, idx+1)
                        sparse_pauli = (*sparse_pauli, (idx, op))
                    pauli_dict[sparse_pauli] = coefficient
                else:
                    if fmt is not None and fmt != "dense":
                        raise ValueError("Inconsistent Pauli string file format.")
                    fmt = "dense"
                    coefficient = complex(coef_str)
                    if numq != 0 and len(pauli) != numq:
                        raise ValueError("Inconsistent dense Pauli string length.")
                    numq = len(pauli)
                    pauli_dict[pauli] = coefficient
        # Validate that all coefficients are real (Hamiltonians must be Hermitian)
        # and convert to float
        for pauli in pauli_dict:
            coef = pauli_dict[pauli]
            # Use relative tolerance (scale-invariant)
            if abs(coef.imag) > abs(coef) * 1e-8:
                imag_ratio_percent = abs(coef.imag) / abs(coef) * 100
                raise ValueError(
                    f"Hamiltonian must be Hermitian (real coefficients). "
                    f"Found coefficient {coef} where imaginary part is "
                    f"{imag_ratio_percent:.4g}% of magnitude (max allowed: 1e-6%).")
            pauli_dict[pauli] = coef.real
        if fmt == "dense":
            return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, dense=pauli_dict))
        elif fmt == "sparse":
            return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, sparse=pauli_dict))
        else:
            raise ValueError(f"Invalid Pauli format: \"{fmt}\".")
    elif extension == "json":
        with open(filename, 'r') as file:
            data = json.load(file)

        # Validate structure
        if "n_qubits" not in data:
            raise ValueError("JSON Pauli file must contain 'n_qubits' field.")
        if "terms" not in data:
            raise ValueError("JSON Pauli file must contain 'terms' field.")

        numq = data["n_qubits"]
        pauli_dict = dict()

        for term in data["terms"]:
            if "ops" not in term or "coeff" not in term:
                raise ValueError("Each term must have 'ops' and 'coeff' fields.")

            # Convert ops list to sparse tuple format
            # ops is list of [index, operator] pairs
            sparse_pauli = tuple((idx, op) for idx, op in term["ops"])
            coefficient = term["coeff"]

            # Convert to complex for validation
            if isinstance(coefficient, (int, float)):
                coefficient = complex(coefficient)
            elif isinstance(coefficient, complex):
                pass
            else:
                raise ValueError(f"Coefficient must be numeric, got {type(coefficient)}.")

            pauli_dict[sparse_pauli] = coefficient

        # Validate Hermitian and convert to float
        for pauli in pauli_dict:
            coef = pauli_dict[pauli]
            if abs(coef.imag) > abs(coef) * 1e-8:
                imag_ratio_percent = abs(coef.imag) / abs(coef) * 100
                raise ValueError(
                    f"Hamiltonian must be Hermitian (real coefficients). "
                    f"Found coefficient {coef} where imaginary part is "
                    f"{imag_ratio_percent:.4g}% of magnitude (max allowed: 1e-6%).")
            pauli_dict[pauli] = coef.real

        return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, sparse=pauli_dict))
    else:
        raise ValueError(
            f"Invalid file extension for loading a Pauli string file: \"{extension}\".")

# -------------------------------------------------------------------------------------------------

# TODO: Scott and I have both spent time chasing down the types of different things for a variety
#       of reasons.  If we can get this code to the point where it always returns the same type
#       regardless of the options passed in, then we should annotate the return type.  If it turns
#       out multiple different Python types can be returned from this function because of reliance
#       on duck typing, then any annotation would at best be a comment listing the different types,
#       and that may or may not be as useful.
#           `-> tuple[???,???]`
def get_physical_hamiltonian(config_hamiltonian: HamiltonianConfiguration):

    logger.info("Beginning `get_physical_hamiltonian()` function.")

    if config_hamiltonian.source == "numpy":
        return load_numpy(config_hamiltonian)
    elif config_hamiltonian.source == "LCPS":
        return load_LCPS(config_hamiltonian)
    elif config_hamiltonian.source == "hdf5":
        return load_hdf5(config_hamiltonian)
    elif config_hamiltonian.source == "pauli":
        return load_pauli(config_hamiltonian)
    else:
        raise ValueError(f"Invalid Hamiltonian source \"{config_hamiltonian.source}\".")
