from common.bosons_binary import BosonicBinaryEncoding
from common.MixedFermionBosonOperator import MixedFermionBosonOperator
from qre_types import GeneralConfiguration, HamiltonianConfiguration, value

from functools import cache, reduce
import h5py
import numpy as np
from openfermion import InteractionOperator, QubitOperator, count_qubits, bravyi_kitaev, \
                        jordan_wigner, binary_code_transform
from pyscf import scf, gto, lib, ao2mo
import scipy.constants as sc

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

# TODO: These are generally useful utilities.  Where should they live?
def sparse_to_dense_pauli(sparse_pauli, num_qubits):
    dense_pauli = ["I",] * num_qubits
    for idx, op in sparse_pauli:
        dense_pauli[idx] = op
    return "".join(dense_pauli)
def dense_to_sparse_pauli(dense_pauli):
    sparse_pauli = tuple()
    for idx, op in enumerate(dense_pauli):
        if op in ["X", "Y", "Z"]:
            sparse_pauli = (*sparse_pauli, (idx,op))
        elif op != "I":
            raise ValueError(f"Invalid character in dense pauli string: \"{op}\".")
    return sparse_pauli

# -------------------------------------------------------------------------------------------------

class LinearCombinationOfPauliStrings:
    def __init__(self, **kwargs):
        self._nq = None
        self._format = None
        self._data = None
        self._nq = kwargs["num_qubits"]
        for f in [ "dense", "sparse" ]:
            if f in kwargs:
                if self._format is not None:
                    raise ValueError(
                        "Too many formats provided to LinearCombinationOfPauliStrings.")
                self._format = f
                self._data = kwargs[f]
                assert isinstance(self._data, dict)
        if self._format is None:
            raise ValueError("No data provided to LinearCombinationOfPauliStrings.")
    def num_qubits(self):
        raise NotImplementedError()
    def get_dense_pauli_strings(self):
        if self._format == "dense":
            return self._data
        elif self._format == "sparse":
            return {sparse_to_dense_pauli(pauli, self._nq) : coef
                    for pauli, coef in self._data.items()}
        else:
            raise ValueError("Invalid data format \"{self._format}\".")
    def get_sparse_pauli_strings(self):
        if self._format == "dense":
            return {dense_to_sparse_pauli(pauli) : coef for pauli, coef in self._data.items()}
        elif self._format == "sparse":
            return self._data
        else:
            raise ValueError("Invalid data format \"{self._format}\".")
    def energy_shift(self, shift):
        all_identity = tuple()
        if self._format == "dense":
            all_identity = sparse_to_dense_pauli(all_identity, self._nq)
        identity_coefficient = self._data.get(all_identity, 0.0) + shift
        self._data[all_identity] = identity_coefficient

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
            if isinstance(self._H, LinearCombinationOfPauliStrings):
                return self._H.get_dense_pauli_strings()
            else:
                as_tuples = self.get_all_pauli_strings(return_as="tuples")
                Nq = self.num_qubits()
                # TODO: Why is the all-identity string excluded?  Isn't that a bug?
                return {sparse_to_dense_pauli(pauli, Nq) : coef
                    for pauli, coef in as_tuples.items() if pauli != ()}
        else:
            raise ValueError("  ".join([
                "The value of return_as must be \"tuples\" or \"strings\".",
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
            config_general: GeneralConfiguration,
            config_hamiltonian: HamiltonianConfiguration):
        config_general.log("Computing initial energy bounds.")
        pauli_sum = self.get_all_pauli_strings()
        config_general.log_verbose(f"-- number of Pauli strings = {len(pauli_sum)}")
        energy_shift = pauli_sum[()] # the encoding only lists non-identity matrices, so () = I
        dE = sum(abs(coefficient) for coefficient in pauli_sum.values()) - abs(energy_shift)
        Elo0 = energy_shift - dE
        Ehi0 = energy_shift + dE
        config_general.log_verbose(f"-- energy shift = {energy_shift}")
        config_general.log_verbose(f"-- computed bounds = [{Elo0}, {Ehi0})")
        EloU = value(config_hamiltonian.lower_bound, float('-inf'))
        EhiU = value(config_hamiltonian.upper_bound, float('inf'))
        Elo1 = max(Elo0, EloU)
        Ehi1 = min(Ehi0, EhiU)
        config_general.log_verbose(f"-- limited bounds = [{Elo1}, {Ehi1})")
        if config_hamiltonian.exact_energy_lower_bound:
            assert config_hamiltonian.lower_bound is not None
            Elo1 = config_hamiltonian.lower_bound
        if config_hamiltonian.exact_energy_upper_bound:
            assert config_hamiltonian.upper_bound is not None
            Ehi1 = config_hamiltonian.upper_bound
        config_general.log_verbose(f"-- initial bounds = [{Elo1}, {Ehi1})")
        return (Elo1, Ehi1)

# -------------------------------------------------------------------------------------------------

def _verify_and_construct_second_quantization(
        config_general, config_hamiltonian, f0, f1, f2, bs, fb):
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
        H.set_fermionic_mapping(config_hamiltonian.fermion_to_qubit_transform)
        config_general.log(
                f"fermionic second-quantization Hamiltonian uses {H.num_qubits()} qubits.")
        return H
    else:
        assert bs is not None
        assert fb is not None
        assert Nb > 0
        H = Hamiltonian(MixedFermionBosonOperator(f0, f1, f2, bs, fb))
        H.set_fermionic_mapping(config_hamiltonian.fermion_to_qubit_transform)
        H.set_bosonic_mapping(
                config_hamiltonian.boson_to_qubit_transform,
                config_hamiltonian.max_bosons_per_state)
        config_general.log(" ".join([
            "mixed fermionic-bosonic second-quantization Hamiltonian",
            f"uses {H.num_qubits()} qubits."]))
        return H

# -------------------------------------------------------------------------------------------------

def load_hdf5(
        config_general: GeneralConfiguration,
        config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    config_general.log(
            f"Loading second-quantization Hamiltonian from HDF5 file \"{filename}\".")
    data = h5py.File(filename)
    f0 = 0      # Currently don't support constant terms in HDF5
    f1 = data["1e"]
    f2 = data["2e"]
    bs = None   # Currently don't support bosons in HDF5
    fb = None   # Currently don't support bosons in HDF5
    return _verify_and_construct_second_quantization(
            config_general, config_hamiltonian, f0, f1, f2, bs, fb)

# -------------------------------------------------------------------------------------------------

def load_numpy(
        config_general: GeneralConfiguration,
        config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    config_general.log(
            f"Loading second-quantization Hamiltonian from file \"{filename}\".")
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
    return _verify_and_construct_second_quantization(
            config_general, config_hamiltonian, f0, f1, f2, bs, fb)

# -------------------------------------------------------------------------------------------------

def load_pauli(
        config_general: GeneralConfiguration,
        config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    config_general.log(
            f"Loading Pauli string Hamiltonian from file \"{filename}\".")
    extension = filename[filename.rfind('.')+1:]
    if extension in [ "txt", "dat" ]:
        fmt = None
        numq = 0
        pauli_dict = dict()
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line[0] == "#":
                    continue
                idx = line.find(' ')
                coef_str = line[:idx].strip()
                pauli = line[idx+1:].strip()
                if pauli[0] == '[':
                    if fmt is not None and fmt != "sparse":
                        raise ValueError("Inconsistent Pauli string file format.")
                    fmt = "sparse"
                    coefficient = np.complex128(coef_str[1:-1])
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
                    coefficient = np.complex128(coef_str)
                    if numq != 0 and len(pauli) != numq:
                        raise ValueError("Inconsistent dense Pauli string length.")
                    numq = len(pauli)
                    pauli_dict[pauli] = coefficient
        if fmt == "dense":
            return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, dense=pauli_dict))
        elif fmt == "sparse":
            return Hamiltonian(LinearCombinationOfPauliStrings(num_qubits=numq, sparse=pauli_dict))
        else:
            raise ValueError(f"Invalid Pauli format: \"{fmt}\".")
    elif extension == "json":
        raise NotImplementedError("JSON Pauli string file not yet implemented.")
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
def get_physical_hamiltonian(
        config_general: GeneralConfiguration,
        config_hamiltonian: HamiltonianConfiguration):

    config_general.log("Beginning `get_physical_hamiltonian()` function.")

    if config_hamiltonian.source == "numpy":
        return load_numpy(config_general, config_hamiltonian)
    elif config_hamiltonian.source == "LCPS":
        return load_LCPS(config_general, config_hamiltonian)
    elif config_hamiltonian.source == "hdf5":
        return load_hdf5(config_general, config_hamiltonian)
    elif config_hamiltonian.source == "pauli":
        return load_pauli(config_general, config_hamiltonian)
    else:
        raise ValueError(f"Invalid Hamiltonian source \"{config_hamiltonian.source}\".")
