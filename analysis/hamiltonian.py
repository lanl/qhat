import h5py
import json
import logging
import numpy as np
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from openfermion import (
    InteractionOperator,
    QubitOperator,
    bravyi_kitaev,
    jordan_wigner,
)
from openfermion.chem.molecular_data import spinorb_from_spatial

from pyscf import ao2mo, tools

from qhat.analysis.config_types import HamiltonianConfiguration
from qhat.analysis.utils import value
from qhat.common.bosons_binary import BosonicBinaryEncoding
from qhat.common.hermiticity import (
    DEFAULT_RTOL,
    check_hermitian_one_body,
    check_hermitian_two_body,
    real_pauli_coefficients,
    real_pauli_file_coefficients,
    require_real_scalar,
)
from qhat.common.MixedFermionBosonOperator import MixedFermionBosonOperator
from qhat.common.terms.base_term import BaseTerm
from qhat.common.terms.single_pauli_term import SinglePauliTerm

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

def _validate_fermionic_tensors(constant, one_body, two_body):
    """Check second-quantization tensors; return (real constant, number of spin orbitals)."""
    if len(one_body.shape) != 2 or one_body.shape[0] != one_body.shape[1]:
        raise ValueError(f"One-body tensor must be square, got shape {one_body.shape}.")
    num_orbitals = one_body.shape[0]
    if tuple(two_body.shape) != (num_orbitals,) * 4:
        raise ValueError(
            f"Two-body tensor must have shape {(num_orbitals,) * 4}, got {two_body.shape}.")
    constant = require_real_scalar("constant energy term", constant)
    check_hermitian_one_body("one-body tensor", one_body)
    check_hermitian_two_body("two-body tensor", two_body)
    return constant, num_orbitals

def _pauli_terms_from_qubit_operator(qubit_operator, num_qubits):
    # The qubit mappings return complex coefficients (with round-off imaginary parts) when the
    # integrals are complex, even for a Hermitian Hamiltonian; make them real here.
    pauli_strings = real_pauli_coefficients(qubit_operator.terms)
    return [SinglePauliTerm(ps, c, num_qubits) for ps, c in pauli_strings.items()]

# -------------------------------------------------------------------------------------------------

class Hamiltonian:
    def __init__(self, terms, num_qubits, core_operator=None):
        """
        Construct a Hamiltonian from a collection of terms and a number of qubits (optionally also
        a core operator).

        The Hamiltonian is assumed to be the sum of the terms.  While addition of operators is
        commutative, the order does matter because some methods will transform the terms into
        time-evolution operators that multiply together, and time-evolution operators are not
        commutative under multiplication.
        """
        self._terms = list(terms)
        self._num_qubits = num_qubits
        self._core_operator = core_operator
        self._energy_shift_total = 0.0
        if not self._terms:
            raise ValueError("A Hamiltonian must have at least one term.")
        for term in self._terms:
            # All terms must actually be terms
            if not isinstance(term, BaseTerm):
                raise TypeError("All terms in the Hamiltonian must derive from BaseTerm.")
            # Ensure the number of qubits matches
            # -- They don't _actually_ need to match, because each term could use a subset of the
            #    qubits.  But to simplify the logic, we're just going to require that the number of
            #    qubits match and the consistency check and stitching together of terms is trivial.
            if term.num_qubits() != self._num_qubits:
                raise ValueError(f"Mismatch in number of qubits.  Term {term} has "
                                 f"{term.num_qubits()} qubits but should have {self._num_qubits}.")
    def num_qubits(self):
        """
        The number of qubits that this Hamiltonian operates on.
        """
        return self._num_qubits
    def terms(self):
        """
        Returns an ordered container of the terms in the Hamiltonian.
        """
        return self._terms
    # TODO: Rename to just `core_operator` for consistency with num_qubits() and terms()?
    def get_core_operator(self):
        """
        Returns the core operator, which may be needed by downstream algorithms (such as block
        encodings from pyLIQTR).
        """
        # TODO: Do we really need to store the core_operator, or can it be constructed on demand?
        if self._core_operator is None:
            raise ValueError(
                "This Hamiltonian was not generated with a core operator, so requesting the core "
                "operator is an invalid operation.")
        return self._core_operator
    def to_pauli_strings(self):
        # TODO: This is a temporary bridge method to be removed once Pauli-based analyses
        #       (eigenvalue bounds, nested-commutator bounds, etc) are updated to use the term
        #       interface.
        pauli_strings = {}
        for term in self._terms:
            term_dict = term.to_pauli_strings()
            for pstr, coeff in term_dict.items():
                pauli_strings[pstr] = pauli_strings.get(pstr, 0) + coeff
        return pauli_strings
    def energy_shift(self, dE):
        """
        Shifts the energy spectrum of the Hamiltonian by shifting the identity term.
        """
        dE = require_real_scalar("energy shift", dE)
        # Track the cumulative shift for analysis corrections
        self._energy_shift_total += dE
        # Update the identity term
        # TODO: This suggests that an identity term type would be useful, especially if we restrict
        #       all other term types to prevent them from being identity.
        for i, term in enumerate(self._terms):
            if isinstance(term, SinglePauliTerm) and term.pauli_string.to_sparse() == ():
                new_coefficient = term.coefficient + dE
                # Shifting to centre the spectrum cancels the identity coefficient; drop the term
                # when only round-off remains (but never leave the Hamiltonian with no terms).
                cancelled = (abs(new_coefficient)
                             <= DEFAULT_RTOL * max(abs(term.coefficient), abs(dE)))
                if cancelled and len(self._terms) > 1:
                    del self._terms[i]
                else:
                    self._terms[i] = SinglePauliTerm(term.pauli_string, new_coefficient)
                break
        else:
            if dE != 0.0:
                self._terms.insert(0, SinglePauliTerm((), dE, self._num_qubits))
        # TODO: If we can reconstruct the core operator on demand, then we don't need to update it
        #       here, preventing inconsistencies and eliminating the need to worry if it somehow
        #       was set to an invalid value.
        if isinstance(self._core_operator, InteractionOperator):
            t0 = self._core_operator.constant + dE
            t1 = self._core_operator.one_body_tensor
            t2 = self._core_operator.two_body_tensor
            self._core_operator = InteractionOperator(t0, t1, t2)
        elif isinstance(self._core_operator, MixedFermionBosonOperator):
            self._core_operator.energy_shift(dE)
        elif self._core_operator is not None:
            raise TypeError(
                "Unable to apply an energy shift a Hamiltonian with core_operator of type "
                f"\"{type(self._core_operator)}\".")
    # TODO: Rename to just `energy_shift` for consistency with `terms()` and `num_qubits()`?  But
    #       `energy_shift()` already exists as the function that updates the energy shift.
    def get_energy_shift(self) -> float:
        """
        Get the cumulative energy shift applied to this Hamiltonian.

        This is used by analysis functions to correct eigenvalues and errors
        for meaningful comparisons between exact and approximate representations.

        Returns:
            The total energy shift applied via energy_shift() calls
        """
        return self._energy_shift_total
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

        # TODO: We'll need to add a to_matrix method to each term to support this feature
        # Get Pauli strings as dense format (full string per term)
        pauli_dict = {pstr.to_dense(): coeff for pstr, coeff in self.to_pauli_strings().items()}
        num_qubits = self.num_qubits()

        logger.verbose(f"Converting Hamiltonian to matrix for {num_qubits} qubits")
        logger.verbose(f"Hamiltonian has {len(pauli_dict)} Pauli terms")

        # Create operator (dense or sparse based on memory threshold)
        return create_hamiltonian_operator(
            pauli_dict, num_qubits, memory_threshold_gb,
            force_dense=force_dense,
            force_sparse=force_sparse
        )
    def compute_initial_energy_bounds(
            self,
            config_hamiltonian: HamiltonianConfiguration):
        # TODO: This is one of the uses of to_pauli_strings that we need to address in the future
        logger.info("Computing initial energy bounds.")
        pauli_sum = {pstr.to_sparse() : coeff for pstr, coeff in self.to_pauli_strings().items()}
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
        Elo1 = require_real_scalar("lower energy bound", Elo1)
        Ehi1 = require_real_scalar("upper energy bound", Ehi1)
        logger.verbose(f"-- initial bounds = [{Elo1}, {Ehi1})")
        return (Elo1, Ehi1)
    @classmethod
    def from_pauli_strings(cls, pauli_strings, num_qubits):
        # TODO: This assumes that pauli_strings is a dict of {pauli_string : coefficient, ...}, but
        #       different types for the Pauli strings are supported because the SinglePauliString
        #       constructor handles that for us.  Do we want to support the tuple((pauli_string,
        #       coefficient), ...) format as well?
        #    -- The main complication there is deduplication without overwriting a previous value.
        #    -- Currently the way to get a duplicate value is if you use the sparse format and
        #       switch the order (e.g.: ((3,'Z'),(0,'X')) vs ((0,'X'),(3,'Z')) would be the same
        #       Pauli string but hash to different dict keys) or to mix different formats (e.g.,
        #       tuple vs str vs PauliString could represent the same Pauli string but would hash to
        #       different dict keys).
        pauli_terms = [SinglePauliTerm(ps, c, num_qubits) for ps, c in pauli_strings.items()]
        # Different keys can represent the same Pauli string (see above).  This washes out once
        # SinglePauliTerm canonicalizes to PauliString, so combine coefficients on that key.
        combined = {}
        for term in pauli_terms:
            combined[term.pauli_string] = combined.get(term.pauli_string, 0.0) + term.coefficient
        return cls((SinglePauliTerm(ps, c) for ps, c in combined.items()), num_qubits)
    # TODO: Revisit whether from_second_quantization and from_mixed_fermion_boson should merge.
    #       Mixed fermion-boson input is still second quantization, just with additional tensors.
    #       One option is a single from_second_quantization with optional bosonic kwargs (or a
    #       grouped bosonic-data argument).  Revisit once the Hamiltonians we need are better
    #       understood.
    @classmethod
    def from_second_quantization(cls, constant, one_body, two_body, fermion_mapping):
        """
        Construct a Hamiltonian from fermionic second-quantization tensors in the OpenFermion
        convention.  `fermion_mapping` is a name in `fermionic_mapping` or a callable that maps an
        InteractionOperator to a QubitOperator.
        """
        constant, _ = _validate_fermionic_tensors(constant, one_body, two_body)
        if isinstance(fermion_mapping, str):
            fermion_mapping = fermionic_mapping[fermion_mapping]
        operator = InteractionOperator(constant, one_body, two_body)
        terms = _pauli_terms_from_qubit_operator(fermion_mapping(operator), operator.n_qubits)
        return cls(terms, operator.n_qubits, core_operator=operator)
    @classmethod
    def from_mixed_fermion_boson(
            cls, constant, one_body, two_body, boson_number_coefficient, interaction,
            fermion_mapping, boson_mapping, max_bosons_per_state):
        """
        Construct a Hamiltonian from mixed fermion-boson second-quantization tensors.  The
        fermionic tensors are in the OpenFermion convention, and `interaction` has shape
        (Nf, Nf, Nb).  `fermion_mapping` is a name in `fermionic_mapping` or a callable;
        `boson_mapping` is a name in `bosonic_mapping` (which uses `max_bosons_per_state`) or an
        encoding object.
        """
        constant, num_orbitals = _validate_fermionic_tensors(constant, one_body, two_body)
        if (len(interaction.shape) != 3
                or tuple(interaction.shape[:2]) != (num_orbitals, num_orbitals)
                or interaction.shape[2] == 0):
            raise ValueError(
                f"Fermion-boson interaction tensor must have shape ({num_orbitals}, "
                f"{num_orbitals}, Nb) with Nb > 0, got {interaction.shape}.")
        boson_number_coefficient = require_real_scalar(
            "bosonic number coefficient", boson_number_coefficient)
        for n in range(interaction.shape[2]):
            check_hermitian_one_body(
                f"fermion-boson interaction tensor [:,:,{n}]", interaction[:, :, n])
        if isinstance(fermion_mapping, str):
            fermion_mapping = fermionic_mapping[fermion_mapping]
        if isinstance(boson_mapping, str):
            if max_bosons_per_state is None:
                raise ValueError("Maximum bosons per state is required for a bosonic mapping.")
            boson_mapping = bosonic_mapping[boson_mapping](max_bosons_per_state)
        operator = MixedFermionBosonOperator(
            constant, one_body, two_body, boson_number_coefficient, interaction)
        operator.set_fermionic_encoding(fermion_mapping)
        operator.set_bosonic_encoding(boson_mapping)
        terms = _pauli_terms_from_qubit_operator(
            operator.generate_qubit_operator(), operator.num_qubits())
        return cls(terms, operator.num_qubits(), core_operator=operator)

# -------------------------------------------------------------------------------------------------

def _convert_to_openfermion(f1in, f2in, integral_convention):
    if integral_convention == "openfermion":
        # no conversion needed
        return (f1in, f2in)
    elif integral_convention == "chemist_spatial" or integral_convention == "chemist_spin":
        # absorb factor of 0.5 and reorder indices chemist -> physicist
        f2trans = 0.5 * f2in.transpose(0, 2, 3, 1)
        # transform from spatial to spin orbitals, if needed
        if integral_convention == "chemist_spatial":
            f1, f2 = spinorb_from_spatial(f1in, f2trans)
        else:
            f1, f2 = (f1in, f2trans)
        return (f1, f2)
    else:
        raise ValueError(f"Invalid integral_convention '{integral_convention}'")

# -------------------------------------------------------------------------------------------------

def _construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb, integral_convention="openfermion"):
    f1, f2 = _convert_to_openfermion(f1, f2, integral_convention)
    fermion_mapping = config_hamiltonian.fermion_to_qubit_transform
    if bs is None and fb is None:
        # pure-fermion case
        H = Hamiltonian.from_second_quantization(f0, f1, f2, fermion_mapping)
        logger.info(f"Fermionic second-quantization Hamiltonian uses {H.num_qubits()} qubits.")
        logger.info(f"Mapped fermionic operators to qubit operators using {fermion_mapping} method.")
        return H
    else:
        # mixed fermion-boson case (currently no pure-boson case)
        if bs is None or fb is None:
            raise ValueError(
                    "Second-quantization Hamiltonians with bosons must provide both the "
                    "bosonic scalar term and the fermion-boson interaction term.")
        if config_hamiltonian.max_bosons_per_state is None:
            raise ValueError(
                "Hamiltonian includes bosons but user did not specify maximum bosons per state.")
        boson_mapping = config_hamiltonian.boson_to_qubit_transform
        max_bosons = config_hamiltonian.max_bosons_per_state
        H = Hamiltonian.from_mixed_fermion_boson(
                f0, f1, f2, bs, fb, fermion_mapping, boson_mapping, max_bosons)
        logger.info(
            f"Mixed fermionic-bosonic second-quantization Hamiltonian uses {H.num_qubits()} qubits.")
        logger.info(f"Mapped fermionic operators to qubit operators using {fermion_mapping} method.")
        logger.info(f"Mapped bosonic operators to qubit operators using {boson_mapping} method "
                    f"with {max_bosons} maximum bosons per state.")
        return H

# -------------------------------------------------------------------------------------------------

def load_hdf5(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading second-quantization Hamiltonian from HDF5 file \"{filename}\".")
    logger.info(f"Ignoring integral_convention field and using default (openfermion).")
    data = h5py.File(filename)
    f0 = 0      # Currently don't support constant terms in HDF5
    f1 = data["1e"]
    f2 = data["2e"]
    bs = None   # Currently don't support bosons in HDF5
    fb = None   # Currently don't support bosons in HDF5
    return _construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb)

# -------------------------------------------------------------------------------------------------

_no_argument = object()

def _get_any_of(data: Mapping[str, np.ndarray], keys: list[str], as_scalar: bool = False, default: Any = _no_argument) -> Any:
    """Return data[k] for any k in keys for which that value is defined.
    Raise KeyError on failure unless a default value is provided.  If
    as_scalar is True, extract a scalar from a 0-D array."""

    # Try each key in turn.
    for k in keys:
        try:
            v = data[k]
            return v[()] if as_scalar else v
        except KeyError:
            pass

    # No keys were found.  Return None if the value is optional, otherwise
    # a KeyError.
    if default is not _no_argument:
        return default
    raise KeyError(f"failed to find any of the keys {keys}")

# -------------------------------------------------------------------------------------------------

def load_numpy(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    integral = config_hamiltonian.integral_convention
    logger.info(f"Loading second-quantization Hamiltonian from file \"{filename}\".")
    logger.info(f"Using integral_convention field '{integral}'")
    data = np.load(filename)
    f0 = _get_any_of(data, ["constant", "ecore"], as_scalar=True, default=0.0)
    f1 = _get_any_of(data, ["one_body", "h1e"])
    f2 = _get_any_of(data, ["two_body", "eri"])
    bs = _get_any_of(data, ["bosonic_scalar"], default=None)
    fb = _get_any_of(data, ["fb_interaction"], default=None)
    return _construct_second_quantization(config_hamiltonian, f0, f1, f2, bs, fb, integral)

# -------------------------------------------------------------------------------------------------

def load_fci(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading second-quantization Hamiltonian from file \"{filename}\".")
    logger.info(f"Ignoring integral_convention field and using chemist_spatial.")
    data = tools.fcidump.read(filename)

    n = data['NORB']
    f0 = data['ECORE']

    # check that f1 has correct shape
    f1 = data['H1']
    if f1.shape != (n, n):
        raise ValueError(f"Incorrect f1 shape {f1.shape}. "
                         f"f1 shape should be ({n}, {n}).")

    # reader returns a packed f2 - unpack it here
    f2packed = data['H2']
    f2 = ao2mo.restore('s1', f2packed, n)

    # check that f2 has correct shape
    if f2.shape != (n, n, n, n):
        raise ValueError(f"Incorrect f2 shape {f2.shape}. "
                         f"f2 shape should be ({n}, {n}, {n}, {n}).")

    return _construct_second_quantization(config_hamiltonian, f0, f1, f2, None, None, "chemist_spatial")

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

    pauli_dict = real_pauli_file_coefficients(pauli_dict)

    logger.info(f"Loaded {len(pauli_dict)} Pauli terms on {numq} qubits.")

    return Hamiltonian.from_pauli_strings(pauli_dict, numq)

# -------------------------------------------------------------------------------------------------

def load_pauli(config_hamiltonian: HamiltonianConfiguration):
    filename = config_hamiltonian.filename
    logger.info(f"Loading Pauli string Hamiltonian from file \"{filename}\".")
    extension = Path(filename).suffix[1:]  # Remove leading dot.

    # Check for HamLib HDF5 format
    if extension in ["h5", "hdf5"]:
        return load_hamlib_hdf5(config_hamiltonian)
    elif extension in ["txt", "dat"]:
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
        pauli_dict = real_pauli_file_coefficients(pauli_dict)
        return Hamiltonian.from_pauli_strings(pauli_dict, numq)
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

        pauli_dict = real_pauli_file_coefficients(pauli_dict)

        return Hamiltonian.from_pauli_strings(pauli_dict, numq)
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
    elif config_hamiltonian.source == "hdf5":
        return load_hdf5(config_hamiltonian)
    elif config_hamiltonian.source == "pauli":
        return load_pauli(config_hamiltonian)
    elif config_hamiltonian.source == "fci":
        return load_fci(config_hamiltonian)
    else:
        raise ValueError(f"Invalid Hamiltonian source \"{config_hamiltonian.source}\".")
