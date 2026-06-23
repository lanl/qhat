"""JW/BK Trotter-isomorphism comparison tools for QHAT.

This package adds reusable comparison utilities around QHAT's existing Pauli
matrix and Trotterization helpers.  OpenFermion-dependent functions are loaded
only when requested so importing the package remains lightweight.
"""

from importlib import import_module

from qhat.trotter_isomorphism.comparison import (
    bitstring_overlap_comparison,
    compare_jw_bk_operator_norm_errors,
    compare_state_dependent_trotter_error,
    fidelity_from_overlap_amplitude,
    infidelity_from_overlap_amplitude,
    normalize_state,
    parse_bitstring,
    trotter_overlap_amplitude,
)
from qhat.trotter_isomorphism.pauli_keys import (
    anticommutes,
    clean_qubit_operator,
    dense_string_to_pauli_key,
    format_pauli_key,
    pauli_key,
    pauli_key_to_dense_string,
    pauli_key_to_matrix,
    pauli_weight,
    qubit_operator_to_matrix,
)
from qhat.trotter_isomorphism.product_formula import (
    dense_trotter_unitary,
    exact_unitary,
    first_order_trotter_error_matrix,
    first_order_trotter_error_norm,
    first_order_trotter_unitary,
    second_order_trotter_error_matrix,
    second_order_trotter_error_norm,
    second_order_trotter_unitary,
    trotter_error_matrix,
    trotter_error_norm,
)
from qhat.trotter_isomorphism.states import bitstring_state, bits_from_index, index_from_bits

__version__ = "0.1.0"

_LAZY_EXPORTS = {
    "JWBKAnalysis": ("qhat.trotter_isomorphism.jw_bk", "JWBKAnalysis"),
    "build_clifford_from_z_map": (
        "qhat.trotter_isomorphism.jw_bk",
        "build_clifford_from_z_map",
    ),
    "build_jw_bk_analysis": ("qhat.trotter_isomorphism.jw_bk", "build_jw_bk_analysis"),
    "build_matched_term_lists": (
        "qhat.trotter_isomorphism.jw_bk",
        "build_matched_term_lists",
    ),
    "count_qubits_in_fermion_operator": (
        "qhat.trotter_isomorphism.jw_bk",
        "count_qubits_in_fermion_operator",
    ),
    "find_jw_bk_isomorphism": (
        "qhat.trotter_isomorphism.jw_bk",
        "find_jw_bk_isomorphism",
    ),
    "verify_graph_isomorphism": (
        "qhat.trotter_isomorphism.jw_bk",
        "verify_graph_isomorphism",
    ),
    "fermion_operator_to_majorana_operator_dict": (
        "qhat.trotter_isomorphism.majorana",
        "fermion_operator_to_majorana_operator_dict",
    ),
    "fermion_term_to_majorana_terms": (
        "qhat.trotter_isomorphism.majorana",
        "fermion_term_to_majorana_terms",
    ),
    "majorana_monomial_as_fermion_operator": (
        "qhat.trotter_isomorphism.majorana",
        "majorana_monomial_as_fermion_operator",
    ),
    "single_majorana_as_fermion_operator": (
        "qhat.trotter_isomorphism.majorana",
        "single_majorana_as_fermion_operator",
    ),
}

__all__ = [
    "__version__",
    "JWBKAnalysis",
    "anticommutes",
    "bitstring_state",
    "bitstring_overlap_comparison",
    "bits_from_index",
    "build_clifford_from_z_map",
    "build_jw_bk_analysis",
    "build_matched_term_lists",
    "compare_jw_bk_operator_norm_errors",
    "compare_state_dependent_trotter_error",
    "clean_qubit_operator",
    "count_qubits_in_fermion_operator",
    "dense_string_to_pauli_key",
    "dense_trotter_unitary",
    "exact_unitary",
    "fermion_operator_to_majorana_operator_dict",
    "fermion_term_to_majorana_terms",
    "fidelity_from_overlap_amplitude",
    "find_jw_bk_isomorphism",
    "first_order_trotter_error_matrix",
    "first_order_trotter_error_norm",
    "first_order_trotter_unitary",
    "format_pauli_key",
    "index_from_bits",
    "infidelity_from_overlap_amplitude",
    "majorana_monomial_as_fermion_operator",
    "normalize_state",
    "parse_bitstring",
    "pauli_key",
    "pauli_key_to_dense_string",
    "pauli_key_to_matrix",
    "pauli_weight",
    "qubit_operator_to_matrix",
    "second_order_trotter_error_matrix",
    "second_order_trotter_error_norm",
    "second_order_trotter_unitary",
    "single_majorana_as_fermion_operator",
    "trotter_error_matrix",
    "trotter_error_norm",
    "trotter_overlap_amplitude",
    "verify_graph_isomorphism",
]


def __getattr__(name):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _LAZY_EXPORTS[name]
    module = import_module(module_name)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value
