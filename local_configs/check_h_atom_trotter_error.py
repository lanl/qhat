import pickle
import numpy as np
from scipy.linalg import norm

from openfermion.transforms import get_fermion_operator

from qhat.trotter_isomorphism.jw_bk import (
    build_clifford_from_z_map,
    build_jw_bk_analysis,
    verify_graph_isomorphism,
)
from qhat.trotter_isomorphism.comparison import (
    compare_jw_bk_operator_norm_errors,
    bitstring_overlap_comparison,
)
from qhat.trotter_isomorphism.pauli_keys import format_pauli_key


active_pickle = "/tmp/qhat_h_atom_smoke/h_atom_as-001-001.pickle"

time = 1.0
num_steps = 1
method = 1

# Hydrogen active space has one occupied spin orbital and one vacant spin orbital.
# In JW occupation-bit notation, use |10>.
psi0_bitstring = "10"

with open(active_pickle, "rb") as f:
    active_hamiltonian = pickle.load(f)

fermion_op = get_fermion_operator(active_hamiltonian)

analysis = build_jw_bk_analysis(
    fermion_op,
    n_qubits=active_hamiltonian.n_qubits,
    atol=1e-10,
    remove_identity=False,
    build_matrices=True,
)

C_jw_to_bk, B, A = build_clifford_from_z_map(
    analysis.jw_to_bk,
    analysis.n_qubits,
)

print("Hydrogen atom JW/BK Trotter-error smoke test")
print("============================================")
print(f"active pickle = {active_pickle}")
print(f"n_qubits      = {analysis.n_qubits}")
print(f"time          = {time}")
print(f"num_steps     = {num_steps}")
print(f"method        = {method}")
print(f"psi0 JW       = |{psi0_bitstring}>")
print()

print("Graph-isomorphism check")
print("----------------------")
print(f"JW vertices = {len(analysis.jw_to_bk)}")
print(f"BK vertices = {len(analysis.bk_to_jw)}")
print(f"verified    = {verify_graph_isomorphism(analysis.jw_to_bk)}")
print()

print("Matched JW/BK term order")
print("------------------------")
for (jw_key, jw_coeff), (bk_key, bk_coeff) in zip(
    analysis.jw_terms,
    analysis.bk_terms,
):
    print(
        f"JW {format_pauli_key(jw_key):8s} {jw_coeff.real:+.16e}"
        f"   -->   BK {format_pauli_key(bk_key):8s} {bk_coeff.real:+.16e}"
    )
print()

print("Clifford sanity check")
print("---------------------")
unitarity_error = norm(
    C_jw_to_bk.conj().T @ C_jw_to_bk - np.eye(2**analysis.n_qubits),
    ord=2,
)
hamiltonian_conjugation_error = norm(
    C_jw_to_bk @ analysis.H_jw @ C_jw_to_bk.conj().T - analysis.H_bk,
    ord=2,
)
print(f"||C^dag C - I||_2         = {unitarity_error:.16e}")
print(f"||C H_JW C^dag - H_BK||_2 = {hamiltonian_conjugation_error:.16e}")
print()

operator_result = compare_jw_bk_operator_norm_errors(
    H_jw=analysis.H_jw,
    H_bk=analysis.H_bk,
    jw_terms=analysis.jw_terms,
    bk_terms=analysis.bk_terms,
    time=time,
    num_steps=num_steps,
    n_qubits=analysis.n_qubits,
    method=method,
)

print("Operator-norm Trotter error")
print("---------------------------")
print(f"JW error = {operator_result['error_jw']:.16e}")
print(f"BK error = {operator_result['error_bk']:.16e}")
print(f"|diff|   = {operator_result['error_difference']:.16e}")
print()

state_result = bitstring_overlap_comparison(
    psi0=psi0_bitstring,
    C_jw_to_bk=C_jw_to_bk,
    H_jw=analysis.H_jw,
    H_bk=analysis.H_bk,
    jw_terms=analysis.jw_terms,
    bk_terms=analysis.bk_terms,
    time=time,
    num_steps=num_steps,
    n_qubits=analysis.n_qubits,
    method=method,
)

print("State-dependent Trotter overlap")
print("--------------------------------")
print(f"JW amplitude      = {state_result['amp_jw']}")
print(f"BK amplitude      = {state_result['amp_bk']}")
print(f"JW fidelity       = {state_result['fidelity_jw']:.16e}")
print(f"BK fidelity       = {state_result['fidelity_bk']:.16e}")
print(f"JW infidelity     = {state_result['infidelity_jw']:.16e}")
print(f"BK infidelity     = {state_result['infidelity_bk']:.16e}")
print(f"|amplitude diff|  = {state_result['amplitude_difference']:.16e}")
print(f"|fidelity diff|   = {state_result['fidelity_difference']:.16e}")
print(f"|infidelity diff| = {state_result['infidelity_difference']:.16e}")
