import pickle

from openfermion.transforms import get_fermion_operator

from qhat.trotter_isomorphism.jw_bk import (
    build_jw_bk_analysis,
    verify_graph_isomorphism,
)
from qhat.trotter_isomorphism.majorana import (
    fermion_operator_to_majorana_operator_dict,
)
from qhat.trotter_isomorphism.pauli_keys import format_pauli_key


active_pickle = "/tmp/qhat_h_atom_smoke/h_atom_as-001-001.pickle"

with open(active_pickle, "rb") as f:
    active_hamiltonian = pickle.load(f)

fermion_op = get_fermion_operator(active_hamiltonian)

print("Hydrogen atom Majorana JW/BK smoke test")
print("======================================")
print(f"active pickle = {active_pickle}")
print(f"n_qubits      = {active_hamiltonian.n_qubits}")
print()

print("FermionOperator")
print("---------------")
print(fermion_op)
print()

majorana_terms = fermion_operator_to_majorana_operator_dict(
    fermion_op,
    atol=1e-10,
)

print("Majorana expansion")
print("------------------")
print(f"number of Majorana monomials = {len(majorana_terms)}")
for label, coeff in sorted(majorana_terms.items()):
    print(f"gamma{label}: {coeff}")
print()

analysis = build_jw_bk_analysis(
    fermion_op,
    n_qubits=active_hamiltonian.n_qubits,
    atol=1e-10,
    remove_identity=False,
    build_matrices=True,
)

print("Graph-isomorphism check")
print("----------------------")
print(f"JW vertices = {len(analysis.jw_to_bk)}")
print(f"BK vertices = {len(analysis.bk_to_jw)}")
print(f"verified    = {verify_graph_isomorphism(analysis.jw_to_bk)}")
print()

print("Majorana-induced JW -> BK Pauli map")
print("-----------------------------------")
for jw_key in sorted(analysis.jw_to_bk, key=format_pauli_key):
    bk_key = analysis.jw_to_bk[jw_key]
    print(
        f"{format_pauli_key(jw_key):10s} -> {format_pauli_key(bk_key)}"
    )
print()

print("Matched JW/BK terms used for Trotter ordering")
print("---------------------------------------------")
for (jw_key, jw_coeff), (bk_key, bk_coeff) in zip(
    analysis.jw_terms,
    analysis.bk_terms,
):
    print(
        f"JW {format_pauli_key(jw_key):10s} coeff={jw_coeff.real:+.16e}"
        f"   -->   BK {format_pauli_key(bk_key):10s} coeff={bk_coeff.real:+.16e}"
    )
print()

print("Dense-matrix sanity checks")
print("--------------------------")
print(f"H_jw shape = {analysis.H_jw.shape}")
print(f"H_bk shape = {analysis.H_bk.shape}")
print("done")
