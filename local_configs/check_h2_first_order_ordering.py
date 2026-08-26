import pickle
import random

from openfermion.transforms import get_fermion_operator

from qhat.trotter_isomorphism.jw_bk import (
    build_jw_bk_analysis,
    build_clifford_from_z_map,
    verify_graph_isomorphism,
)
from qhat.trotter_isomorphism.comparison import (
    compare_jw_bk_operator_norm_errors,
    bitstring_overlap_comparison,
)
from qhat.trotter_isomorphism.pauli_keys import format_pauli_key


active_pickle = "/tmp/qhat_h2_sto3g/h2_sto3g_as-002-002.pickle"

time = 1.0
num_steps = 1
method = 1
psi0_bitstring = "1100"

random_seed = 12345
num_random_orders = 10

PRINT_ORDER_DETAILS = True
REPORT_PATH = "local_configs/h2_ordering_report.txt"

with open(active_pickle, "rb") as f:
    active_hamiltonian = pickle.load(f)

fermion_op = get_fermion_operator(active_hamiltonian)
n_qubits = active_hamiltonian.n_qubits


def evaluate_order(order_name, jw_order):
    """Evaluate JW/BK first-order Trotter errors for one JW Pauli ordering."""

    analysis = build_jw_bk_analysis(
        fermion_op,
        n_qubits=n_qubits,
        atol=1e-10,
        remove_identity=False,
        build_matrices=True,
        jw_order=jw_order,
    )

    C_jw_to_bk, B, A = build_clifford_from_z_map(
        analysis.jw_to_bk,
        analysis.n_qubits,
    )

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

    return {
        "name": order_name,
        "analysis": analysis,
        "jw_error": operator_result["error_jw"],
        "bk_error": operator_result["error_bk"],
        "error_difference": operator_result["error_difference"],
        "jw_infidelity": state_result["infidelity_jw"],
        "bk_infidelity": state_result["infidelity_bk"],
        "infidelity_difference": state_result["infidelity_difference"],
    }

def format_matched_order(result):
    """Return a readable JW -> BK matched ordering report for one order case."""

    analysis = result["analysis"]

    lines = []
    lines.append("")
    lines.append(f"Matched JW/BK order for {result['name']}")
    lines.append("-" * (27 + len(result["name"])))
    lines.append(
        f"JW operator error = {result['jw_error']:.16e}"
    )
    lines.append(
        f"BK operator error = {result['bk_error']:.16e}"
    )
    lines.append(
        f"|JW-BK|           = {result['error_difference']:.16e}"
    )
    lines.append(
        f"JW infidelity     = {result['jw_infidelity']:.16e}"
    )
    lines.append(
        f"BK infidelity     = {result['bk_infidelity']:.16e}"
    )
    lines.append(
        f"|infid diff|      = {result['infidelity_difference']:.16e}"
    )
    lines.append("")

    lines.append(
        f"{'idx':>3s}  "
        f"{'JW Pauli string':25s} "
        f"{'JW coeff':>18s}   "
        f"{'BK Pauli string':25s} "
        f"{'BK coeff':>18s}"
    )
    lines.append(
        f"{'-'*3}  "
        f"{'-'*25} "
        f"{'-'*18}   "
        f"{'-'*25} "
        f"{'-'*18}"
    )

    for j, ((jw_key, jw_coeff), (bk_key, bk_coeff)) in enumerate(
        zip(analysis.jw_terms, analysis.bk_terms)
    ):
        lines.append(
            f"{j:3d}  "
            f"{format_pauli_key(jw_key):25s} "
            f"{jw_coeff.real:+18.10e}   "
            f"{format_pauli_key(bk_key):25s} "
            f"{bk_coeff.real:+18.10e}"
        )

    return "\n".join(lines)

# Build the default analysis once to get the default JW order.
default_analysis = build_jw_bk_analysis(
    fermion_op,
    n_qubits=n_qubits,
    atol=1e-10,
    remove_identity=False,
    build_matrices=True,
)

default_order = list(default_analysis.jw_order)

# Keep identity fixed at the front.
# The identity only contributes a global phase, so its position should not
# affect the Trotter error. Keeping it fixed makes the experiment cleaner.
identity_terms = [p for p in default_order if p == ()]
nonidentity_terms = [p for p in default_order if p != ()]

orders = []

orders.append(("default_sorted_order", identity_terms + nonidentity_terms))
orders.append(("reversed_nonidentity_order", identity_terms + list(reversed(nonidentity_terms))))

rng = random.Random(random_seed)

for k in range(num_random_orders):
    shuffled = list(nonidentity_terms)
    rng.shuffle(shuffled)
    orders.append((f"random_order_{k:02d}", identity_terms + shuffled))


print("H2/STO-3G first-order Pauli-ordering experiment")
print("==============================================")
print(f"active pickle     = {active_pickle}")
print(f"n_qubits          = {n_qubits}")
print(f"time              = {time}")
print(f"num_steps         = {num_steps}")
print(f"method            = {method}")
print(f"psi0 JW           = |{psi0_bitstring}>")
print(f"random seed       = {random_seed}")
print(f"random order count = {num_random_orders}")
print()

print("Default JW order")
print("----------------")
for j, key in enumerate(default_order):
    print(f"{j:2d}: {format_pauli_key(key)}")
print()

results = []
report_sections = []

for order_name, jw_order in orders:
    result = evaluate_order(order_name, jw_order)
    results.append(result)

    order_report = format_matched_order(result)
    report_sections.append(order_report)

    if PRINT_ORDER_DETAILS:
        print(order_report)

print("Summary")
print("-------")
print(
    f"{'order':30s} "
    f"{'JW op error':>20s} "
    f"{'BK op error':>20s} "
    f"{'|JW-BK|':>12s} "
    f"{'JW infidelity':>20s} "
    f"{'BK infidelity':>20s} "
    f"{'|infid diff|':>14s}"
)

for result in results:
    print(
        f"{result['name']:30s} "
        f"{result['jw_error']:20.16e} "
        f"{result['bk_error']:20.16e} "
        f"{result['error_difference']:12.3e} "
        f"{result['jw_infidelity']:20.16e} "
        f"{result['bk_infidelity']:20.16e} "
        f"{result['infidelity_difference']:14.3e}"
    )

with open(REPORT_PATH, "w") as f:
    f.write("H2/STO-3G first-order Pauli-ordering report\n")
    f.write("===========================================\n")
    f.write(f"active pickle       = {active_pickle}\n")
    f.write(f"n_qubits            = {n_qubits}\n")
    f.write(f"time                = {time}\n")
    f.write(f"num_steps           = {num_steps}\n")
    f.write(f"method              = {method}\n")
    f.write(f"psi0 JW             = |{psi0_bitstring}>\n")
    f.write(f"random seed         = {random_seed}\n")
    f.write(f"random order count  = {num_random_orders}\n")
    f.write("\n")
    f.write("Important note:\n")
    f.write(
        "The Jordan-Wigner qubit/orbital mapping is fixed. "
        "Only the Pauli-string order inside the first-order Trotter product is shuffled.\n"
    )
    f.write(
        "For each shuffled JW order, the BK order is generated by applying the JW->BK map term-by-term.\n"
    )
    f.write("\n")

    for section in report_sections:
        f.write(section)
        f.write("\n\n")

    f.write("Summary\n")
    f.write("-------\n")
    f.write(
        f"{'order':30s} "
        f"{'JW op error':>20s} "
        f"{'BK op error':>20s} "
        f"{'|JW-BK|':>12s} "
        f"{'JW infidelity':>20s} "
        f"{'BK infidelity':>20s} "
        f"{'|infid diff|':>14s}\n"
    )

    for result in results:
        f.write(
            f"{result['name']:30s} "
            f"{result['jw_error']:20.16e} "
            f"{result['bk_error']:20.16e} "
            f"{result['error_difference']:12.3e} "
            f"{result['jw_infidelity']:20.16e} "
            f"{result['bk_infidelity']:20.16e} "
            f"{result['infidelity_difference']:14.3e}\n"
        )

print()
print(f"Full JW/BK ordering report written to: {REPORT_PATH}")