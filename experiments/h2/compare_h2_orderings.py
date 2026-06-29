import os
import sys
import csv
from collections import defaultdict

import networkx as nx
import numpy as np
from scipy.linalg import expm

from qhat.analysis.configuration import load_configuration
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.ordering import (
    reorder_paulis,
    create_commutativity_graph,
)


PAULI = {
    "I": np.array([[1, 0], [0, 1]], dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def load_terms_from_config(config_file):
    saved_argv = sys.argv[:]
    try:
        sys.argv = [sys.argv[0], config_file]
        state = load_configuration()
    finally:
        sys.argv = saved_argv

    ham = get_physical_hamiltonian(state.config_hamiltonian)
    terms = ham.get_all_pauli_strings(return_as="strings")
    return dict(terms)


def pauli_matrix(pauli_string):
    out = np.array([[1]], dtype=complex)
    for p in pauli_string:
        out = np.kron(out, PAULI[p])
    return out


def hamiltonian_matrix(terms):
    n = len(next(iter(terms)))
    dim = 2**n
    H = np.zeros((dim, dim), dtype=complex)

    for p, c in terms.items():
        H += complex(c) * pauli_matrix(p)

    return H


def first_order_trotter(ordered_terms, t=1.0, steps=1):
    n = len(next(iter(ordered_terms)))
    dim = 2**n
    dt = t / steps

    U_step = np.eye(dim, dtype=complex)

    for p, c in ordered_terms.items():
        P = pauli_matrix(p)
        U_step = expm(-1j * dt * complex(c) * P) @ U_step

    U = np.eye(dim, dtype=complex)
    for _ in range(steps):
        U = U_step @ U

    return U


def print_greedy_coloring_stats(terms):
    term_list = list(terms.items())
    G = create_commutativity_graph(term_list, pauli_string_format="explicit")
    coloring = nx.greedy_color(G)

    groups = defaultdict(list)
    for node, color in coloring.items():
        groups[color].append(term_list[node])

    print()
    print("Pauli noncommutation graph")
    print("--------------------------")
    print(f"Pauli terms: {len(term_list)}")
    print(f"Graph edges: {G.number_of_edges()}")
    print(f"Greedy color groups: {len(groups)}")
    print(f"Group sizes: {[len(groups[k]) for k in sorted(groups)]}")

    print()
    print("Greedy color groups")
    print("-------------------")
    for color in sorted(groups):
        print(f"color {color}:")
        for p, c in groups[color]:
            print(f"  {complex(c).real:+.12f}  {p}")


def main():
    config_file = sys.argv[1] if len(sys.argv) > 1 else "experiments/h2/config_analysis_none.py"

    terms = load_terms_from_config(config_file)

    methods = [
        ("none", None),
        ("magnitude", "magnitude"),
        ("lexicographical", "lexicographical"),
        ("group_evolve_greedy", "group_evolve_greedy"),
    ]

    print(f"Loaded from: {config_file}")
    print(f"Number of Pauli terms: {len(terms)}")
    print(f"Number of qubits: {len(next(iter(terms)))}")

    print_greedy_coloring_stats(terms)

    csv_path = "experiments/h2/h2_ordered_pauli_terms.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "position", "pauli_string", "coefficient_real", "coefficient_imag"])

        for label, method in methods:
            ordered = reorder_paulis(terms, method)
            for k, (p, c) in enumerate(ordered.items()):
                c = complex(c)
                writer.writerow([label, k, p, c.real, c.imag])

    print()
    print(f"Wrote ordered Pauli terms to: {csv_path}")

    H = hamiltonian_matrix(terms)
    U_exact = expm(-1j * H)

    print()
    print("First-order Trotter error, t = 1.0")
    print("----------------------------------")
    print("method, steps, spectral_error, frobenius_error")

    for label, method in methods:
        ordered = reorder_paulis(terms, method)

        for steps in [1, 2, 4, 8, 16, 32]:
            U_trot = first_order_trotter(ordered, t=1.0, steps=steps)
            D = U_trot - U_exact
            spectral = np.linalg.norm(D, ord=2)
            frob = np.linalg.norm(D, ord="fro")
            print(f"{label}, {steps}, {spectral:.12e}, {frob:.12e}")


if __name__ == "__main__":
    main()
