import sys
import csv
from collections import defaultdict

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

try:
    from openfermion import FermionOperator, normal_ordered
except Exception as exc:
    raise SystemExit(
        "Could not import OpenFermion. Try:\n"
        "  python3.11 -m pip install openfermion\n"
    ) from exc


TOL = 1e-10

# For commutation/noncommutation, this factor does not change the graph.
# It only rescales two-body coefficients.
TWO_BODY_FACTOR = 1.0


def term_to_label(term):
    if term == ():
        return "I"

    parts = []
    for mode, action in term:
        if action == 1:
            parts.append(f"a{mode}^")
        else:
            parts.append(f"a{mode}")
    return " ".join(parts)


def is_zero_fermion_operator(op, tol=TOL):
    op = normal_ordered(op)
    return all(abs(c) < tol for c in op.terms.values())


def add_normal_ordered_terms(aggregate, op):
    op = normal_ordered(op)

    for term, coeff in op.terms.items():
        if abs(coeff) < TOL:
            continue
        aggregate[term] += coeff


def load_fermionic_terms(tensor_path):
    data = np.load(tensor_path)

    print("Tensor file:", tensor_path)
    print("NPZ keys:", list(data.keys()))

    if "one_body" not in data or "two_body" not in data:
        raise KeyError(
            "Expected keys 'one_body' and 'two_body' in tensor file."
        )

    one_body = data["one_body"]
    two_body = data["two_body"]

    n_orbitals = one_body.shape[0]

    aggregate = defaultdict(complex)

    # One-body terms: h[p,q] a_p^ a_q
    for p in range(n_orbitals):
        for q in range(n_orbitals):
            coeff = one_body[p, q]
            if abs(coeff) < TOL:
                continue

            op = FermionOperator(((p, 1), (q, 0)), coeff)
            add_normal_ordered_terms(aggregate, op)

    # Two-body terms: h[p,q,r,s] a_p^ a_q^ a_r a_s
    for p in range(n_orbitals):
        for q in range(n_orbitals):
            for r in range(n_orbitals):
                for s in range(n_orbitals):
                    coeff = TWO_BODY_FACTOR * two_body[p, q, r, s]
                    if abs(coeff) < TOL:
                        continue

                    op = FermionOperator(
                        ((p, 1), (q, 1), (r, 0), (s, 0)),
                        coeff,
                    )
                    add_normal_ordered_terms(aggregate, op)

    terms = []
    for term, coeff in aggregate.items():
        if abs(coeff) < TOL:
            continue
        if term == ():
            continue
        terms.append((term, coeff, FermionOperator(term, coeff)))

    return terms


def build_fermionic_noncommutation_graph(terms):
    G = nx.Graph()

    for i, (term, coeff, op) in enumerate(terms):
        G.add_node(
            i,
            label=term_to_label(term),
            coefficient_real=float(np.real(coeff)),
            coefficient_imag=float(np.imag(coeff)),
        )

    for i in range(len(terms)):
        for j in range(i + 1, len(terms)):
            op_i = terms[i][2]
            op_j = terms[j][2]

            comm = normal_ordered(op_i * op_j - op_j * op_i)

            if not is_zero_fermion_operator(comm):
                G.add_edge(i, j)

    return G


def draw_graph(G, coloring, out_png, out_graphml):
    groups = defaultdict(list)
    for node, color in coloring.items():
        groups[color].append(node)

    pos = {}
    x_spacing = 4.0
    y_spacing = 1.2

    for color in sorted(groups):
        nodes = groups[color]
        n = len(nodes)
        for k, node in enumerate(nodes):
            x = color * x_spacing
            y = (n - 1) * y_spacing / 2 - k * y_spacing
            pos[node] = (x, y)

    labels = {
        node: f"{G.nodes[node]['label']}\n{G.nodes[node]['coefficient_real']:+.4f}"
        for node in G.nodes()
    }

    node_colors = [coloring[node] for node in G.nodes()]

    plt.figure(figsize=(16, 10))

    nx.draw_networkx_edges(
        G,
        pos,
        alpha=0.45,
        width=1.0,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        cmap=plt.cm.Set3,
        edgecolors="black",
        linewidths=1.0,
        node_size=1800,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        labels=labels,
        font_size=7,
    )

    plt.title(
        "H2 Fermionic Noncommutation Graph\n"
        "nodes = fermionic monomials, edges = nonzero fermionic commutator, colors = greedy commuting groups"
    )
    plt.axis("off")
    plt.tight_layout()

    plt.savefig(out_png, dpi=300)
    nx.write_graphml(G, out_graphml)


def main():
    tensor_path = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "experiments/h2/h2_0.7414_sto3g_as-002-002.tensors.npz"
    )

    terms = load_fermionic_terms(tensor_path)
    G = build_fermionic_noncommutation_graph(terms)
    coloring = nx.greedy_color(G)

    groups = defaultdict(list)
    for node, color in coloring.items():
        groups[color].append(node)

    print()
    print("Fermionic terms:", len(terms))
    print("Noncommutation edges:", G.number_of_edges())
    print("Greedy color groups:", len(groups))
    print("Group sizes:", [len(groups[c]) for c in sorted(groups)])

    print()
    print("Fermionic commuting groups:")
    for color in sorted(groups):
        print(f"\ncolor {color}:")
        for node in groups[color]:
            label = G.nodes[node]["label"]
            coeff = G.nodes[node]["coefficient_real"]
            print(f"  {coeff:+.12f}  {label}")

    csv_path = "experiments/h2/h2_fermionic_terms_colored.csv"
    png_path = "experiments/h2/h2_fermionic_noncommutation_graph.png"
    graphml_path = "experiments/h2/h2_fermionic_noncommutation_graph.graphml"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "node",
                "color",
                "fermionic_term",
                "coefficient_real",
                "coefficient_imag",
            ]
        )

        for node in G.nodes():
            writer.writerow(
                [
                    node,
                    coloring[node],
                    G.nodes[node]["label"],
                    G.nodes[node]["coefficient_real"],
                    G.nodes[node]["coefficient_imag"],
                ]
            )

    draw_graph(G, coloring, png_path, graphml_path)

    print()
    print("Wrote:")
    print(" ", csv_path)
    print(" ", png_path)
    print(" ", graphml_path)


if __name__ == "__main__":
    main()
