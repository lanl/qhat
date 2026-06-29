import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import networkx as nx

from qhat.analysis.configuration import load_configuration
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.ordering import create_commutativity_graph


def load_terms_from_config(config_file):
    saved_argv = sys.argv[:]
    try:
        sys.argv = [sys.argv[0], config_file]
        state = load_configuration()
    finally:
        sys.argv = saved_argv

    ham = get_physical_hamiltonian(state.config_hamiltonian)
    return dict(ham.get_all_pauli_strings(return_as="strings"))


def main():
    config_file = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "experiments/h2/config_analysis_none.py"
    )

    terms = load_terms_from_config(config_file)
    term_list = list(terms.items())

    # Despite the function name, QHAT adds edges between NONCOMMUTING Pauli strings.
    G = create_commutativity_graph(term_list, pauli_string_format="explicit")

    coloring = nx.greedy_color(G)

    groups = defaultdict(list)
    for node, color in coloring.items():
        groups[color].append(node)

    print("Number of Pauli terms:", len(term_list))
    print("Number of graph edges:", G.number_of_edges())
    print("Number of greedy color groups:", len(groups))
    print("Group sizes:", [len(groups[c]) for c in sorted(groups)])

    # Store useful node attributes for GraphML / external visualization.
    for i, (pauli, coeff) in enumerate(term_list):
        G.nodes[i]["pauli"] = pauli
        G.nodes[i]["coefficient"] = float(coeff.real if hasattr(coeff, "real") else coeff)
        G.nodes[i]["color_group"] = int(coloring[i])
        G.nodes[i]["label"] = f"{pauli}\n{float(coeff):+.4f}"

    # Make a clean layered layout:
    # each color group gets its own vertical column.
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
        i: f"{pauli}\n{float(coeff):+.4f}"
        for i, (pauli, coeff) in enumerate(term_list)
    }

    node_colors = [coloring[i] for i in G.nodes()]
    node_sizes = [
        1400 + 9000 * abs(float(coeff))
        for _, coeff in term_list
    ]

    plt.figure(figsize=(14, 9))

    nx.draw_networkx_edges(
        G,
        pos,
        width=1.2,
        alpha=0.55,
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=plt.cm.Set3,
        edgecolors="black",
        linewidths=1.0,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        labels=labels,
        font_size=8,
    )

    plt.title(
        "H2/STO-3G JW Pauli Noncommutation Graph\n"
        "nodes = Pauli terms, edges = anticommutation, colors = greedy commuting groups"
    )
    plt.axis("off")
    plt.tight_layout()

    png_path = "experiments/h2/h2_noncommutation_graph.png"
    pdf_path = "experiments/h2/h2_noncommutation_graph.pdf"
    graphml_path = "experiments/h2/h2_noncommutation_graph.graphml"

    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    nx.write_graphml(G, graphml_path)

    print()
    print("Wrote:")
    print(" ", png_path)
    print(" ", pdf_path)
    print(" ", graphml_path)


if __name__ == "__main__":
    main()
