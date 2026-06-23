import networkx as nx
from collections import defaultdict

#TODO: The type of pauli_strings is a generator object; is it necessary for the output to also be a generator object,
#  or is an explicit list fine?
# Some of this stuff is a bit messy and different from the code I originally wrote since my own code had a Pauli Term class
# instead of just a generator/list of raw data. This stuff could possibly be cleaned up by 
# My old code is spread across multiple disjoint files, some of which assume a string format and others a tuple format;
# might be good to revisit this decision. Currently, I assume a string format.
def reorder_paulis(pauli_strings, ordering_method):
    # note that this pauli_strings is a generateor, this "exhausts" the generator;
    #   meaning that if it referenced again later, it will produce an empty list
    pauli_string_list = list(pauli_strings.items())

    #TODO: not a very thorough validation and it's a bit ugly
    if not isinstance(pauli_string_list[0][0], str):
        raise Exception("This method currently only accepts pauli strings written in 'string' format (e.g. XIIYZIX)")

    if ordering_method == None:
        return dict(pauli_string_list)
    elif ordering_method == "magnitude":
        return magnitude(pauli_string_list)
    elif ordering_method == "lexicographical":
        return lexicographical(pauli_string_list)
    elif ordering_method == "group_evolve_xyz":
        return group_evolve_xyz(pauli_string_list)
    elif ordering_method == "group_evolve_greedy":
        return group_evolve_greedy(pauli_string_list)
    else:
        raise Exception(f"The Trotter ordering method {ordering_method} is not currently supported.")

def magnitude(terms):
    """
    Sort by descending |coeff|
    """

    def key(t):
        coeff = t[1]
        return abs(coeff)

    return dict(sorted(terms, key=key, reverse=True))

def lexicographical(terms):
    """
    Sort lexicographically by Pauli string in dense format (e.g. XIIYIZ)
    We assume that I < X < Y < Z in the lexicographical ordering
    """

    def key(t):
        pauli_string = t[0]
        return pauli_string

    return dict(sorted(terms, key=key))

def group_evolve_xyz(terms):
    Xs = []
    Ys = []
    Zs = []
    for term in terms:
        pauli_string = term[0]
        pauli_types = set(pauli_string) #throw out duplicates to see which pauli types exist (I, X, Y, Z)
        pauli_types.discard("I") #throw out the identity I if it exists in the string
        if len(pauli_types) > 1:
            raise Exception(f"Cannot use this method, group_evolve_xyz can only be used if every pauli term has at most \
                            one non-identity pauli type, but this Hamiltonian has the string {pauli_string}")

        if len(pauli_types) == 0:
            Xs.append(term)
        else:
            pauli_type = list(pauli_types)[0] #extract the one pauli type
        
        if pauli_type == "X":
            Xs.append(term)
        elif pauli_types == "Y":
            Ys.append(term)
        elif pauli_type == "Z":
            Zs.append(term)
        else:
            raise Exception(f"Unsupported Pauli type: {pauli_type}. The only allowable Pauli types are I, X, Y, Z.")
        
    return dict(Xs + Ys + Zs)

def group_evolve_greedy(pauli_string_list):
    G = create_commutativity_graph(pauli_string_list, pauli_string_format = "explicit")
    coloring = nx.greedy_color(G)

    pauli_string_groups = defaultdict(list)
    num_groups = max(coloring.values())+1
    for index, color in coloring.items():
        pauli_string_groups[color].append(pauli_string_list[index])

    pauli_string_list_reordered = []
    for i in range(num_groups):
        pauli_string_list_reordered += pauli_string_groups[i]

    return dict(pauli_string_list_reordered)



def is_pauli_commuting(P1, P2, pauli_string_format = "dict"):
    if pauli_string_format == "dict":
        anti_count = True #entrywise, how many pauli matrices anticommute
        for i in (P1.keys() & P2.keys()): #identity commutes with everything
            if P1[i] != P2[i]:
                anti_count = not anti_count
    elif pauli_string_format == "explicit":
        if len(P1) != len(P2):
            raise Exception("Lengths of pauli strings (including identities) must be the same!")
        anti_count = True
        for i in range(len(P1)):
            if P1[i] != "I" and P2[i] != "I" and P1[i] != P2[i]:
                anti_count = not anti_count
            
    return anti_count

def create_commutativity_graph(pauli_strings, pauli_string_format = "dict"):
    G = nx.Graph()


    for i in range(len(pauli_strings)):
        G.add_node(i, pauli=pauli_strings[i])

    edges = []
    for i in range(len(pauli_strings)):
        for j in range(i+1, len(pauli_strings)):
            P1 = pauli_strings[i][0] #the [1] extracts the actual string, discards the coeff
            P2 = pauli_strings[j][0]

            if not is_pauli_commuting(P1, P2, pauli_string_format):
                edges.append((i,j))
    
    G.add_edges_from(edges)

    return G

def chromatic_number_greedy(G, strategy="largest_first", interchange=False):
    coloring = nx.greedy_color(G, strategy=strategy, interchange = interchange)
    return max(coloring.values())+1


def chromatic_number_and_coloring_gurobi(G: nx.Graph, time_limit=None, verbose=False):
    """
    Returns (c, coloring) where:
      - c is the chromatic number
      - coloring is a dict: node -> color in {0, ..., c-1}
    """
    # Bounds to speed up the search
    lb = len(nx.approximation.max_clique(G))  # ω(G) lower bound
    greedy = nx.algorithms.coloring.greedy_color(G, strategy="largest_first")
    ub = max(greedy.values()) + 1   # greedy upper bound

    # One global Gurobi environment for params
    env = gp.Env(empty=True)
    if not verbose:
        env.setParam("OutputFlag", 0)
    if time_limit is not None:
        env.setParam("TimeLimit", time_limit)
    env.start()

    def solve_with_K(K: int):
        """Solve model with K available colors. Returns (feasible, coloring_dict_or_None)."""
        m = gp.Model(env=env)
        V = list(G.nodes())

        x = m.addVars(V, range(K), vtype=GRB.BINARY, name="x")  # node-color assignment
        y = m.addVars(range(K), vtype=GRB.BINARY, name="y")     # color used

        # Each node gets exactly one color
        for v in V:
            m.addConstr(gp.quicksum(x[v, k] for k in range(K)) == 1, name=f"onecolor[{v}]")

        # Adjacent nodes cannot share a color
        for u, v in G.edges():
            for k in range(K):
                m.addConstr(x[u, k] + x[v, k] <= 1, name=f"edge[{u},{v},{k}]")

        # Link: if node v uses color k => that color is marked used
        for v in V:
            for k in range(K):
                m.addConstr(x[v, k] <= y[k], name=f"link[{v},{k}]")

        # Symmetry breaking: use colors in order (helps a lot)
        for k in range(K - 1):
            m.addConstr(y[k] >= y[k + 1], name=f"sym[{k}]")

        # Objective: minimize number of used colors within these K
        m.setObjective(gp.quicksum(y[k] for k in range(K)), GRB.MINIMIZE)

        m.optimize()

        if m.Status != GRB.OPTIMAL:
            return False, None

        # Extract chosen color per node
        raw_color = {}
        for v in V:
            chosen = None
            for k in range(K):
                if x[v, k].X > 0.5:
                    chosen = k
                    break
            if chosen is None:
                # Shouldn't happen if optimal, but be safe
                return False, None
            raw_color[v] = chosen

        # Compress colors to 0..c-1 in case some low-probability gap occurs
        # (with symmetry breaking it usually won't, but this guarantees the spec)
        used = sorted(set(raw_color.values()))
        remap = {old: new for new, old in enumerate(used)}
        coloring = {v: remap[k] for v, k in raw_color.items()}

        return True, coloring

    # Binary search on K to find the minimum feasible K
    best_coloring = None
    lo, hi = lb, ub
    while lo <= hi:
        mid = (lo + hi) // 2
        feasible, coloring = solve_with_K(mid)
        if feasible:
            best_coloring = coloring
            hi = mid - 1
        else:
            lo = mid + 1

    if best_coloring is None:
        raise RuntimeError("No coloring found (unexpected).")

    c = max(best_coloring.values()) + 1 if best_coloring else 0
    return c, best_coloring
