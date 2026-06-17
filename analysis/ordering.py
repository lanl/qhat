# TODO: This function currently expects pauli_strings to be a mapping from
# dense Pauli-string keys to coefficients, e.g. {"XII": 0.5}. Future cleanup
# may add support for generator inputs or a dedicated Pauli term type.
# Some of this stuff is a bit messy and different from the code I originally wrote since my own code had a Pauli Term class
# instead of just a generator/list of raw data. This stuff could possibly be cleaned up by
# My old code is spread across multiple disjoint files, some of which assume a string format and others a tuple format;
# might be good to revisit this decision. Currently, I assume a string format.


def reorder_paulis(pauli_strings, ordering_method):
    """
    Reorder Pauli-string terms for Trotterized time evolution.

    The input is expected to be a mapping from dense Pauli strings to coefficients, for example {"XII": 0.5, "ZZI": -1.2}.
    The returned dictionary contains the same Pauli terms and coefficients, but inserted in the order selected by ordering_method.

    If ordering_method is None, the original input order is preserved.
    """
    pauli_string_list = list(pauli_strings.items())

    # Nothing to reorder; return an empty mapping so callers can still use .items().
    if not pauli_string_list:
        return {}

    # TODO: not a very thorough validation and it's a bit ugly
    if not isinstance(pauli_string_list[0][0], str):
        raise Exception(
            "This method currently only accepts pauli strings written in 'string' format (e.g. XIIYZIX)"
        )

    if ordering_method is None:
        return dict(pauli_string_list)
    elif ordering_method == "magnitude":
        return magnitude(pauli_string_list)
    elif ordering_method == "lexicographical":
        return lexicographical(pauli_string_list)
    elif ordering_method == "group_evolve_xyz":
        return group_evolve_xyz(pauli_string_list)
    else:
        raise Exception(
            f"The Trotter ordering method {ordering_method} is not currently supported."
        )


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
        pauli_types = set(pauli_string)  # throw out duplicates to see which pauli types exist (I, X, Y, Z)
        pauli_types.discard("I") # throw out the identity I if it exists in the string
        if len(pauli_types) > 1:
            raise Exception(
                f"Cannot use this method, group_evolve_xyz can only be used if every pauli term has at most \
                            one non-identity pauli type, but this Hamiltonian has the string {pauli_string}"
            )

        if len(pauli_types) == 0:
            Xs.append(term)
            continue
        else:
            pauli_type = list(pauli_types)[0]  # extract the one pauli type

        if pauli_type == "X":
            Xs.append(term)
        elif pauli_type == "Y":
            Ys.append(term)
        elif pauli_type == "Z":
            Zs.append(term)
        else:
            raise Exception(
                f"Unsupported Pauli type: {pauli_type}. The only allowable Pauli types are I, X, Y, Z."
            )

    return dict(Xs + Ys + Zs)
