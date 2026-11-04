"""Trotter coefficient and sequence helpers local to the comparison package."""

from math import cbrt


def get_trotterization_coefficients(method):
    """Get coefficients for a named ramped Trotterization method."""
    if isinstance(method, str) or isinstance(method, int):
        m = method
        if isinstance(method, str):
            m = method.lower()

        if m in (1, "first order"):
            return (1.0,)
        if m in (2, "second order"):
            return (0.5, 0.5)
        if m in (3, "third order", "ruth 1983"):
            return (7.0 / 24.0, 3.0 / 8.0, 3.0 / 8.0, -25.0 / 24.0, 1.0)
        if m == "symmetrized ruth 1983":
            return (
                7.0 / 48.0,
                3.0 / 16.0,
                3.0 / 16.0,
                -25.0 / 48.0,
                0.5,
                0.5,
                -25.0 / 48.0,
                3.0 / 16.0,
                3.0 / 16.0,
                7.0 / 48.0,
            )
        if m in (4, "fourth order", "suzuki 1990"):
            s2 = 0.5 / (4.0 - cbrt(4.0))
            k = 0.5 - 4.0 * s2
            return (s2, s2, s2, s2, k, k, s2, s2, s2, s2)
        if m in (8, "eighth order", "morales 2022", "morales 2025"):
            b1 = 0.12783360986284110837857554950443
            b2 = 0.56148845266356446893590729572808
            b3 = -0.38400573301491401473462588779099
            b4 = 0.15982762208609923217390166127256
            b5 = -0.40049110428180105319963667975074
            b6 = 0.18669648149540687549831902999911
            b7 = 0.26020394234904150277316667709864
            b8 = 0.29137384767986663096528500968049
            k = 0.5 - (b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8)
            return (
                b1 / 2,
                b1 / 2,
                b2 / 2,
                b2 / 2,
                b3 / 2,
                b3 / 2,
                b4 / 2,
                b4 / 2,
                b5 / 2,
                b5 / 2,
                b6 / 2,
                b6 / 2,
                b7 / 2,
                b7 / 2,
                b8 / 2,
                b8 / 2,
                k,
                k,
                b8 / 2,
                b8 / 2,
                b7 / 2,
                b7 / 2,
                b6 / 2,
                b6 / 2,
                b5 / 2,
                b5 / 2,
                b4 / 2,
                b4 / 2,
                b3 / 2,
                b3 / 2,
                b2 / 2,
                b2 / 2,
                b1 / 2,
                b1 / 2,
            )
        raise ValueError(f"Unknown Trotter method \"{method}\".")

    return tuple(method)


def expand_ramped_trotterization(num_terms, coefficients, num_steps, combine_terms=True):
    """Expand ramped Trotterization into ``(term_index, coefficient)`` pairs."""
    if num_terms == 0:
        raise ValueError("Must have at least one term")
    if num_steps == 0:
        raise ValueError("Must have at least one step")
    if len(coefficients) == 0:
        raise ValueError("Must have at least one coefficient")

    result = []
    ascending = True

    for _ in range(num_steps):
        for coeff in coefficients:
            indices = range(num_terms) if ascending else range(num_terms - 1, -1, -1)
            result.extend((idx, coeff) for idx in indices)
            ascending = not ascending

    if not combine_terms:
        return result

    combined = []
    current_idx, current_coeff = result[0]

    for idx, coeff in result[1:]:
        if idx == current_idx:
            current_coeff += coeff
        else:
            combined.append((current_idx, current_coeff))
            current_idx = idx
            current_coeff = coeff

    combined.append((current_idx, current_coeff))
    return combined
