"""
SymUCCSD: Symmetry-reduced Unitary Coupled Cluster Singles and Doubles

This module implements symmetry reduction for Hamiltonian generation based on
molecular point-group symmetry, as described in:
    Cao et al., "Progress toward larger molecular simulation on a quantum computer:
    Simulating a system with up to 28 qubits accelerated by point-group symmetry"
    Physical Review A 105, 062452 (2022)

The key idea: Only excitation operators that preserve the irreducible representation (irrep)
of the reference state contribute to the Hamiltonian. This can reduce operator count by 26-82%
depending on the molecular symmetry.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)


class PointGroupTables:
    """
    Point group character and product tables for Abelian groups.

    Supports: C1, Cs, Ci, C2, C2v, C2h, D2, D2h

    The product tables define how irreps combine under direct product:
        irrep_A ⊗ irrep_B = irrep_C

    This is used to determine the symmetry of excitation operators.
    """

    # Define product tables for supported Abelian groups
    # Each table maps (irrep1, irrep2) -> product_irrep

    TABLES = {
        # C1: Only identity, single irrep "A"
        "C1": {
            ("A", "A"): "A",
        },

        # Cs: One mirror plane
        # Irreps: A' (symmetric), A" (antisymmetric)
        "Cs": {
            ("A'", "A'"): "A'",
            ("A'", 'A"'): 'A"',
            ('A"', "A'"): 'A"',
            ('A"', 'A"'): "A'",
        },

        # Ci: Inversion center
        # Irreps: Ag (symmetric), Au (antisymmetric)
        "Ci": {
            ("Ag", "Ag"): "Ag",
            ("Ag", "Au"): "Au",
            ("Au", "Ag"): "Au",
            ("Au", "Au"): "Ag",
        },

        # C2: Single rotation axis
        # Irreps: A (symmetric), B (antisymmetric)
        "C2": {
            ("A", "A"): "A",
            ("A", "B"): "B",
            ("B", "A"): "B",
            ("B", "B"): "A",
        },

        # C2v: C2 axis + 2 mirror planes
        # Irreps: A1, A2, B1, B2
        "C2v": {
            ("A1", "A1"): "A1", ("A1", "A2"): "A2", ("A1", "B1"): "B1", ("A1", "B2"): "B2",
            ("A2", "A1"): "A2", ("A2", "A2"): "A1", ("A2", "B1"): "B2", ("A2", "B2"): "B1",
            ("B1", "A1"): "B1", ("B1", "A2"): "B2", ("B1", "B1"): "A1", ("B1", "B2"): "A2",
            ("B2", "A1"): "B2", ("B2", "A2"): "B1", ("B2", "B1"): "A2", ("B2", "B2"): "A1",
        },

        # C2h: C2 axis + horizontal mirror plane + inversion
        # Irreps: Ag, Bg, Au, Bu
        "C2h": {
            ("Ag", "Ag"): "Ag", ("Ag", "Bg"): "Bg", ("Ag", "Au"): "Au", ("Ag", "Bu"): "Bu",
            ("Bg", "Ag"): "Bg", ("Bg", "Bg"): "Ag", ("Bg", "Au"): "Bu", ("Bg", "Bu"): "Au",
            ("Au", "Ag"): "Au", ("Au", "Bg"): "Bu", ("Au", "Au"): "Ag", ("Au", "Bu"): "Bg",
            ("Bu", "Ag"): "Bu", ("Bu", "Bg"): "Au", ("Bu", "Au"): "Bg", ("Bu", "Bu"): "Ag",
        },

        # D2: Three perpendicular C2 axes
        # Irreps: A, B1, B2, B3
        "D2": {
            ("A", "A"): "A",   ("A", "B1"): "B1", ("A", "B2"): "B2", ("A", "B3"): "B3",
            ("B1", "A"): "B1", ("B1", "B1"): "A", ("B1", "B2"): "B3", ("B1", "B3"): "B2",
            ("B2", "A"): "B2", ("B2", "B1"): "B3", ("B2", "B2"): "A", ("B2", "B3"): "B1",
            ("B3", "A"): "B3", ("B3", "B1"): "B2", ("B3", "B2"): "B1", ("B3", "B3"): "A",
        },

        # D2h: Three C2 axes + 3 mirror planes + inversion (highest Abelian subgroup of many)
        # Irreps: Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u
        "D2h": {
            # g ⊗ g = g, g ⊗ u = u, u ⊗ g = u, u ⊗ u = g
            ("Ag", "Ag"): "Ag",   ("Ag", "B1g"): "B1g", ("Ag", "B2g"): "B2g", ("Ag", "B3g"): "B3g",
            ("Ag", "Au"): "Au",   ("Ag", "B1u"): "B1u", ("Ag", "B2u"): "B2u", ("Ag", "B3u"): "B3u",

            ("B1g", "Ag"): "B1g", ("B1g", "B1g"): "Ag", ("B1g", "B2g"): "B3g", ("B1g", "B3g"): "B2g",
            ("B1g", "Au"): "B1u", ("B1g", "B1u"): "Au", ("B1g", "B2u"): "B3u", ("B1g", "B3u"): "B2u",

            ("B2g", "Ag"): "B2g", ("B2g", "B1g"): "B3g", ("B2g", "B2g"): "Ag", ("B2g", "B3g"): "B1g",
            ("B2g", "Au"): "B2u", ("B2g", "B1u"): "B3u", ("B2g", "B2u"): "Au", ("B2g", "B3u"): "B1u",

            ("B3g", "Ag"): "B3g", ("B3g", "B1g"): "B2g", ("B3g", "B2g"): "B1g", ("B3g", "B3g"): "Ag",
            ("B3g", "Au"): "B3u", ("B3g", "B1u"): "B2u", ("B3g", "B2u"): "B1u", ("B3g", "B3u"): "Au",

            ("Au", "Ag"): "Au",   ("Au", "B1g"): "B1u", ("Au", "B2g"): "B2u", ("Au", "B3g"): "B3u",
            ("Au", "Au"): "Ag",   ("Au", "B1u"): "B1g", ("Au", "B2u"): "B2g", ("Au", "B3u"): "B3g",

            ("B1u", "Ag"): "B1u", ("B1u", "B1g"): "Au", ("B1u", "B2g"): "B3u", ("B1u", "B3g"): "B2u",
            ("B1u", "Au"): "B1g", ("B1u", "B1u"): "Ag", ("B1u", "B2u"): "B3g", ("B1u", "B3u"): "B2g",

            ("B2u", "Ag"): "B2u", ("B2u", "B1g"): "B3u", ("B2u", "B2g"): "Au", ("B2u", "B3g"): "B1u",
            ("B2u", "Au"): "B2g", ("B2u", "B1u"): "B3g", ("B2u", "B2u"): "Ag", ("B2u", "B3u"): "B1g",

            ("B3u", "Ag"): "B3u", ("B3u", "B1g"): "B2u", ("B3u", "B2g"): "B1u", ("B3u", "B3g"): "Au",
            ("B3u", "Au"): "B3g", ("B3u", "B1u"): "B2g", ("B3u", "B2u"): "B1g", ("B3u", "B3u"): "Ag",
        },

        # Non-Abelian groups - map to their highest Abelian subgroup
        # Users can extend this for full non-Abelian support
        "D∞h": None,  # Use D2h as approximation
        "Dooh": None,  # Use D2h as approximation
        "C∞v": None,  # Use C2v as approximation
        "Coov": None,  # Use C2v as approximation
    }

    # Irrep name mappings for non-Abelian groups to their Abelian subgroups
    IRREP_MAPPINGS = {
        # Dooh (D∞h) → D2h mapping
        # PySCF uses A1g, A2g, E1gx/E1gy, A2u, A1u, E1ux/E1uy for Dooh
        # Map to D2h: Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u
        "Dooh": {
            "A1g": "Ag", "A2g": "B1g",
            "E1gx": "B2g", "E1gy": "B3g",
            "A2u": "B1u", "A1u": "Au",
            "E1ux": "B2u", "E1uy": "B3u",
        },
        # Coov (C∞v) → C2v mapping
        "Coov": {
            "A1": "A1", "A2": "A2",
            "E1x": "B1", "E1y": "B2",
        },
    }

    @classmethod
    def get_abelian_subgroup(cls, point_group: str) -> str:
        """Map non-Abelian groups to their highest Abelian subgroup."""
        mapping = {
            "D∞h": "D2h", "Dooh": "D2h",
            "C∞v": "C2v", "Coov": "C2v",
            "D3": "C2", "D3h": "C2h", "D3d": "C2h",
            "D4": "D2", "D4h": "D2h",
            "D6": "D2", "D6h": "D2h",
            "T": "D2", "Td": "D2", "Th": "D2h",
            "O": "D2", "Oh": "D2h",
        }
        return mapping.get(point_group, point_group)

    @classmethod
    def map_irrep_name(cls, point_group: str, irrep_name: str) -> str:
        """Map irrep name from non-Abelian group to Abelian subgroup."""
        if point_group in cls.IRREP_MAPPINGS:
            mapping = cls.IRREP_MAPPINGS[point_group]
            return mapping.get(irrep_name, irrep_name)
        return irrep_name

    @classmethod
    def multiply(cls, point_group: str, irrep1: str, irrep2: str) -> str:
        """
        Compute the direct product of two irreps in the given point group.

        Args:
            point_group: Point group name (e.g., "D2h", "C2v", "Dooh")
            irrep1: First irrep label
            irrep2: Second irrep label

        Returns:
            Product irrep label
        """
        # Map to Abelian subgroup if needed
        original_pg = point_group
        pg = cls.get_abelian_subgroup(point_group)

        # Map irrep names if needed
        mapped_irrep1 = cls.map_irrep_name(original_pg, irrep1)
        mapped_irrep2 = cls.map_irrep_name(original_pg, irrep2)

        if pg not in cls.TABLES:
            raise ValueError(f"Unsupported point group: {point_group} (mapped to {pg})")

        table = cls.TABLES[pg]
        key = (mapped_irrep1, mapped_irrep2)

        if key not in table:
            raise ValueError(f"Invalid irrep pair for {pg}: {mapped_irrep1} ⊗ {mapped_irrep2} " +
                           f"(original: {irrep1} ⊗ {irrep2})")

        result = table[key]

        # Map result back to original group notation if needed
        # (For now, keep in Abelian notation for consistency)
        return result

    @classmethod
    def supported_groups(cls) -> List[str]:
        """Return list of supported point groups."""
        return [pg for pg, table in cls.TABLES.items() if table is not None]


def irrep_id_to_name(point_group: str, irrep_id: int) -> str:
    """
    Convert PySCF irrep ID to irrep name for a given point group.

    PySCF stores irreps as integer IDs. This function converts them to
    standard Mulliken notation.

    Args:
        point_group: Point group name from PySCF
        irrep_id: Integer irrep ID from PySCF

    Returns:
        Irrep name (e.g., "Ag", "B1u")
    """
    try:
        import pyscf.symm
        return pyscf.symm.irrep_id2name(point_group, irrep_id)
    except Exception as e:
        logger.warning(f"Could not convert irrep ID {irrep_id} for group {point_group}: {e}")
        return f"irrep_{irrep_id}"


def get_reference_irrep(occupied_mo_irreps: List[str], point_group: str) -> str:
    """
    Compute the irrep of the reference (Hartree-Fock) state.

    The reference state is a Slater determinant of occupied orbitals.
    Its irrep is the direct product of all occupied orbital irreps.

    Args:
        occupied_mo_irreps: List of irrep names for occupied MOs
        point_group: Point group name

    Returns:
        Irrep name of the reference state
    """
    if not occupied_mo_irreps:
        raise ValueError("No occupied orbitals provided")

    # Start with the first irrep
    result = occupied_mo_irreps[0]

    # Multiply with each subsequent occupied orbital irrep
    for irrep in occupied_mo_irreps[1:]:
        result = PointGroupTables.multiply(point_group, result, irrep)

    logger.debug(f"Reference state irrep: {result} (from {occupied_mo_irreps})")
    return result


def get_excitation_irrep(
    orbital_irreps: List[str],
    point_group: str,
    i: Optional[int] = None,
    a: Optional[int] = None,
    j: Optional[int] = None,
    b: Optional[int] = None
) -> str:
    """
    Compute the irrep of an excitation operator.

    For single excitation i→a: irrep(a) ⊗ irrep(i)
    For double excitation i,j→a,b: irrep(a) ⊗ irrep(b) ⊗ irrep(j) ⊗ irrep(i)

    Note: Order matters for proper convention matching

    Args:
        orbital_irreps: List of irrep names for all MOs
        point_group: Point group name
        i, j: Occupied orbital indices (creation part)
        a, b: Virtual orbital indices (annihilation part)

    Returns:
        Irrep name of the excitation
    """
    if a is not None and i is not None:
        if j is not None and b is not None:
            # Double excitation: a†b†ji
            exc_irrep = PointGroupTables.multiply(point_group, orbital_irreps[a], orbital_irreps[b])
            exc_irrep = PointGroupTables.multiply(point_group, exc_irrep, orbital_irreps[j])
            exc_irrep = PointGroupTables.multiply(point_group, exc_irrep, orbital_irreps[i])
        else:
            # Single excitation: a†i
            exc_irrep = PointGroupTables.multiply(point_group, orbital_irreps[a], orbital_irreps[i])
        return exc_irrep
    else:
        raise ValueError("Must specify at least i and a for excitation")


def filter_hamiltonian_tensors(
    constant: float,
    one_body: np.ndarray,
    two_body: np.ndarray,
    mo_irreps: np.ndarray,
    point_group: str,
    num_occupied: int
) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Filter Hamiltonian tensors to keep only symmetry-preserving terms.

    This is the core symmetry reduction operation. For each tensor element,
    we check if the corresponding excitation preserves the reference irrep.
    If not, we zero it out.

    Args:
        constant: Zero-body term (nuclear repulsion + core energy)
        one_body: One-body tensor h[p,q] (spatial orbitals)
        two_body: Two-body tensor g[p,q,r,s] (spatial orbitals, physicist convention)
        mo_irreps: Array of irrep IDs for each MO (from PySCF)
        point_group: Point group name
        num_occupied: Number of occupied spatial orbitals

    Returns:
        Tuple of (filtered_constant, filtered_one_body, filtered_two_body)
    """
    # Convert irrep IDs to names
    irrep_names = [irrep_id_to_name(point_group, iid) for iid in mo_irreps]

    # Get reference irrep from occupied orbitals (doubled for spin)
    occupied_irreps = irrep_names[:num_occupied]
    ref_irrep = get_reference_irrep(occupied_irreps, point_group)

    logger.info(f"Reference state irrep: {ref_irrep}")
    logger.info(f"Filtering Hamiltonian with {point_group} symmetry")

    # Constant term always preserved
    filtered_constant = constant

    # Filter one-body tensor
    filtered_one_body = one_body.copy()
    one_body_kept = 0
    one_body_total = 0

    for p in range(one_body.shape[0]):
        for q in range(one_body.shape[1]):
            if abs(one_body[p, q]) < 1e-12:
                continue  # Already zero

            one_body_total += 1

            # One-body term h[p,q] corresponds to excitation q→p (a†_p a_q)
            try:
                exc_irrep = get_excitation_irrep(irrep_names, point_group, i=q, a=p)
                if exc_irrep == ref_irrep:
                    one_body_kept += 1
                else:
                    filtered_one_body[p, q] = 0.0
            except Exception as e:
                logger.debug(f"Could not determine irrep for h[{p},{q}]: {e}")
                # Keep term if we can't determine (conservative)
                one_body_kept += 1

    # Filter two-body tensor
    filtered_two_body = two_body.copy()
    two_body_kept = 0
    two_body_total = 0

    for p in range(two_body.shape[0]):
        for q in range(two_body.shape[1]):
            for r in range(two_body.shape[2]):
                for s in range(two_body.shape[3]):
                    if abs(two_body[p, q, r, s]) < 1e-12:
                        continue  # Already zero

                    two_body_total += 1

                    # Two-body term g[p,q,r,s] corresponds to a†_p a†_q a_s a_r
                    try:
                        exc_irrep = get_excitation_irrep(
                            irrep_names, point_group, i=r, a=p, j=s, b=q
                        )
                        if exc_irrep == ref_irrep:
                            two_body_kept += 1
                        else:
                            filtered_two_body[p, q, r, s] = 0.0
                    except Exception as e:
                        logger.debug(f"Could not determine irrep for g[{p},{q},{r},{s}]: {e}")
                        # Keep term if we can't determine (conservative)
                        two_body_kept += 1

    # Report statistics
    if one_body_total > 0:
        one_body_reduction = 100 * (1 - one_body_kept / one_body_total)
        logger.info(f"One-body terms: kept {one_body_kept}/{one_body_total} " +
                   f"({one_body_reduction:.1f}% reduction)")

    if two_body_total > 0:
        two_body_reduction = 100 * (1 - two_body_kept / two_body_total)
        logger.info(f"Two-body terms: kept {two_body_kept}/{two_body_total} " +
                   f"({two_body_reduction:.1f}% reduction)")

    total_kept = one_body_kept + two_body_kept
    total = one_body_total + two_body_total
    if total > 0:
        total_reduction = 100 * (1 - total_kept / total)
        logger.info(f"Total reduction: {total_reduction:.1f}% " +
                   f"({total_kept}/{total} terms kept)")

    return filtered_constant, filtered_one_body, filtered_two_body
