"""
Unit tests for symmetry_reduction.py
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from symmetry_reduction import (
    PointGroupTables,
    get_reference_irrep,
    get_excitation_irrep,
    filter_hamiltonian_tensors,
    irrep_id_to_name,
)


class TestPointGroupTables:
    """Test point group product table operations."""

    def test_c1_products(self):
        """C1 has only one irrep, so all products give A."""
        assert PointGroupTables.multiply("C1", "A", "A") == "A"

    def test_cs_products(self):
        """Cs has A' and A" irreps."""
        assert PointGroupTables.multiply("Cs", "A'", "A'") == "A'"
        assert PointGroupTables.multiply("Cs", "A'", 'A"') == 'A"'
        assert PointGroupTables.multiply("Cs", 'A"', "A'") == 'A"'
        assert PointGroupTables.multiply("Cs", 'A"', 'A"') == "A'"

    def test_ci_products(self):
        """Ci has Ag and Au irreps."""
        assert PointGroupTables.multiply("Ci", "Ag", "Ag") == "Ag"
        assert PointGroupTables.multiply("Ci", "Ag", "Au") == "Au"
        assert PointGroupTables.multiply("Ci", "Au", "Ag") == "Au"
        assert PointGroupTables.multiply("Ci", "Au", "Au") == "Ag"

    def test_c2_products(self):
        """C2 has A and B irreps."""
        assert PointGroupTables.multiply("C2", "A", "A") == "A"
        assert PointGroupTables.multiply("C2", "A", "B") == "B"
        assert PointGroupTables.multiply("C2", "B", "A") == "B"
        assert PointGroupTables.multiply("C2", "B", "B") == "A"

    def test_c2v_products(self):
        """C2v has A1, A2, B1, B2 irreps."""
        # Test a few key products
        assert PointGroupTables.multiply("C2v", "A1", "A1") == "A1"
        assert PointGroupTables.multiply("C2v", "A1", "B1") == "B1"
        assert PointGroupTables.multiply("C2v", "B1", "B1") == "A1"
        assert PointGroupTables.multiply("C2v", "B1", "B2") == "A2"
        assert PointGroupTables.multiply("C2v", "A2", "A2") == "A1"

    def test_d2h_products(self):
        """D2h has 8 irreps. Test key products."""
        # Total symmetric
        assert PointGroupTables.multiply("D2h", "Ag", "Ag") == "Ag"

        # g x g = g
        assert PointGroupTables.multiply("D2h", "B1g", "B1g") == "Ag"
        assert PointGroupTables.multiply("D2h", "B1g", "B2g") == "B3g"

        # g x u = u
        assert PointGroupTables.multiply("D2h", "Ag", "Au") == "Au"
        assert PointGroupTables.multiply("D2h", "B1g", "B1u") == "Au"

        # u x u = g
        assert PointGroupTables.multiply("D2h", "Au", "Au") == "Ag"
        assert PointGroupTables.multiply("D2h", "B1u", "B1u") == "Ag"

    def test_non_abelian_mapping(self):
        """Non-Abelian groups should map to Abelian subgroups."""
        assert PointGroupTables.get_abelian_subgroup("Dooh") == "D2h"
        assert PointGroupTables.get_abelian_subgroup("Coov") == "C2v"
        assert PointGroupTables.get_abelian_subgroup("D2h") == "D2h"  # Already Abelian

    def test_unsupported_group(self):
        """Unsupported groups should raise ValueError."""
        with pytest.raises(ValueError):
            PointGroupTables.multiply("XYZ", "A", "A")

    def test_invalid_irrep_pair(self):
        """Invalid irrep combinations should raise ValueError."""
        with pytest.raises(ValueError):
            PointGroupTables.multiply("C2v", "A1", "XYZ")


class TestReferenceIrrep:
    """Test reference state irrep calculation."""

    def test_c2v_reference(self):
        """Test reference irrep for C2v system."""
        # Example: H2O with three occupied orbitals (A1, A1, B1)
        occupied = ["A1", "A1", "B1"]
        ref_irrep = get_reference_irrep(occupied, "C2v")
        # A1 ⊗ A1 = A1, A1 ⊗ B1 = B1
        assert ref_irrep == "B1"

    def test_d2h_reference(self):
        """Test reference irrep for D2h system (should be Ag for closed shell)."""
        # BeH2 example: three occupied (Ag, B1u, Ag)
        occupied = ["Ag", "B1u", "Ag"]
        ref_irrep = get_reference_irrep(occupied, "D2h")
        # Ag ⊗ B1u = B1u, B1u ⊗ Ag = B1u
        assert ref_irrep == "B1u"

    def test_single_orbital(self):
        """Single occupied orbital."""
        occupied = ["Ag"]
        ref_irrep = get_reference_irrep(occupied, "D2h")
        assert ref_irrep == "Ag"

    def test_empty_list(self):
        """Empty list should raise ValueError."""
        with pytest.raises(ValueError):
            get_reference_irrep([], "C2v")


class TestExcitationIrrep:
    """Test excitation operator irrep calculation."""

    def test_single_excitation_d2h(self):
        """Single excitation in D2h."""
        # Example orbitals
        irreps = ["Ag", "B1u", "Ag", "B2u", "B3g"]
        # Excitation 0→3: Ag→B2u
        exc_irrep = get_excitation_irrep(irreps, "D2h", i=0, a=3)
        # B2u ⊗ Ag = B2u
        assert exc_irrep == "B2u"

    def test_double_excitation_d2h(self):
        """Double excitation in D2h."""
        # Example orbitals
        irreps = ["Ag", "B1u", "Ag", "B2u", "B3g"]
        # Excitation 0,1→3,4: Ag,B1u → B2u,B3g
        exc_irrep = get_excitation_irrep(irreps, "D2h", i=0, j=1, a=3, b=4)
        # B2u ⊗ B3g = B1u, B1u ⊗ B1u = Ag, Ag ⊗ Ag = Ag
        # Actually: B2u ⊗ B3g ⊗ B1u ⊗ Ag
        expected = PointGroupTables.multiply("D2h", "B2u", "B3g")
        expected = PointGroupTables.multiply("D2h", expected, "B1u")
        expected = PointGroupTables.multiply("D2h", expected, "Ag")
        assert exc_irrep == expected

    def test_missing_indices(self):
        """Missing indices should raise ValueError."""
        irreps = ["Ag", "B1u"]
        with pytest.raises(ValueError):
            get_excitation_irrep(irreps, "D2h")


class TestHamiltonianFiltering:
    """Test full Hamiltonian tensor filtering."""

    def test_filter_preserves_shape(self):
        """Filtering should preserve tensor shapes."""
        # Create dummy tensors
        constant = 1.0
        one_body = np.random.rand(4, 4)
        two_body = np.random.rand(4, 4, 4, 4)
        mo_irreps = np.array([0, 0, 1, 1])  # Mock irrep IDs

        f_const, f_one, f_two = filter_hamiltonian_tensors(
            constant, one_body, two_body, mo_irreps, "C2", num_occupied=2
        )

        assert f_const == constant
        assert f_one.shape == one_body.shape
        assert f_two.shape == two_body.shape

    def test_filter_zeros_non_preserving(self):
        """Non-symmetry-preserving terms should be zeroed."""
        # Simple C2 case: A (id=0), B (id=1)
        # Reference: A ⊗ A = A
        constant = 1.0
        one_body = np.ones((4, 4))
        two_body = np.ones((4, 4, 4, 4))
        mo_irreps = np.array([0, 0, 1, 1])  # orbitals 0,1 are A; 2,3 are B

        f_const, f_one, f_two = filter_hamiltonian_tensors(
            constant, one_body, two_body, mo_irreps, "C2", num_occupied=2
        )

        # Check that only A-preserving terms survive
        # Excitation 0→2 (A→B): A ⊗ B = B ≠ A, should be zeroed
        assert f_one[2, 0] == 0.0

        # Excitation 0→0 (A→A): A ⊗ A = A, should be kept
        assert f_one[0, 0] == 1.0

    def test_filter_reduction_statistics(self, caplog):
        """Test that reduction statistics are logged."""
        import logging
        caplog.set_level(logging.INFO)

        constant = 1.0
        one_body = np.random.rand(4, 4)
        two_body = np.random.rand(4, 4, 4, 4)
        mo_irreps = np.array([0, 0, 1, 1])

        filter_hamiltonian_tensors(
            constant, one_body, two_body, mo_irreps, "C2", num_occupied=2
        )

        # Check that statistics were logged
        assert "One-body terms:" in caplog.text
        assert "Two-body terms:" in caplog.text
        assert "Total reduction:" in caplog.text


class TestIrrepConversion:
    """Test PySCF irrep ID to name conversion."""

    def test_irrep_id_to_name_d2h(self):
        """Test conversion for D2h group."""
        # PySCF D2h irrep IDs: 0=Ag, 1=B1g, 2=B2g, 3=B3g, 4=Au, 5=B1u, 6=B2u, 7=B3u
        assert irrep_id_to_name("D2h", 0) == "Ag"
        assert irrep_id_to_name("D2h", 5) == "B1u"

    def test_irrep_id_to_name_c2v(self):
        """Test conversion for C2v group."""
        # PySCF C2v irrep IDs: 0=A1, 1=A2, 2=B1, 3=B2
        assert irrep_id_to_name("C2v", 0) == "A1"
        assert irrep_id_to_name("C2v", 2) == "B1"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
