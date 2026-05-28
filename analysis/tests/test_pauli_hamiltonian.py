"""
Tests for Pauli string Hamiltonian loading and LinearCombinationOfPauliStrings class.

Tests cover:
- LinearCombinationOfPauliStrings class functionality
- Utility functions (dense/sparse conversion)
- File loading (dense and sparse formats)
- Hermitian validation
- Integration with Hamiltonian API
"""

import os
import tempfile
import pytest
import numpy as np

from qhat.analysis.config_types import (
    GeneralConfiguration,
    GeneralConfigurationUser,
    HamiltonianConfiguration,
)
from qhat.analysis.hamiltonian import (
    Hamiltonian,
    LinearCombinationOfPauliStrings,
    dense_to_sparse_pauli,
    load_pauli,
    sparse_to_dense_pauli,
)


# ==================================================================================
# Shared fixtures
# ==================================================================================

@pytest.fixture(scope="session")
def general_config():
    """Create a single GeneralConfiguration for the entire test session."""
    user_config = GeneralConfigurationUser()
    return GeneralConfiguration(user_config)


# ==================================================================================
# Test: Utility functions
# ==================================================================================

class TestUtilityFunctions:
    """Test sparse_to_dense_pauli and dense_to_sparse_pauli conversions."""

    def test_sparse_to_dense_identity(self):
        """Test converting identity (empty sparse) to dense."""
        assert sparse_to_dense_pauli(tuple(), 4) == "IIII"
        assert sparse_to_dense_pauli(tuple(), 1) == "I"

    def test_sparse_to_dense_single_op(self):
        """Test converting single Pauli operator."""
        assert sparse_to_dense_pauli(((0, 'X'),), 4) == "XIII"
        assert sparse_to_dense_pauli(((2, 'Z'),), 4) == "IIZI"
        assert sparse_to_dense_pauli(((3, 'Y'),), 5) == "IIIYI"

    def test_sparse_to_dense_multi_op(self):
        """Test converting multiple Pauli operators."""
        assert sparse_to_dense_pauli(((0, 'X'), (2, 'Z')), 4) == "XIZI"
        assert sparse_to_dense_pauli(((1, 'Y'), (3, 'Z')), 5) == "IYIZI"

    def test_dense_to_sparse_identity(self):
        """Test converting dense identity to sparse."""
        assert dense_to_sparse_pauli("IIII") == tuple()
        assert dense_to_sparse_pauli("I") == tuple()

    def test_dense_to_sparse_single_op(self):
        """Test converting single operator from dense."""
        assert dense_to_sparse_pauli("XIII") == ((0, 'X'),)
        assert dense_to_sparse_pauli("IIZI") == ((2, 'Z'),)

    def test_dense_to_sparse_multi_op(self):
        """Test converting multiple operators from dense."""
        assert dense_to_sparse_pauli("XIZI") == ((0, 'X'), (2, 'Z'))
        assert dense_to_sparse_pauli("XYZ") == ((0, 'X'), (1, 'Y'), (2, 'Z'))

    def test_dense_to_sparse_invalid_char(self):
        """Test that invalid characters raise ValueError."""
        with pytest.raises(ValueError, match="Invalid character"):
            dense_to_sparse_pauli("XQZI")

    def test_round_trip_conversion(self):
        """Test that converting back and forth preserves the Pauli string."""
        test_cases = [
            ("IIII", 4),
            ("XIII", 4),
            ("XYZI", 4),
            ("XYZIXYZ", 7),
        ]
        for dense, nq in test_cases:
            sparse = dense_to_sparse_pauli(dense)
            dense_back = sparse_to_dense_pauli(sparse, nq)
            assert dense == dense_back


# ==================================================================================
# Test: LinearCombinationOfPauliStrings class
# ==================================================================================

class TestLinearCombinationOfPauliStrings:
    """Test the LinearCombinationOfPauliStrings class."""

    def test_dense_format_init(self):
        """Test initialization with dense format."""
        data = {"IIII": 1.0, "XIII": 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, dense=data)
        assert lcps.num_qubits() == 4
        assert lcps._format == "dense"
        assert len(lcps._data) == 2

    def test_sparse_format_init(self):
        """Test initialization with sparse format."""
        data = {tuple(): 1.0, ((0, 'X'),): 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, sparse=data)
        assert lcps.num_qubits() == 4
        assert lcps._format == "sparse"
        assert len(lcps._data) == 2

    def test_no_data_raises_error(self):
        """Test that initialization without data raises ValueError."""
        with pytest.raises(ValueError, match="No data provided"):
            LinearCombinationOfPauliStrings(num_qubits=4)

    def test_multiple_formats_raises_error(self):
        """Test that providing both formats raises ValueError."""
        with pytest.raises(ValueError, match="Too many formats"):
            LinearCombinationOfPauliStrings(
                num_qubits=4,
                dense={"IIII": 1.0},
                sparse={tuple(): 1.0}
            )

    def test_get_dense_from_dense(self):
        """Test getting dense format when stored as dense."""
        data = {"IIII": 1.0, "XIII": 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, dense=data)
        result = lcps.get_dense_pauli_strings()
        assert result == data

    def test_get_sparse_from_sparse(self):
        """Test getting sparse format when stored as sparse."""
        data = {tuple(): 1.0, ((0, 'X'),): 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, sparse=data)
        result = lcps.get_sparse_pauli_strings()
        assert result == data

    def test_get_sparse_from_dense(self):
        """Test converting dense to sparse."""
        dense_data = {"IIII": 1.0, "XIII": 0.5, "IXZI": 0.3}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, dense=dense_data)
        sparse_result = lcps.get_sparse_pauli_strings()

        expected = {
            tuple(): 1.0,
            ((0, 'X'),): 0.5,
            ((1, 'X'), (2, 'Z')): 0.3
        }
        assert sparse_result == expected

    def test_get_dense_from_sparse(self):
        """Test converting sparse to dense."""
        sparse_data = {tuple(): 1.0, ((0, 'X'),): 0.5, ((1, 'X'), (2, 'Z')): 0.3}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, sparse=sparse_data)
        dense_result = lcps.get_dense_pauli_strings()

        expected = {"IIII": 1.0, "XIII": 0.5, "IXZI": 0.3}
        assert dense_result == expected

    def test_energy_shift_dense(self):
        """Test energy shift with dense format."""
        data = {"IIII": 1.0, "XIII": 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, dense=data)
        lcps.energy_shift(2.5)
        assert lcps._data["IIII"] == 3.5
        assert lcps._data["XIII"] == 0.5

    def test_energy_shift_sparse(self):
        """Test energy shift with sparse format."""
        data = {tuple(): 1.0, ((0, 'X'),): 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, sparse=data)
        lcps.energy_shift(2.5)
        assert lcps._data[tuple()] == 3.5
        assert lcps._data[((0, 'X'),)] == 0.5

    def test_energy_shift_adds_identity(self):
        """Test that energy shift adds identity term if not present."""
        data = {"XIII": 0.5}
        lcps = LinearCombinationOfPauliStrings(num_qubits=4, dense=data)
        lcps.energy_shift(2.5)
        assert lcps._data["IIII"] == 2.5
        assert lcps._data["XIII"] == 0.5


# ==================================================================================
# Test: File loading
# ==================================================================================

class TestFileLoading:
    """Test loading Pauli string Hamiltonians from files."""

    @pytest.fixture
    def config(self, general_config):
        """Create configuration objects for testing."""
        config_ham = HamiltonianConfiguration()
        return general_config, config_ham

    def test_load_dense_format(self, config, tmp_path):
        """Test loading dense format Pauli file."""
        # Create test file
        test_file = tmp_path / "test_dense.txt"
        test_file.write_text("""# Test Hamiltonian
2.5 IIII
0.5 XIII
-0.3 IXII
1.2 XXII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 4
        pauli_dict = H.get_all_pauli_strings(return_as="strings")
        assert len(pauli_dict) == 4
        assert pauli_dict["IIII"] == 2.5
        assert pauli_dict["XIII"] == 0.5
        assert pauli_dict["IXII"] == -0.3
        assert pauli_dict["XXII"] == 1.2

    def test_load_sparse_format(self, config, tmp_path):
        """Test loading sparse format Pauli file."""
        # Create test file
        test_file = tmp_path / "test_sparse.dat"
        test_file.write_text("""# Test Hamiltonian
(2.5) []+
(0.5) [X0]+
(-0.3) [X1]+
(1.2) [X0 X1]+
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 2
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        assert len(pauli_dict) == 4
        assert pauli_dict[tuple()] == 2.5
        assert pauli_dict[((0, 'X'),)] == 0.5
        assert pauli_dict[((1, 'X'),)] == -0.3
        assert pauli_dict[((0, 'X'), (1, 'X'))] == 1.2

    def test_load_empty_lines_and_comments(self, config, tmp_path):
        """Test that empty lines and comments are handled correctly."""
        test_file = tmp_path / "test_comments.txt"
        test_file.write_text("""# This is a comment
1.0 IIII

# Another comment

0.5 XIII
# Yet another comment
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        pauli_dict = H.get_all_pauli_strings(return_as="strings")
        assert len(pauli_dict) == 2

    def test_inconsistent_qubit_count_dense(self, config, tmp_path):
        """Test that inconsistent qubit counts raise error in dense format."""
        test_file = tmp_path / "test_inconsistent.txt"
        test_file.write_text("""1.0 IIII
0.5 XIII
0.3 IXI
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="Inconsistent dense Pauli string length"):
            load_pauli(config_ham)

    def test_mixed_format_raises_error(self, config, tmp_path):
        """Test that mixing dense and sparse formats raises error."""
        test_file = tmp_path / "test_mixed.txt"
        test_file.write_text("""1.0 IIII
(0.5) [X0]+
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="Inconsistent Pauli string file format"):
            load_pauli(config_ham)

    def test_load_json_format(self, config, tmp_path):
        """Test loading JSON format Pauli file."""
        import json
        test_file = tmp_path / "test.json"
        data = {
            "n_qubits": 4,
            "terms": [
                {"ops": [], "coeff": 2.5},
                {"ops": [[0, "X"]], "coeff": 0.5},
                {"ops": [[1, "X"]], "coeff": -0.3},
                {"ops": [[0, "X"], [1, "X"]], "coeff": 1.2}
            ]
        }
        test_file.write_text(json.dumps(data))

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 4
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        assert len(pauli_dict) == 4
        assert pauli_dict[tuple()] == 2.5
        assert pauli_dict[((0, 'X'),)] == 0.5
        assert pauli_dict[((1, 'X'),)] == -0.3
        assert pauli_dict[((0, 'X'), (1, 'X'))] == 1.2

    def test_json_missing_n_qubits(self, config, tmp_path):
        """Test that JSON without n_qubits raises error."""
        import json
        test_file = tmp_path / "test.json"
        data = {"terms": [{"ops": [], "coeff": 1.0}]}
        test_file.write_text(json.dumps(data))

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="must contain 'n_qubits'"):
            load_pauli(config_ham)

    def test_json_missing_terms(self, config, tmp_path):
        """Test that JSON without terms raises error."""
        import json
        test_file = tmp_path / "test.json"
        data = {"n_qubits": 4}
        test_file.write_text(json.dumps(data))

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="must contain 'terms'"):
            load_pauli(config_ham)

    def test_json_term_missing_fields(self, config, tmp_path):
        """Test that JSON terms must have ops and coeff."""
        import json
        test_file = tmp_path / "test.json"
        data = {
            "n_qubits": 4,
            "terms": [{"ops": [[0, "X"]]}]  # missing coeff
        }
        test_file.write_text(json.dumps(data))

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="must have 'ops' and 'coeff'"):
            load_pauli(config_ham)


# ==================================================================================
# Test: Hermitian validation
# ==================================================================================

class TestHermitianValidation:
    """Test validation that Hamiltonian coefficients are real."""

    @pytest.fixture
    def config(self, general_config):
        """Create configuration objects for testing."""
        config_ham = HamiltonianConfiguration()
        return general_config, config_ham

    def test_real_coefficients_pass(self, config, tmp_path):
        """Test that real coefficients pass validation."""
        test_file = tmp_path / "test_real.txt"
        test_file.write_text("""1.0 IIII
0.5 XIII
-0.3 IXII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        # Should not raise
        assert H.num_qubits() == 4

    def test_small_imaginary_passes(self, config, tmp_path):
        """Test that negligible imaginary parts pass (< 1e-8 relative)."""
        test_file = tmp_path / "test_small_imag.txt"
        # 1e-9 relative imaginary part should pass
        test_file.write_text("""(1.0+1e-9j) IIII
(0.5+5e-10j) XIII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        # Should not raise
        assert H.num_qubits() == 4

    def test_large_imaginary_fails(self, config, tmp_path):
        """Test that significant imaginary parts fail validation."""
        test_file = tmp_path / "test_large_imag.txt"
        test_file.write_text("""1.0 IIII
(0.5+0.1j) XIII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="Hamiltonian must be Hermitian"):
            load_pauli(config_ham)

    def test_relative_tolerance_small_coefficients(self, config, tmp_path):
        """Test that relative tolerance works for small coefficients."""
        test_file = tmp_path / "test_small_coef.txt"
        # 1e-12 + 1e-11j has 90% imaginary component, should fail
        test_file.write_text("""(1e-12+1e-11j) IIII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="Hamiltonian must be Hermitian"):
            load_pauli(config_ham)

    def test_converted_to_float(self, config, tmp_path):
        """Test that validated coefficients are converted to float."""
        test_file = tmp_path / "test_float.txt"
        test_file.write_text("""1.0 IIII
0.5 XIII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        pauli_dict = H.get_all_pauli_strings(return_as="strings")
        # Check that values are float, not complex
        for val in pauli_dict.values():
            assert isinstance(val, float)


# ==================================================================================
# Test: Hamiltonian integration
# ==================================================================================

class TestHamiltonianIntegration:
    """Test integration with Hamiltonian class and analysis pipeline."""

    @pytest.fixture
    def config(self, general_config):
        """Create configuration objects for testing."""
        config_ham = HamiltonianConfiguration()
        return general_config, config_ham

    @pytest.fixture
    def sample_hamiltonian(self, config, tmp_path):
        """Create a sample Hamiltonian for testing."""
        test_file = tmp_path / "test_ham.txt"
        test_file.write_text("""2.5 IIII
0.5 XIII
-0.3 IXII
0.7 IIZI
1.2 XXII
-0.8 YYII
0.6 ZZII
0.4 IXZI
0.25 XYZI
-0.15 ZYXI
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        return load_pauli(config_ham), config_gen, config_ham

    def test_get_all_pauli_strings_tuples(self, sample_hamiltonian):
        """Test getting Pauli strings as tuples."""
        H, _, _ = sample_hamiltonian
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")

        assert len(pauli_dict) == 10
        assert pauli_dict[tuple()] == 2.5
        assert pauli_dict[((0, 'X'),)] == 0.5

    def test_get_all_pauli_strings_strings(self, sample_hamiltonian):
        """Test getting Pauli strings as strings."""
        H, _, _ = sample_hamiltonian
        pauli_dict = H.get_all_pauli_strings(return_as="strings")

        assert len(pauli_dict) == 10
        assert pauli_dict["IIII"] == 2.5
        assert pauli_dict["XIII"] == 0.5

    def test_energy_bounds(self, sample_hamiltonian):
        """Test energy bounds computation."""
        H, config_gen, config_ham = sample_hamiltonian
        bounds = H.compute_initial_energy_bounds(config_ham)

        # Bounds should be [E_0 - sum|c_i|, E_0 + sum|c_i|]
        # where E_0 is the identity coefficient
        assert len(bounds) == 2
        assert bounds[0] < bounds[1]

        # Verify calculation
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        identity = pauli_dict[tuple()]
        total_abs = sum(abs(c) for c in pauli_dict.values())
        dE = total_abs - abs(identity)

        assert abs(bounds[0] - (identity - dE)) < 1e-10
        assert abs(bounds[1] - (identity + dE)) < 1e-10

    def test_energy_shift(self, sample_hamiltonian):
        """Test energy shift functionality."""
        H, _, _ = sample_hamiltonian

        pauli_before = H.get_all_pauli_strings(return_as="tuples")
        identity_before = pauli_before[tuple()]

        H.energy_shift(10.0)

        pauli_after = H.get_all_pauli_strings(return_as="tuples")
        identity_after = pauli_after[tuple()]

        assert abs(identity_after - identity_before - 10.0) < 1e-10

    def test_get_grouped_terms(self, sample_hamiltonian):
        """Test getting grouped terms as QubitOperators."""
        H, _, _ = sample_hamiltonian
        grouped = H.get_grouped_terms()

        # Should return list of QubitOperator objects
        assert len(grouped) == 10
        from openfermion import QubitOperator
        for term in grouped:
            assert isinstance(term, QubitOperator)

    def test_get_core_operator(self, sample_hamiltonian):
        """Test getting core operator."""
        H, _, _ = sample_hamiltonian
        core = H.get_core_operator()

        assert isinstance(core, LinearCombinationOfPauliStrings)
        assert core.num_qubits() == 4


# ==================================================================================
# Test: Edge cases
# ==================================================================================

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.fixture
    def config(self, general_config):
        """Create configuration objects for testing."""
        config_ham = HamiltonianConfiguration()
        return general_config, config_ham

    def test_single_qubit(self, config, tmp_path):
        """Test with single qubit system."""
        test_file = tmp_path / "test_1q.txt"
        test_file.write_text("""1.0 I
0.5 X
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 1

    def test_large_qubit_count(self, config, tmp_path):
        """Test with larger qubit systems."""
        test_file = tmp_path / "test_10q.txt"
        test_file.write_text(f"""1.0 {"I" * 10}
0.5 {"X" + "I" * 9}
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 10

    def test_only_identity(self, config, tmp_path):
        """Test Hamiltonian with only identity term."""
        test_file = tmp_path / "test_identity.txt"
        test_file.write_text("""5.0 IIII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        pauli_dict = H.get_all_pauli_strings(return_as="strings")
        assert len(pauli_dict) == 1
        assert pauli_dict["IIII"] == 5.0

    def test_no_identity_term(self, config, tmp_path):
        """Test Hamiltonian without identity term."""
        test_file = tmp_path / "test_no_identity.txt"
        test_file.write_text("""0.5 XIII
0.3 IXII
""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        # Identity should be 0.0 (or not present)
        identity = pauli_dict.get(tuple(), 0.0)
        assert identity == 0.0

    def test_json_complex_hamiltonian(self, config, tmp_path):
        """Test loading a more complex JSON Hamiltonian."""
        import json
        test_file = tmp_path / "test_complex.json"
        data = {
            "n_qubits": 3,
            "terms": [
                {"ops": [[0, "X"], [1, "Y"], [2, "Z"]], "coeff": 0.5},
                {"ops": [[0, "Z"]], "coeff": 1.5},
                {"ops": [[1, "Z"]], "coeff": 1.5},
                {"ops": [[2, "Z"]], "coeff": 1.5}
            ]
        }
        test_file.write_text(json.dumps(data))

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))
        H = load_pauli(config_ham)

        assert H.num_qubits() == 3
        pauli_dict = H.get_all_pauli_strings(return_as="strings")
        assert pauli_dict["XYZ"] == 0.5
        assert pauli_dict["ZII"] == 1.5

    def test_invalid_extension(self, config, tmp_path):
        """Test that invalid file extension raises ValueError."""
        test_file = tmp_path / "test.xyz"
        test_file.write_text("""1.0 IIII""")

        config_gen, config_ham = config
        config_ham.load_pauli_strings(str(test_file))

        with pytest.raises(ValueError, match="Invalid file extension"):
            load_pauli(config_ham)

    def test_load_json_test_file(self, config):
        """Test loading the data_pauli.json file."""
        test_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "data_pauli.json"
        )

        config_gen, config_ham = config
        config_ham.load_pauli_strings(test_file)
        H = load_pauli(config_ham)

        # Verify basic properties
        assert H.num_qubits() == 4
        pauli_dict = H.get_all_pauli_strings(return_as="tuples")
        assert len(pauli_dict) == 10

        # Verify specific terms
        assert pauli_dict[tuple()] == 2.5
        assert pauli_dict[((0, 'X'),)] == 0.5
        assert pauli_dict[((0, 'X'), (1, 'X'))] == 1.2
        assert pauli_dict[((0, 'X'), (1, 'Y'), (2, 'Z'))] == 0.25
