"""The `multi_control_method` unitary option (plan D1): TOML, validation, plumbing and logging."""

import logging

import pytest

from qhat.analysis.algorithm import build_algorithm, log_multi_control_method
from qhat.analysis.config_types import AlgorithmConfiguration, UnitaryConfiguration
from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings
from qhat.analysis.unitary import encode_ramped_trotter
from qhat.common.trotter_flattened import Trotterization


@pytest.fixture
def hamiltonian():
    pauli_dict = {((0, 'X'), (1, 'Y')): 0.5, ((1, 'Y'), (2, 'Z')): 0.3, ((0, 'Z'),): 0.2}
    return Hamiltonian(LinearCombinationOfPauliStrings(sparse=pauli_dict, num_qubits=3))


def unitary_config(**kwargs):
    config = UnitaryConfiguration()
    config.encode_ramped_trotter(energy_error=1e-2, ordering_method=None, **kwargs)
    return config


class TestUnitaryConfiguration:
    def test_unset_is_not_in_toml(self):
        config = unitary_config()
        assert config.multi_control_method is None
        assert "multi_control_method" not in config._generate_TOML_table()

    @pytest.mark.parametrize("method", ["ancilla", "direct"])
    def test_set_is_in_toml(self, method):
        table = unitary_config(multi_control_method=method)._generate_TOML_table()
        assert table["multi_control_method"] == method

    def test_set_to_default_is_in_toml(self):
        """The TOML records what the user set, even when it equals the default."""
        table = unitary_config(multi_control_method="ancilla")._generate_TOML_table()
        assert table["multi_control_method"] == "ancilla"

    def test_unknown_value_raises_listing_valid_ones(self):
        with pytest.raises(ValueError, match=r"\['ancilla', 'direct'\]"):
            unitary_config(multi_control_method="toffoli")


class TestEncodeRampedTrotter:
    def test_default_is_ancilla(self, hamiltonian):
        trotter = encode_ramped_trotter(unitary_config(), hamiltonian, tevol_hbar=0.7)
        assert isinstance(trotter, Trotterization)
        assert trotter.multi_control_method == "ancilla"

    def test_direct_reaches_the_bloq(self, hamiltonian):
        trotter = encode_ramped_trotter(
            unitary_config(multi_control_method="direct"), hamiltonian, tevol_hbar=0.7)
        assert trotter.multi_control_method == "direct"
        assert trotter.controlled().controlled().multi_control_method == "direct"

    def test_original_with_method_raises(self, hamiltonian):
        config = unitary_config(trotter_implementation="original", multi_control_method="direct")
        with pytest.raises(ValueError, match="flattened"):
            encode_ramped_trotter(config, hamiltonian, tevol_hbar=0.7)


class TestLogging:
    def trotter(self, **kwargs):
        return Trotterization((("XZ", 0.5), ("ZY", 0.3)), (0.5, 0.5), 0.7, 3, **kwargs)

    @pytest.mark.parametrize("method", ["ancilla", "direct"])
    def test_ignored_with_one_control(self, caplog, method):
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            log_multi_control_method(self.trotter(multi_control_method=method), 1)
        assert f'not used (1 control per unitary; "{method}" is ignored), 0 ancillas' in caplog.text

    @pytest.mark.parametrize("method, ancillas", [("ancilla", 2), ("direct", 0)])
    def test_method_and_ancillas_with_three_controls(self, caplog, method, ancillas):
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            log_multi_control_method(self.trotter(multi_control_method=method), 3)
        assert f'"{method}" (3 controls), {ancillas} ancillas added' in caplog.text

    def test_controlled_time_evolution_logs(self, caplog):
        config = AlgorithmConfiguration()
        config.method = "controlled time evolution"
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            controlled = build_algorithm(config, self.trotter())
        assert controlled.ctrl_values == (1,)
        assert "multi-control method: not used" in caplog.text

    def test_uncontrolled_time_evolution_does_not_log(self, caplog):
        config = AlgorithmConfiguration()
        config.method = "time evolution"
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            build_algorithm(config, self.trotter())
        assert "multi-control method" not in caplog.text


# -------------------------------------------------------------------------------------------------
# The `controlled_rotation_method` unitary option (plan T10), by the same rules
# -------------------------------------------------------------------------------------------------

class TestControlledRotationMethodOption:
    def test_unset_is_not_in_toml(self):
        config = unitary_config()
        assert config.controlled_rotation_method is None
        assert "controlled_rotation_method" not in config._generate_TOML_table()

    @pytest.mark.parametrize("method", ["textbook", "ancilla"])
    def test_set_is_in_toml(self, method):
        """Recorded whenever the user set it, the default included."""
        table = unitary_config(controlled_rotation_method=method)._generate_TOML_table()
        assert table["controlled_rotation_method"] == method

    def test_unknown_value_raises_listing_valid_ones(self):
        with pytest.raises(ValueError, match=r"\['textbook', 'ancilla'\]"):
            unitary_config(controlled_rotation_method="toffoli")

    def test_default_is_textbook(self, hamiltonian):
        trotter = encode_ramped_trotter(unitary_config(), hamiltonian, tevol_hbar=0.7)
        assert trotter.controlled_rotation_method == "textbook"

    def test_ancilla_reaches_the_bloq(self, hamiltonian):
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="ancilla"), hamiltonian, tevol_hbar=0.7)
        assert trotter.controlled_rotation_method == "ancilla"
        assert trotter.controlled().controlled_rotation_method == "ancilla"

    def test_original_with_method_raises(self, hamiltonian):
        config = unitary_config(trotter_implementation="original",
                                controlled_rotation_method="ancilla")
        with pytest.raises(ValueError, match="flattened"):
            encode_ramped_trotter(config, hamiltonian, tevol_hbar=0.7)

    def test_conflict_with_direct_and_several_controls_raises(self, hamiltonian):
        """Two settings that can't both be honoured: the build raises instead of choosing."""
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="ancilla", multi_control_method="direct"),
            hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = "controlled time evolution"
        config.control_values = [1, 0]
        with pytest.raises(ValueError, match="multi_control_method 'direct'"):
            build_algorithm(config, trotter)

    def test_textbook_qpe_with_ancilla_raises(self, hamiltonian):
        """TextbookQPE wraps U in Qualtran's Controlled, which loses the merged phase (R4)."""
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="ancilla"), hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = "QPE: qualtran textbook"
        config.num_phase_qubits = 3
        with pytest.raises(ValueError, match="QPE: QHAT flexible"):
            build_algorithm(config, trotter)

    @pytest.mark.parametrize("method", ["QPE: QHAT flexible", "controlled time evolution"])
    def test_other_controlled_methods_with_ancilla_work(self, hamiltonian, method):
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="ancilla"), hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = method
        config.num_phase_qubits = 3
        build_algorithm(config, trotter)

    def test_textbook_qpe_with_textbook_rotation_works(self, hamiltonian):
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="textbook"), hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = "QPE: qualtran textbook"
        config.num_phase_qubits = 3
        build_algorithm(config, trotter)

    @pytest.mark.parametrize("method, line", [
        ("ancilla", 'controlled-rotation method: "ancilla" (AND ancilla per rotation, 1 rotation '
                    '+ 4 T each), 1 ancilla added'),
        ("textbook", 'controlled-rotation method: "textbook" (2 rotations each), 0 ancillas')])
    def test_log_line_for_controlled_algorithms(self, caplog, hamiltonian, method, line):
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method=method), hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = "controlled time evolution"
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            build_algorithm(config, trotter)
        assert line in caplog.text

    def test_not_logged_for_plain_time_evolution(self, caplog, hamiltonian):
        trotter = encode_ramped_trotter(
            unitary_config(controlled_rotation_method="ancilla"), hamiltonian, tevol_hbar=0.7)
        config = AlgorithmConfiguration()
        config.method = "time evolution"
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            build_algorithm(config, trotter)
        assert "controlled-rotation method" not in caplog.text
