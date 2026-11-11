"""The `algorithm.controlled_power_method` setting (plan T11): sign-controlled QPE from the pipeline."""

import logging

import pytest

from qhat.analysis.algorithm import build_algorithm
from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.analysis.qubit_summary import log_qubit_summary
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.pauli_string_evolution import PauliStringEvolution
from qhat.common.trotter_flattened import Trotterization

# Three qubits, non-commuting terms, an odd time and an even step count
TERMS = (("XZI", 0.5), ("IZY", 0.3), ("ZXX", 0.2))


def trotter(method="second order", num_steps=4, **kwargs):
    return Trotterization.from_method(TERMS, method, 0.7, num_steps, **kwargs)


def make_config(method="QPE: QHAT flexible", **kwargs):
    config = AlgorithmConfiguration()
    config.method = method
    config.num_phase_qubits = 3
    for key, value in kwargs.items():
        setattr(config, key, value)
    return config


class TestSetting:
    def test_unset_is_controlled_and_not_in_toml(self):
        config = make_config()
        assert build_algorithm(config, trotter()).controlled_power_method == "controlled"
        assert "controlled_power_method" not in config._generate_TOML_table()

    def test_set_is_in_toml(self):
        table = make_config(controlled_power_method="sign")._generate_TOML_table()
        assert table["controlled_power_method"] == "sign"

    @pytest.mark.parametrize("value", ["sign", "Sign"])
    def test_sign_reaches_the_bloq(self, value):
        u = trotter()
        qpe = build_algorithm(make_config(controlled_power_method=value), u)
        assert qpe == FlexibleQPE(u, 3, controlled_power_method="sign")

    def test_explicit_controlled(self):
        u = trotter()
        assert (build_algorithm(make_config(controlled_power_method="controlled"), u)
                == FlexibleQPE(u, 3))

    def test_sized_from_requirements(self):
        """Sizing from energy_error still works, and keeps the method."""
        config = make_config(controlled_power_method="sign", num_phase_qubits=None,
                             energy_error=0.05, probability_of_failure=0.1)
        u = trotter()
        qpe = build_algorithm(config, u)
        assert qpe.controlled_power_method == "sign"
        default = build_algorithm(make_config(num_phase_qubits=None, energy_error=0.05,
                                              probability_of_failure=0.1), u)
        assert qpe.num_ancilla_qubits == default.num_ancilla_qubits

    def test_unknown_value_raises(self):
        with pytest.raises(ValueError, match="controlled_power_method must be one of"):
            build_algorithm(make_config(controlled_power_method="signed"), trotter())

    @pytest.mark.parametrize("method", ["QPE: qualtran textbook", "time evolution",
                                        "controlled time evolution"])
    def test_other_methods_raise(self, method):
        with pytest.raises(ValueError, match="only used by \"QPE: QHAT flexible\""):
            build_algorithm(make_config(method, controlled_power_method="sign"), trotter())

    def test_other_methods_raise_even_for_controlled(self):
        """Set means set: "controlled" with another method also raises."""
        with pytest.raises(ValueError, match="only used by"):
            build_algorithm(make_config("time evolution", controlled_power_method="controlled"),
                            trotter())

    def test_non_trotter_unitary_raises(self):
        with pytest.raises(ValueError, match="flattened Trotterization"):
            build_algorithm(make_config(controlled_power_method="sign"),
                            PauliStringEvolution("XY", coefficient=1.0, time=0.5))

    @pytest.mark.parametrize("method", ["first order", "third order"])
    def test_non_symmetric_trotter_raises(self, method):
        with pytest.raises(ValueError, match="symmetric"):
            build_algorithm(make_config(controlled_power_method="sign"), trotter(method))

    def test_odd_step_count_raises(self):
        with pytest.raises(ValueError, match="even number of Trotter steps"):
            build_algorithm(make_config(controlled_power_method="sign"), trotter(num_steps=3))

    def test_conflicts_with_rotation_ancilla(self):
        with pytest.raises(ValueError, match="conflicts with unitary.controlled_rotation_method"):
            build_algorithm(make_config(controlled_power_method="sign"),
                            trotter(controlled_rotation_method="ancilla"))

    def test_controlled_with_rotation_ancilla_works(self):
        u = trotter(controlled_rotation_method="ancilla")
        qpe = build_algorithm(make_config(controlled_power_method="controlled"), u)
        assert qpe == FlexibleQPE(u, 3)


class TestLog:
    def build_log(self, caplog, config, unitary):
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            build_algorithm(config, unitary)
        return " ".join(caplog.text.split())

    def test_sign(self, caplog):
        text = self.build_log(caplog, make_config(controlled_power_method="sign"), trotter())
        assert "controlled-power method: \"sign\"" in text
        assert "controlled-rotation method: not used" in text

    def test_default(self, caplog):
        text = self.build_log(caplog, make_config(), trotter())
        assert "controlled-power method: \"controlled\"" in text
        assert "controlled-rotation method: \"textbook\"" in text

    def test_qubit_summary(self, caplog):
        u = trotter()
        config = make_config(controlled_power_method="sign")
        algorithm = build_algorithm(config, u)
        caplog.clear()
        with caplog.at_level(logging.INFO, logger="qhat.analysis.qubit_summary"):
            log_qubit_summary(config, u, algorithm)
        text = " ".join(caplog.text.split())
        assert "rotation ancillas: 0 (sign-controlled powers" in text
        assert "total: 6" in text
        assert "WARNING" not in text


class TestDriverEndToEnd:
    """Run the driver on a real input file and check the log, the TOML and the T count."""

    def run_driver(self, monkeypatch, tmp_path, algorithm_lines, steps=4):
        import sys
        from pathlib import Path
        from qhat.analysis import driver
        data = Path(__file__).parent / "data_pauli.json"
        tmp_path.mkdir(parents=True, exist_ok=True)
        config_file = tmp_path / "config.py"
        config_file.write_text(
            f'general.output_directory = "{tmp_path}"\n'
            f'hamiltonian.load_pauli_strings("{data}")\n'
            'unitary.encode_ramped_trotter(energy_error=1e-1, trotter_order="fourth order", '
            f'trotter_steps={steps})\n'
            'algorithm.method = "QPE: QHAT flexible"\nalgorithm.num_phase_qubits = 3\n'
            f'{algorithm_lines}\n'
            'analysis.resource_estimator = ["qualtran", "pyliqtr"]\n')
        monkeypatch.setattr(sys, 'argv', ['driver.py', str(config_file)])
        driver.run()
        return (tmp_path / "analysis.log").read_text(), list(tmp_path.glob("*.toml"))[0].read_text()

    @staticmethod
    def t_counts(toml):
        return {line.split("=")[1].strip() for line in toml.splitlines()
                if line.strip().startswith("T_count =")}

    def test_sign_quarters_the_t_count(self, monkeypatch, tmp_path):
        log, toml = self.run_driver(monkeypatch, tmp_path / "sign",
                                    'algorithm.controlled_power_method = "sign"')
        _, toml_default = self.run_driver(monkeypatch, tmp_path / "default", '')
        assert 'controlled_power_method = "sign"' in toml
        assert "controlled_power_method" not in toml_default
        assert "controlled-power method: \"sign\"" in " ".join(log.split())
        assert "WARNING" not in log
        # Both estimators agree; sign is about a quarter of textbook
        (sign,), (default,) = self.t_counts(toml), self.t_counts(toml_default)
        assert 0.2 < int(sign) / int(default) < 0.3

    def test_odd_steps_raise(self, monkeypatch, tmp_path):
        with pytest.raises(ValueError, match="even number of Trotter steps"):
            self.run_driver(monkeypatch, tmp_path, 'algorithm.controlled_power_method = "sign"',
                            steps=3)
