"""Qubit summary (plan D3) and the `algorithm.control_values` setting that feeds it (plan T6)."""

import logging

import pytest

from qhat.analysis.algorithm import (
    build_algorithm, check_algorithm_settings, control_values)
from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.analysis.qubit_summary import log_qubit_summary
from qhat.analysis.resource_estimation import estimate_resources
from qhat.common.trotter_flattened import Trotterization
from qualtran.resource_counting import get_cost_value, QubitCount

# Three qubits, non-commuting terms, an odd time and several steps
TERMS = (("XZI", 0.5), ("IZY", 0.3), ("ZXX", 0.2))


def trotter(**kwargs):
    return Trotterization(TERMS, (0.5, 0.5), 0.7, 3, **kwargs)


def make_config(method, **kwargs):
    config = AlgorithmConfiguration()
    config.method = method
    for key, value in kwargs.items():
        setattr(config, key, value)
    return config


def summary(caplog, config, unitary=None):
    unitary = unitary or trotter()
    algorithm = build_algorithm(config, unitary)
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="qhat.analysis.qubit_summary"):
        log_qubit_summary(config, unitary, algorithm)
    return " ".join(caplog.text.split()), algorithm    # one space between words


class TestSummaryContent:
    def test_plain_time_evolution(self, caplog):
        text, _ = summary(caplog, make_config("time evolution"))
        assert "system qubits: 3" in text
        assert "phase qubits" not in text
        assert "control qubits" not in text
        assert "multi-control ancillas: 0 (no controls)" in text
        assert "total: 3 (Qualtran QubitCount)" in text
        assert "WARNING" not in text

    def test_flexible_qpe(self, caplog):
        config = make_config("QPE: QHAT flexible", num_phase_qubits=4)
        text, _ = summary(caplog, config)
        assert "system qubits: 3" in text
        assert "phase qubits: 4" in text
        assert "single control per unitary; method not used" in text
        assert "total: 7" in text
        assert "WARNING" not in text

    def test_one_control(self, caplog):
        text, _ = summary(caplog, make_config("controlled time evolution"))
        assert "control qubits: 1" in text
        assert "multi-control ancillas: 0 (single control; method not used)" in text
        assert "total: 4" in text
        assert "WARNING" not in text

    @pytest.mark.parametrize("method, ancillas", [("ancilla", 2), ("direct", 0)])
    def test_three_controls(self, caplog, method, ancillas):
        config = make_config("controlled time evolution", control_values=[1, 0, 1])
        text, algorithm = summary(caplog, config, trotter(multi_control_method=method))
        assert "control qubits: 3" in text
        assert f"multi-control ancillas: {ancillas} (3 controls; method \"{method}\")" in text
        assert f"total: {3 + 3 + ancillas} (Qualtran" in text
        assert get_cost_value(algorithm, QubitCount()) == 3 + 3 + ancillas
        assert "WARNING" not in text

    def test_warns_when_parts_disagree(self, caplog, monkeypatch):
        """If the Qualtran total isn't the sum of the parts, say so."""
        monkeypatch.setattr("qhat.analysis.qubit_summary.get_cost_value",
                            lambda algorithm, key: 5)
        config = make_config("controlled time evolution")
        unitary = trotter()
        algorithm = build_algorithm(config, unitary)
        with caplog.at_level(logging.INFO, logger="qhat.analysis.qubit_summary"):
            log_qubit_summary(config, unitary, algorithm)
        assert "components add up to 4, but Qualtran QubitCount says 5" in caplog.text

    def test_textbook_qpe_total_unavailable(self, caplog):
        """Known issue K1: Qualtran's QubitCount crashes on TextbookQPE; the summary survives."""
        config = make_config("QPE: qualtran textbook", num_phase_qubits=3)
        text, _ = summary(caplog, config)
        assert "system qubits: 3" in text
        assert "phase qubits: 3" in text
        assert "total is unavailable" in text


class TestControlValuesSetting:
    def test_unset_is_one_control_and_not_in_toml(self):
        config = make_config("controlled time evolution")
        assert control_values(config) == (1,)
        assert "control_values" not in config._generate_TOML_table()

    def test_set_is_in_toml(self):
        config = make_config("controlled time evolution", control_values=[1, 0])
        assert list(config._generate_TOML_table()["control_values"]) == [1, 0]

    def test_values_reach_the_bloq_in_order(self):
        config = make_config("controlled time evolution", control_values=[1, 0, 0])
        algorithm = build_algorithm(config, trotter())
        assert algorithm.ctrl_values == (1, 0, 0)
        assert algorithm.num_qubits == 3

    def test_log_reports_the_method(self, caplog):
        config = make_config("controlled time evolution", control_values=[1, 1, 0])
        with caplog.at_level(logging.INFO, logger="qhat.analysis.algorithm"):
            build_algorithm(config, trotter())
        assert '"ancilla" (3 controls), 2 ancillas added' in caplog.text

    @pytest.mark.parametrize("bad", [[2], [1, "a"], [True]])
    def test_bad_entries_raise(self, bad):
        config = make_config("controlled time evolution", control_values=bad)
        with pytest.raises(ValueError, match="control values must be 0 or 1"):
            check_algorithm_settings(config, trotter())

    def test_empty_raises(self):
        config = make_config("controlled time evolution", control_values=[])
        with pytest.raises(ValueError, match="empty"):
            check_algorithm_settings(config, trotter())

    @pytest.mark.parametrize("method", [
        "time evolution", "QPE: QHAT flexible", "QPE: qualtran textbook"])
    def test_other_methods_raise(self, method):
        """QPE wraps the unitary in Qualtran's `Controlled`, which would repeat the AND ladder
        for every term (R2 finding 2), so it must never see several controls."""
        config = make_config(method, num_phase_qubits=3, control_values=[1, 1])
        with pytest.raises(ValueError, match="only used by \"controlled time evolution\""):
            check_algorithm_settings(config, trotter())

    @pytest.mark.parametrize("method, values, match", [
        ("QPE: QHAT flexible", [1, 1], "only used by"),
        ("QPE: qualtran textbook", [1], "only used by"),
        ("time evolution", [0], "only used by"),
        ("controlled time evolution", [], "empty"),
        ("controlled time evolution", [1, 2], "must be 0 or 1")])
    def test_build_algorithm_checks_the_setting(self, method, values, match):
        """`build_algorithm` runs the check itself, so no caller can skip it. Without it, QPE
        would ignore `control_values` while the TOML recorded them (R3 finding 2)."""
        config = make_config(method, num_phase_qubits=3, control_values=values)
        with pytest.raises(ValueError, match=match):
            build_algorithm(config, trotter())

    def test_unitary_without_multi_control_raises(self):
        from qhat.common.pauli_string_evolution import PauliStringEvolution
        config = make_config("controlled time evolution", control_values=[1, 1])
        with pytest.raises(ValueError, match="doesn't support"):
            check_algorithm_settings(config, object())
        # A single control is fine for anything
        check_algorithm_settings(make_config("controlled time evolution", control_values=[1]),
                                 object())


class TestQubitCountInResources:
    @pytest.mark.parametrize("method, ancillas", [("ancilla", 2), ("direct", 0)])
    def test_pyliqtr_reports_the_peak(self, method, ancillas):
        """qubit_count is the peak: 3 system + 3 controls + the AND ancillas."""
        config = make_config("controlled time evolution", control_values=[1, 0, 1])
        algorithm = build_algorithm(config, trotter(multi_control_method=method))
        resources = estimate_resources("pyliqtr", algorithm)["pyliqtr"]
        assert resources["qubit_count"] == get_cost_value(algorithm, QubitCount())
        assert resources["qubit_count"] == 6 + ancillas

    def test_pyliqtr_equals_qualtran_for_qpe(self):
        config = make_config("QPE: QHAT flexible", num_phase_qubits=3)
        algorithm = build_algorithm(config, trotter())
        pyliqtr = estimate_resources("pyliqtr", algorithm)["pyliqtr"]
        qualtran = estimate_resources("qualtran", algorithm)["qualtran"]
        assert pyliqtr["qubit_count"] == qualtran["qubit_count"] == 6

    def test_textbook_qpe_falls_back_to_logical_qubits(self, caplog):
        config = make_config("QPE: qualtran textbook", num_phase_qubits=3)
        algorithm = build_algorithm(config, trotter())
        with caplog.at_level(logging.WARNING):
            resources = estimate_resources("pyliqtr", algorithm)["pyliqtr"]
        assert resources["qubit_count"] == 6
        assert "counts registers only" in caplog.text


class TestDriverEndToEnd:
    """Run the driver on a real input file and check the log and the TOML."""

    def run_driver(self, monkeypatch, tmp_path, algorithm_lines):
        import sys
        from pathlib import Path
        from qhat.analysis import driver
        data = Path(__file__).parent / "data_pauli.json"
        config_file = tmp_path / "config.py"
        config_file.write_text(
            f'general.output_directory = "{tmp_path}"\n'
            f'hamiltonian.load_pauli_strings("{data}")\n'
            'unitary.encode_ramped_trotter(energy_error=1e-1, trotter_order="second order", trotter_steps=3)\n'
            f'{algorithm_lines}\n'
            'analysis.resource_estimator = ["qualtran", "pyliqtr"]\n')
        monkeypatch.setattr(sys, 'argv', ['driver.py', str(config_file)])
        driver.run()
        return (tmp_path / "analysis.log").read_text(), list(tmp_path.glob("*.toml"))[0].read_text()

    def test_three_controls(self, monkeypatch, tmp_path):
        log, toml = self.run_driver(
            monkeypatch, tmp_path,
            'algorithm.method = "controlled time evolution"\n'
            'algorithm.control_values = [1, 0, 1]')
        text = " ".join(log.split())
        assert "control qubits: 3" in text
        assert "multi-control ancillas: 2 (3 controls; method \"ancilla\")" in text
        assert "control_values = [1, 0, 1]" in toml
        # Qualtran and pyLIQTR both report the peak
        assert toml.count("qubit_count") == 2
        counts = {line.split("=")[1].strip() for line in toml.splitlines() if "qubit_count" in line}
        assert len(counts) == 1

    def test_qpe_with_control_values_raises(self, monkeypatch, tmp_path):
        with pytest.raises(ValueError, match="only used by \"controlled time evolution\""):
            self.run_driver(
                monkeypatch, tmp_path,
                'algorithm.method = "QPE: QHAT flexible"\nalgorithm.num_phase_qubits = 3\n'
                'algorithm.control_values = [1, 1]')

    def test_flexible_qpe(self, monkeypatch, tmp_path):
        log, toml = self.run_driver(
            monkeypatch, tmp_path,
            'algorithm.method = "QPE: QHAT flexible"\nalgorithm.num_phase_qubits = 3')
        text = " ".join(log.split())
        assert "phase qubits: 3" in text
        assert "WARNING" not in log
        assert "control_values" not in toml
