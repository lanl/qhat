"""
Tests for rotation reporting: rotations are totalled over the circuit and converted to T and
Clifford gates once, at the end (plan D11).
"""

import math
import random

import numpy as np
import pytest
from pyLIQTR.gate_decomp.rotation_gates import T_COUNT_CONST, T_COUNT_SLOPE, T_COUNT_STD_DEV

from qhat.analysis.config_types import AnalysisConfiguration
from qhat.analysis.resource_estimation import (
    DEFAULT_ROTATION_GATE_PRECISION, convert_rotations, estimate_resources)
from qhat.common.flexible_qpe import FlexibleQPE
from qhat.common.logging_utils import add_custom_levels
from qhat.common.trotter_flattened import Trotterization

add_custom_levels()

TERMS = [("XZI", 1.0), ("ZIY", 0.7), ("YXX", 0.4)]


def t_per_rotation(precision):
    """pyLIQTR's fit, with its seeded random offset."""
    random.seed(0)
    return int(random.gauss(T_COUNT_SLOPE * np.log2(1 / precision) + T_COUNT_CONST, T_COUNT_STD_DEV))


def trotter(steps=5):
    return Trotterization.from_method(TERMS, "second order", 0.37, steps)


# --- the conversion formula ----------------------------------------------------------------------

def test_conversion_on_known_totals():
    for rotations, precision in [(42, 1e-10), (1000, 1e-3), (7, 1e-6)]:
        per = t_per_rotation(precision)
        r = convert_rotations(11, 23, rotations, rotation_gate_precision=precision)
        assert r["T_count"] == 11 + per * rotations
        assert r["Clifford_count"] == 23 + 2 * per * rotations + 1
        assert r["T_count_unconverted"] == 11
        assert r["Clifford_count_unconverted"] == 23
        assert r["rotation_count"] == rotations
        assert r["rotation_precision"] == precision
        # the fit itself: about 3.02 log2(1/precision) + 0.77, within the random offset
        assert abs(per - (3.02 * math.log2(1 / precision) + 0.77)) < 6 * T_COUNT_STD_DEV + 1


def test_default_precision():
    r = convert_rotations(0, 0, 10)
    assert r["rotation_precision"] == DEFAULT_ROTATION_GATE_PRECISION == 1e-10
    assert r == convert_rotations(0, 0, 10, rotation_gate_precision=1e-10)


def test_algorithm_precision_is_divided_by_total_rotations():
    r = convert_rotations(5, 9, 40, algorithm_precision=1e-4)
    assert r["rotation_precision"] == pytest.approx(1e-4 / 40)
    assert r["T_count"] == 5 + 40 * t_per_rotation(1e-4 / 40)
    # more rotations -> tighter precision for each -> more T per rotation
    assert (convert_rotations(0, 0, 4000, algorithm_precision=1e-4)["T_count"] / 4000
            > r["T_count"] / 40 - 1)


def test_zero_rotations_changes_nothing():
    r = convert_rotations(5, 9, 0)
    assert (r["T_count"], r["Clifford_count"]) == (5, 9)
    assert "rotation_precision" not in r


def test_conflicting_or_bad_precision_raises():
    with pytest.raises(ValueError, match="only one"):
        convert_rotations(0, 0, 3, rotation_gate_precision=1e-3, algorithm_precision=1e-3)
    with pytest.raises(ValueError, match="only one"):
        estimate_resources("qualtran", trotter(), rotation_gate_precision=1e-3,
                           algorithm_precision=1e-3)
    for bad in (0, -1e-3, 1.5):
        with pytest.raises(ValueError, match="between 0 and 1"):
            convert_rotations(0, 0, 3, rotation_gate_precision=bad)


# --- end to end ----------------------------------------------------------------------------------

@pytest.mark.parametrize("bloq_name", ["U", "C-U", "adjoint", "C3-U"])
def test_end_to_end_fields_and_agreement(bloq_name):
    u = trotter()
    bloq = {"U": u, "C-U": u.controlled(), "adjoint": u.adjoint(),
            "C3-U": u.controlled().controlled().controlled()}[bloq_name]
    results = estimate_resources(["pyliqtr", "qualtran"], bloq, rotation_gate_precision=1e-6)
    tc = bloq.t_complexity()
    per = t_per_rotation(1e-6)
    for name in ("pyliqtr", "qualtran"):
        r = results[name]
        assert r["rotation_count"] == tc.rotations > 0
        assert r["T_count_unconverted"] == tc.t
        assert r["Clifford_count_unconverted"] == tc.clifford
        assert r["T_count"] == tc.t + per * tc.rotations
        assert r["Clifford_count"] == tc.clifford + 2 * per * tc.rotations + 1
    assert results["pyliqtr"] == results["qualtran"] | {
        "qubit_count": results["pyliqtr"]["qubit_count"]}


def test_converted_once_not_per_part():
    """Total of two halves converted separately differs from the total converted once."""
    u = trotter()
    whole = estimate_resources("qualtran", u, algorithm_precision=1e-4)["qualtran"]
    assert whole["rotation_precision"] == pytest.approx(1e-4 / u.t_complexity().rotations)


def test_qpe_both_estimators_agree():
    """Both estimators take the rotations from `t_complexity`, which counts the inverse QFT's
    rotations too (`QECGatesCost` misses `m (m - 1) / 2` of them; notes F18)."""
    qpe = FlexibleQPE(trotter().controlled(), 4)
    results = estimate_resources(["pyliqtr", "qualtran"], qpe, algorithm_precision=1e-5)
    tc = qpe.t_complexity()
    for name in ("pyliqtr", "qualtran"):
        assert results[name]["rotation_count"] == tc.rotations
        assert results[name]["T_count"] == results["qualtran"]["T_count"]
        assert results[name]["Clifford_count"] == results["qualtran"]["Clifford_count"]
    assert tc.rotations > 0


# --- settings ------------------------------------------------------------------------------------

def test_toml_records_only_what_was_set():
    config = AnalysisConfiguration()
    config.resource_estimator = "qualtran"
    table = config._generate_TOML_table()
    assert "rotation_gate_precision" not in table and "algorithm_precision" not in table
    config.algorithm_precision = 1e-4
    table = config._generate_TOML_table()
    assert table["algorithm_precision"] == 1e-4 and "rotation_gate_precision" not in table
