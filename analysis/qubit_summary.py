"""Qubit summary: where the qubits of the algorithm go, logged once near the results (plan D3).

The total is Qualtran's `QubitCount`, the peak number of qubits in use at any point, including
temporary ancillas. The components are worked out from the settings; if they don't add up to the
total, a warning says so.

Example (QPE with 9 phase qubits on a 12-qubit system):
    Qubit summary:
      system qubits:          12
      phase qubits:            9
      multi-control ancillas:  0  (single control; method not used)
      rotation ancillas:       0  (controlled-rotation method "textbook")
      total:                  21  (Qualtran QubitCount)
"""

import logging

from qualtran.resource_counting import get_cost_value, QubitCount

from qhat.analysis.config_types import AlgorithmConfiguration
from qhat.common import controls

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------------------------------

def num_phase_qubits(algorithm) -> int:
    """Size of the QPE phase register, or 0 if the algorithm has none."""
    for name in ("num_ancilla_qubits", "m_bits"):    # FlexibleQPE, Qualtran's TextbookQPE
        if hasattr(algorithm, name):
            return int(getattr(algorithm, name))
    return 0

def num_system_qubits(unitary) -> int:
    """Qubits the unitary acts on, not counting controls."""
    if hasattr(unitary, "num_qubits"):
        return int(unitary.num_qubits)
    return int(unitary.signature.n_qubits())

# -------------------------------------------------------------------------------------------------

def qubit_summary_components(config_algorithm: AlgorithmConfiguration, unitary, algorithm) -> dict:
    """Qubit counts by role: system, phase, control, ancilla and the method line.

    Returns `{}` for algorithms whose layout isn't broken down (the qubitized QPE variants).
    """
    method = config_algorithm.method.lower()
    if method not in ("time evolution", "controlled time evolution",
                      "qpe: qualtran textbook", "qpe: qhat flexible"):
        return {}
    components = {"system": num_system_qubits(unitary),
                  "phase": num_phase_qubits(algorithm),
                  "control": 0,
                  "ancilla": 0,
                  "ancilla_note": "no controls",
                  "rotation_ancilla": 0,
                  "rotation_note": None}
    if method == "controlled time evolution":
        ctrl_values = tuple(algorithm.ctrl_values) if hasattr(algorithm, "ctrl_values") else (1,)
        multi_method = getattr(algorithm, "multi_control_method", None)
        components["control"] = len(ctrl_values)
        if len(ctrl_values) <= 1:
            components["ancilla_note"] = "single control; method not used"
        elif multi_method is not None:
            components["ancilla"] = controls.num_ancillas(ctrl_values, multi_method)
            components["ancilla_note"] = (
                f"{len(ctrl_values)} controls; method \"{multi_method}\"")
    elif method.startswith("qpe:"):
        # The phase qubits are the controls, one per controlled U^(2^j)
        components["ancilla_note"] = "single control per unitary; method not used"
    rotation_method = getattr(unitary, "controlled_rotation_method", None)
    if (method == "qpe: qhat flexible"
            and getattr(algorithm, "controlled_power_method", None) == "sign"):
        # Sign-controlled powers control no rotation; the phase qubit is part of each term
        components["rotation_note"] = "sign-controlled powers; controlled-rotation method not used"
    elif method != "time evolution" and rotation_method is not None:
        # One ancilla, reused by every controlled rotation (none if the unitary has no rotation)
        components["rotation_ancilla"] = controls.num_rotation_ancillas(
            (1,), rotation_method, has_rotations(unitary))
        components["rotation_note"] = f"controlled-rotation method \"{rotation_method}\""
    return components

def has_rotations(unitary) -> bool:
    """False only if every term of the unitary is all identities (then nothing is rotated)."""
    terms = getattr(unitary, "pauli_terms", None)
    return terms is None or any(set(ps) != {"I"} for ps, _ in terms)

def log_qubit_summary(config_algorithm: AlgorithmConfiguration, unitary, algorithm):
    """Log the qubit summary block (INFO), with a warning if the parts don't add up."""
    try:
        total = int(get_cost_value(algorithm, QubitCount()))
        source = "Qualtran QubitCount"
    except Exception as error:    # e.g. known issue K1: "QPE: qualtran textbook" has no counter
        logger.warning(f"Qubit summary: Qualtran QubitCount failed ({type(error).__name__}); "
                       "the total is unavailable.")
        total = None
        source = None

    parts = qubit_summary_components(config_algorithm, unitary, algorithm)
    lines = ["Qubit summary:"]
    if parts:
        lines.append(f"  system qubits:          {parts['system']:>4}")
        if parts["phase"]:
            lines.append(f"  phase qubits:           {parts['phase']:>4}")
        if parts["control"]:
            lines.append(f"  control qubits:         {parts['control']:>4}")
        lines.append(f"  multi-control ancillas: {parts['ancilla']:>4}  ({parts['ancilla_note']})")
        if parts["rotation_note"]:
            lines.append(f"  rotation ancillas:      {parts['rotation_ancilla']:>4}  "
                         f"({parts['rotation_note']})")
    if total is not None:
        lines.append(f"  total:                  {total:>4}  ({source})")
    logger.info("\n".join(lines))

    if parts and total is not None:
        added = (parts["system"] + parts["phase"] + parts["control"] + parts["ancilla"]
                 + parts["rotation_ancilla"])
        if added != total:
            logger.warning(
                f"Qubit summary: the components add up to {added}, but {source} says {total}.")
