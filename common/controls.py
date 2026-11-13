"""
Control qubits for QHAT's own bloqs

QHAT's evolution bloqs carry their controls in a field, `ctrl_values`, instead of being wrapped
in Qualtran's generic `Controlled`. `ctrl_values` holds one entry per control qubit: 1 means the
bloq acts when that qubit is |1>, 0 when it is |0>. For example, `ctrl_values=(1, 0)` applies U
when the first control is |1> and the second is |0>, and the identity otherwise.

The controls live in a register `ctrl` of shape `(len(ctrl_values),)`, placed before `q`. Each
`.controlled(...)` prepends one value, so `ctrl[0]` is the most recently added control.

With two or more controls, `multi_control_method` picks the circuit:
- "ancilla" (default): AND the controls into one ancilla once (Qualtran `And` / `MultiAnd`, k-1
  ancillas, about 4(k-1) T), run the singly-controlled bloq on it, then uncompute.
- "direct": no ancilla; every rotation is k-controlled (`add_controlled_rz`, 2^k rotations each).
With 0 or 1 controls it is ignored.
"""

import numbers
from collections import Counter
from typing import TYPE_CHECKING, Dict, Iterable, List, Sequence, Set, Tuple

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, CtrlSpec, QBit, Register, Side, Soquet, SoquetT
from qualtran.bloqs.basic_gates import CNOT, GlobalPhase, Rz
from qualtran.bloqs.mcmt import And, MultiAnd
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity

if TYPE_CHECKING:
    from qualtran import AddControlledT
    from qualtran.resource_counting import BloqCountT

MULTI_CONTROL_METHODS = ("ancilla", "direct")
DEFAULT_MULTI_CONTROL_METHOD = "ancilla"


def validate_ctrl_values(owner: str, ctrl_values: Iterable[int]) -> Tuple[int, ...]:
    """Return `ctrl_values` as a tuple of ints, each 0 or 1; raise ValueError for any other value."""
    values = tuple(ctrl_values)
    for value in values:
        if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value not in (0, 1):
            raise ValueError(f"{owner} control values must be 0 or 1, got {values}.")
    return tuple(int(value) for value in values)


def validate_multi_control_method(owner: str, method: str) -> str:
    """Return `method` if it is one of `MULTI_CONTROL_METHODS`; raise ValueError otherwise."""
    if method not in MULTI_CONTROL_METHODS:
        raise ValueError(
            f"{owner} multi_control_method must be one of {MULTI_CONTROL_METHODS}, got {method!r}.")
    return method


def uses_and_ancilla(ctrl_values: Tuple[int, ...], multi_control_method: str) -> bool:
    """True if the bloq ANDs its controls into an ancilla: 2 or more controls and "ancilla"."""
    return len(ctrl_values) >= 2 and multi_control_method == "ancilla"


def num_ancillas(ctrl_values: Tuple[int, ...], multi_control_method: str) -> int:
    """Ancillas the controls add at the peak: k-1 for "ancilla" with k >= 2 controls, else 0."""
    return len(ctrl_values) - 1 if uses_and_ancilla(ctrl_values, multi_control_method) else 0


def peak_qubits(num_qubits: int, ctrl_values: Tuple[int, ...], multi_control_method: str) -> int:
    """`QubitCount` of a controlled bloq whose own decomposition allocates nothing else.

    For example, 10 qubits with 3 controls: 10 + 3 + 2 = 15 for "ancilla", 13 for "direct".
    """
    return num_qubits + len(ctrl_values) + num_ancillas(ctrl_values, multi_control_method)


def ctrl_value_from_spec(ctrl_spec: 'CtrlSpec | None') -> int:
    """The control value (0 or 1) of a single-qubit `CtrlSpec`; `None` means `CtrlSpec()`.

    Raises ValueError for any other `CtrlSpec`.
    """
    if ctrl_spec is None:
        return 1
    if (ctrl_spec.num_ctrl_reg == 1 and ctrl_spec.qdtypes[0] == QBit()
            and ctrl_spec.shapes[0] == () and int(ctrl_spec.cvs[0]) in (0, 1)):
        return int(ctrl_spec.cvs[0])
    raise ValueError(
        f"Only single-qubit controls on |1> or |0> are supported (CtrlSpec() or "
        f"CtrlSpec(cvs=0)), got {ctrl_spec}.")


def ctrl_registers(ctrl_values: Tuple[int, ...]) -> List[Register]:
    """The `ctrl` register for these controls, as a list (empty if there are none)."""
    if not ctrl_values:
        return []
    return [Register('ctrl', QBit(), shape=(len(ctrl_values),), side=Side.THRU)]


def get_ctrl_system(bloq: Bloq, ctrl_spec: 'CtrlSpec | None') -> Tuple[Bloq, 'AddControlledT']:
    """`Bloq.get_ctrl_system` for a bloq with a `ctrl_values` field and a `ctrl` register.

    Returns the same kind of bloq with the new control value prepended, and an adder that wires
    the new control soquet into `ctrl[0]`.
    """
    value = ctrl_value_from_spec(ctrl_spec)
    controlled = attrs.evolve(bloq, ctrl_values=(value,) + bloq.ctrl_values)
    reg_names = [reg.name for reg in bloq.signature]

    def add_controlled(
        bb: BloqBuilder, ctrl_soqs: Sequence[SoquetT], in_soqs: Dict[str, SoquetT]
    ) -> Tuple[Iterable[SoquetT], Iterable[SoquetT]]:
        (new_ctrl,) = ctrl_soqs
        soqs = dict(in_soqs)
        old_ctrls = list(soqs.pop('ctrl', []))
        out = bb.add_d(controlled, ctrl=np.array([new_ctrl] + old_ctrls), **soqs)
        ctrl_out = out.pop('ctrl')
        if old_ctrls:
            out['ctrl'] = ctrl_out[1:]
        return (ctrl_out[0],), [out[name] for name in reg_names]

    return controlled, add_controlled


def _and_ladder(ctrl_values: Tuple[int, ...]) -> Bloq:
    """The bloq that ANDs k >= 2 controls into one target: `And` for 2, `MultiAnd` for more."""
    if len(ctrl_values) == 2:
        return And(*ctrl_values)
    return MultiAnd(cvs=ctrl_values)


def add_with_and_ancilla(
    bb: BloqBuilder, bloq: Bloq, ctrl: SoquetT, **soqs: SoquetT
) -> Dict[str, SoquetT]:
    """Add `bloq` (with k >= 2 controls) as: AND the controls, the bloq on one control, un-AND.

    The AND target is |1> exactly when the controls match `bloq.ctrl_values`, so the inner bloq is
    `bloq` with `ctrl_values=(1,)`. Returns the output soquets, including `ctrl`.
    """
    ladder = _and_ladder(bloq.ctrl_values)
    if isinstance(ladder, And):
        ctrl, target = bb.add(ladder, ctrl=ctrl)
        junk = None
    else:
        ctrl, junk, target = bb.add(ladder, ctrl=ctrl)
    inner = attrs.evolve(bloq, ctrl_values=(1,))
    out = bb.add_d(inner, ctrl=np.array([target]), **soqs)
    (target,) = out.pop('ctrl')
    if junk is None:
        ctrl = bb.add(ladder.adjoint(), ctrl=ctrl, target=target)
    else:
        ctrl = bb.add(ladder.adjoint(), ctrl=ctrl, junk=junk, target=target)
    return {'ctrl': ctrl, **out}


def and_ancilla_call_graph(bloq: Bloq) -> Set['BloqCountT']:
    """Call graph of `add_with_and_ancilla(bloq)`: the AND, its uncompute, and the bloq on one
    control, each once."""
    ladder = _and_ladder(bloq.ctrl_values)
    return {(ladder, 1), (ladder.adjoint(), 1), (attrs.evolve(bloq, ctrl_values=(1,)), 1)}


def and_ancilla_t_complexity(bloq: Bloq) -> TComplexity:
    """T-complexity of `add_with_and_ancilla(bloq)`: the AND and its uncompute, plus the bloq
    on one control."""
    ladder = _and_ladder(bloq.ctrl_values)
    return (t_complexity(ladder) + t_complexity(ladder.adjoint())
            + t_complexity(attrs.evolve(bloq, ctrl_values=(1,))))


def _sign(value: int) -> int:
    """Z eigenvalue of the control state that activates: -1 for |1>, +1 for |0>."""
    return -1 if value == 1 else 1


def add_controlled_rz(
    bb: BloqBuilder, theta: float, target: Soquet, ctrl: Sequence[Soquet],
    ctrl_values: Tuple[int, ...]
) -> Tuple[Soquet, List[Soquet]]:
    """Add Rz(theta) on `target`, controlled on `ctrl` matching `ctrl_values`; exact, no ancilla.

    With projector Pi = prod_j (1 + s_j Z_j) / 2 onto the control pattern (s_j = -1 for |1>, +1
    for |0>), C^k Rz(theta) = exp(-i theta/2 Z_t Pi) = prod over subsets S of the controls of
    exp(-i theta/2^(k+1) prod_{j in S} s_j Z_t Z_S). The 2^k parity rotations are walked in Gray
    code order, one CNOT (control -> target) between rotations and one at the end, so the target
    parity returns to Z_t. That is 2^k Rz and 2^k CNOTs (none for k = 0).

    One control on |1> gives Rz(theta/2) . CNOT . Rz(-theta/2) . CNOT, the textbook circuit.

    Returns the target and control soquets.
    """
    ctrl = list(ctrl)
    k = len(ctrl_values)
    if k == 0:
        return bb.add(Rz(theta), q=target), ctrl
    gray = _gray_code(k)
    for step, (subset, angle) in enumerate(zip(gray, _controlled_rz_angles(theta, ctrl_values))):
        target = bb.add(Rz(angle), q=target)
        flipped = (subset ^ gray[(step + 1) % 2 ** k]).bit_length() - 1
        ctrl[flipped], target = bb.add(CNOT(), ctrl=ctrl[flipped], target=target)
    return target, ctrl


def _gray_code(k: int) -> List[int]:
    """The 2^k subsets of k controls in Gray-code order, as bit masks, e.g. [0, 1, 3, 2]."""
    return [i ^ (i >> 1) for i in range(2 ** k)]


def _controlled_rz_angles(theta: float, ctrl_values: Tuple[int, ...]) -> List[float]:
    """The Rz angles `add_controlled_rz` applies, in order (theta alone for no controls)."""
    k = len(ctrl_values)
    if k == 0:
        return [theta]
    return [float(np.prod([_sign(ctrl_values[j]) for j in range(k) if subset >> j & 1]))
            * theta / 2 ** k for subset in _gray_code(k)]


def controlled_rz_call_graph(theta: float, ctrl_values: Tuple[int, ...]) -> Counter:
    """The gates of `add_controlled_rz`, counted: its Rz rotations and 2^k CNOTs."""
    counts = Counter(Rz(angle) for angle in _controlled_rz_angles(theta, ctrl_values))
    if ctrl_values:
        counts[CNOT()] += 2 ** len(ctrl_values)
    return counts


def controlled_rz_t_complexity(num_ctrls: int) -> TComplexity:
    """T-complexity of `add_controlled_rz`: 2^k rotations and 2^k CNOTs (1 rotation for k = 0)."""
    if num_ctrls == 0:
        return TComplexity(rotations=1)
    return TComplexity(clifford=2 ** num_ctrls, rotations=2 ** num_ctrls)


def add_controlled_phase(
    bb: BloqBuilder, phi: float, ctrl: Sequence[Soquet], ctrl_values: Tuple[int, ...]
) -> List[Soquet]:
    """Add exp(-i phi) on the control pattern `ctrl_values` (identity elsewhere); exact, no ancilla.

    Peeling off the last control, with P = (1 + s Z) / 2 its projector:
    exp(-i phi Pi' P) = exp(-i phi/2 Pi') . exp(-i (s phi/2) Pi' Z), where the second factor is
    Rz(s phi) on the last control, controlled on the others (`add_controlled_rz`), and the first is
    the same problem with one control fewer. With no controls it is a `GlobalPhase`.

    k controls cost 2^k - 1 rotations and 2^k - 2 CNOTs. Returns the control soquets.
    """
    ctrl = list(ctrl)
    for last in reversed(range(len(ctrl_values))):
        ctrl[last], ctrl[:last] = add_controlled_rz(
            bb, _sign(ctrl_values[last]) * phi, ctrl[last], ctrl[:last], ctrl_values[:last])
        phi = phi / 2
    bb.add(GlobalPhase(np.exp(-1j * phi)))
    return ctrl


def controlled_phase_call_graph(phi: float, ctrl_values: Tuple[int, ...]) -> Counter:
    """The gates of `add_controlled_phase`, counted, following the same recursion."""
    counts = Counter()
    for last in reversed(range(len(ctrl_values))):
        counts += controlled_rz_call_graph(_sign(ctrl_values[last]) * phi, ctrl_values[:last])
        phi = phi / 2
    counts[GlobalPhase(np.exp(-1j * phi))] += 1
    return counts


def controlled_phase_t_complexity(num_ctrls: int) -> TComplexity:
    """T-complexity of `add_controlled_phase`: 2^k - 1 rotations and 2^k - 2 CNOTs."""
    total = TComplexity()
    for k in range(num_ctrls):
        total += controlled_rz_t_complexity(k)
    return total


def controlled_matrix(unitary: np.ndarray, ctrl_values: Tuple[int, ...]) -> np.ndarray:
    """The matrix of U controlled on `ctrl_values`: U on the matching block, identity elsewhere.

    The control qubits are the most significant, in order. For example, with one control on |0>
    the result is [[U, 0], [0, I]].
    """
    if not ctrl_values:
        return unitary
    dim = unitary.shape[0]
    block = int(''.join(str(value) for value in ctrl_values), 2)
    result = np.eye(dim * 2 ** len(ctrl_values), dtype=np.result_type(unitary, complex))
    result[block * dim:(block + 1) * dim, block * dim:(block + 1) * dim] = unitary
    return result


def ctrl_repr(ctrl_values: Tuple[int, ...], multi_control_method: str) -> str:
    """The `repr` suffix for the controls, e.g. ", ctrl_values=(1, 0)"; empty if there are none.

    `multi_control_method` is shown only when it isn't the default.
    """
    text = f", ctrl_values={ctrl_values}" if ctrl_values else ""
    if multi_control_method != DEFAULT_MULTI_CONTROL_METHOD:
        text += f", multi_control_method={multi_control_method!r}"
    return text


def ctrl_str(ctrl_values: Tuple[int, ...], text: str) -> str:
    """`text` marked with its controls, e.g. `C[1](text)`; `text` alone if there are none."""
    if not ctrl_values:
        return text
    return f"C[{','.join(str(value) for value in ctrl_values)}]({text})"
