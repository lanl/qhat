"""
Control qubits for QHAT's own bloqs

QHAT's evolution bloqs carry their controls in a field, `ctrl_values`, instead of being wrapped
in Qualtran's generic `Controlled`. `ctrl_values` holds one entry per control qubit: 1 means the
bloq acts when that qubit is |1>, 0 when it is |0>. For example, `ctrl_values=(1, 0)` applies U
when the first control is |1> and the second is |0>, and the identity otherwise.

The controls live in a register `ctrl` of shape `(len(ctrl_values),)`, placed before `q`. Each
`.controlled(...)` prepends one value, so `ctrl[0]` is the most recently added control.
"""

import numbers
from typing import TYPE_CHECKING, Dict, Iterable, List, Sequence, Tuple

import attrs
import numpy as np

from qualtran import Bloq, BloqBuilder, CtrlSpec, QBit, Register, Side, SoquetT

if TYPE_CHECKING:
    from qualtran import AddControlledT

# Until multiple controls are implemented, a bloq can hold at most this many.
MAX_NUM_CONTROLS = 1


def validate_ctrl_values(owner: str, ctrl_values: Iterable[int]) -> Tuple[int, ...]:
    """Return `ctrl_values` as a tuple of ints, each 0 or 1.

    Raises ValueError for any other value, and NotImplementedError for more than one control.
    """
    values = tuple(ctrl_values)
    for value in values:
        if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value not in (0, 1):
            raise ValueError(f"{owner} control values must be 0 or 1, got {values}.")
    if len(values) > MAX_NUM_CONTROLS:
        raise NotImplementedError(
            f"{owner} supports at most {MAX_NUM_CONTROLS} control qubit, got {len(values)}.")
    return tuple(int(value) for value in values)


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


def ctrl_str(ctrl_values: Tuple[int, ...], text: str) -> str:
    """`text` marked with its controls, e.g. `C[1](text)`; `text` alone if there are none."""
    if not ctrl_values:
        return text
    return f"C[{','.join(str(value) for value in ctrl_values)}]({text})"
