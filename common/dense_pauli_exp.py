import attrs
from cirq import LineQubit, PauliStringPhasorGate
from functools import cached_property, reduce
from math import cbrt
import numpy as np
from pyLIQTR.utils.pauli_string_manip import convert_to_dense_pauli_string
from qualtran import Bloq, BloqBuilder, Signature, SoquetT
from qualtran.bloqs.basic_gates import Identity, XGate, YGate, ZGate
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import TComplexity, t_complexity
from qualtran.symbolics import SymbolicFloat
from typing import Dict, Sequence

# =================================================================================================

# The following are from https://github.com/quantumlib/Qualtran/issues/360

def add_from_bloq_register_flat_qubits(
    bb: 'BloqBuilder', cirq_bloq: Bloq, **regs: SoquetT
#) -> Tuple[SoquetT, ...]:
):
    """Helper function to split / join registers for cirq gates expeciting single 'qubits' register.
    Args:
        bb: Bloq builder used during decompostion.
        cirq_bloq: A CirqGateAsBloq wrapped arithmetic gate.
        regs: bloq registers we wish to use as flat list of qubits for cirq gate.
    Returns:
        regs: bloq registers appropriately rejoined following split.
    """
    flat_regs = []
    for _, v in regs.items():
        if v.reg.bitsize == 1:
            flat_regs.append([v])
        else:
            flat_regs.append(bb.split(v))
    qubits = np.concatenate(flat_regs)
    qubits = bb.add(cirq_bloq, q=qubits)
    out_soqs = {}
    start = 0
    for _, v in regs.items():
        if v.reg.bitsize == 1:
            end = start + 1
            out_soqs[v] = qubits[start:end][0]
            start += 1
        else:
            end = start + v.reg.bitsize
            out_soqs[v] = bb.join(qubits[start:end])
            start += v.reg.bitsize
    return tuple(s for _, s in out_soqs.items()) 

def add_from_bloq_registers(
    bb: 'BloqBuilder', cirq_bloq: Bloq, **bloq_regs: SoquetT
#) -> Tuple[SoquetT, ...]:
):
    """Shift from bitsize=n, shape=() to bitsize=1, shape=(n,) and back again"""
    print(f"bloq_regs = {bloq_regs}")
    cirq_regs = {}
    for reg_name, soq in bloq_regs.items():
        cirq_regs[reg_name] = bb.split(soq)
    print(f"cirq_regs = {cirq_regs}")
    cirq_regs = bb.add(cirq_bloq, **cirq_regs)
    out_soqs = {}
    for ix, (reg_name, soq) in enumerate(bloq_regs.items()):
        print(f"cirq_regs[{ix}] = {repr(cirq_regs[ix])} ({type(cirq_regs[ix])})")
        print(np.array(cirq_regs[ix]).shape)
        out_soqs[reg_name] = bb.join(cirq_regs[ix])
    return tuple(s for _, s in out_soqs.items())

# =================================================================================================

def _get_pauli_gate(c):
    if c == 'I':
        return Identity()
    elif c == 'X':
        return XGate()
    elif c == 'Y':
        return YGate()
    elif c == 'Z':
        return ZGate()
    else:
        raise ValueError(f"Invalid character: '{c}' does not represent a Pauli gate.")

# =================================================================================================

class DensePauliStringBase(Bloq):
    pass

# =================================================================================================

@attrs.frozen
class DensePauliExponential(Bloq):
    base: DensePauliStringBase

    @cached_property
    def signature(self) -> Signature:
        return Signature.build(q=self.base.num_qubits)
    def _build_bloq(self):
        # This is an awkward way of working around quirks of Cirq and Qualtran
        base_gate = self.base._cirq_gate()
        theta = base_gate.coefficient
        # If theta has a real component, then exponentiating won't give us a unitary operator
        assert(abs(theta.real) < 1.0e-12)
        # If the eigenvalues of M are all +1 or -1, then exp(i theta M) = M^(-2 theta / pi)
        power = - theta.imag * 2.0 / np.pi
        exp_gate = PauliStringPhasorGate(base_gate.copy(coefficient=1))**power
        return CirqGateAsBloq(exp_gate)
    def build_composite_bloq(self, bb: 'BloqBuilder', **soqs: SoquetT) -> Dict[str, 'SoquetT']:
        # Get the register (name agnostic)
        assert(len(soqs) == 1)
        soq = list(soqs.items())[0]
        reg_name = soq[0]
        register = soq[1]
        # Split the register into individual qubits
        register_s = bb.split(register)
        # Add the internal bloq
        # TODO: This doesn't account for the CirqGateAsBloq operating on less qubits than the full
        #       dense Pauli string is aware of.
        internal_bloq = self._build_bloq()
        register_s = bb.add(internal_bloq, q=register_s)
        # Merge the individual qubits to form the output register
        register = bb.join(register_s)
        return { reg_name : register }
    def _t_complexity_(self) -> TComplexity:
        return t_complexity(self._build_bloq())
    def __pow__(self, power):
        return DensePauliExponential(
                attrs.evolve(
                    self.base,
                    coefficient=self.base.coefficient * power
                )
               )
    def __str__(self):
        return f"exp({self.base.coefficient} * {self.base.string})"
    def __repr__(self):
        return f"exp({self.base.coefficient} * {self.base.string})"

# =================================================================================================

@attrs.frozen
class DensePauliString(DensePauliStringBase):
    """A Qualtran Bloq that represents a dense Pauli string.

    Identity entries are simply no-operations and don't appear in the decomposition.  Non-identity
    entries are represented by Qualtran's Pauli gates (XGate, YGate, ZGate)."""
    condensed: tuple
    coefficient: complex
    num_qubits: int

    @classmethod
    def from_dense_pauli_string(cls, string: str, coefficient: complex):
        for c in string:
            assert(c in 'IXYZ')
        return cls(
                condensed=tuple((string[n],n) for n in range(len(string)) if string[n] != 'I'),
                coefficient=coefficient,
                num_qubits=len(string)
               )
    @cached_property
    def string(self) -> str:
        chars = ['I',] * self.num_qubits
        for c, n in self.condensed:
            chars[n] = c
        return "".join(chars)
    @cached_property
    def num_paulis(self):
        return len(self.condensed)
    @cached_property
    def signature(self) -> Signature:
        return Signature.build(q=self.num_qubits)
    def _cirq_gate(self):
        return convert_to_dense_pauli_string((self.string,self.coefficient))
    def build_composite_bloq(self, bb: 'BloqBuilder', **soqs: SoquetT) -> Dict[str, 'SoquetT']:
        # Get the register (name agnostic)
        assert(len(soqs) == 1)
        soq = list(soqs.items())[0]
        reg_name = soq[0]
        register = soq[1]
        # Split the register into individual qubits
        register_s = bb.split(register)
        # Add all the gates
        for c, n in self.condensed:
            if c == 'I':
                pass
            else:
                register_s[n] = bb.add(_get_pauli_gate(c), q=register_s[n])
        # Merge the individual qubits to form the output register
        register = bb.join(register_s)
        return { reg_name : register }
    def _t_complexity_(self) -> TComplexity:
        # TODO: I'm not convinced this is actually correct.  I'm also not convinced we'll ever want
        #       to compute the TComplexity of a DensePauliString.
        tc = TComplexity()
        for c, _ in self.condensed:
            tc += t_complexity(_get_pauli_gate(c))
        return tc
    def __mul__(self, coef: complex):
        return DensePauliString(self.condensed, coef * self.coefficient, self.num_qubits)
    def __rmul__(self, coef: complex):
        return DensePauliString(self.condensed, coef * self.coefficient, self.num_qubits)
    def exp(self):
        return DensePauliExponential(self)
    def __str__(self):
        return f"{self.coefficient} * {self.string}"
    def __repr__(self):
        return f"{self.coefficient} * {self.string}"

