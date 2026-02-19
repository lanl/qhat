# TODO: This file is nowhere near complete testing.  It's at best a quick sanity-check.

from cirq import LineQubit, X, Y, Z, DensePauliString
from dense_pauli_exp import DensePauliString
from numpy import exp
from pyLIQTR.utils.pauli_string_manip import convert_to_dense_pauli_string
from qualtran.bloqs.basic_gates import XGate, YGate
from qualtran.bloqs.phase_estimation import TextbookQPE
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import t_complexity, TComplexity
from trotter import build_ramped_trotterized_unitary

class BasicBloqs:
    def n_qubits(self):
        return 3
    def yield_PauliLCU_Info(self, return_as):
        assert(return_as == "strings")
        yield ("XII", 0.1)
        yield ("IYI", 0.2)
        yield ("IIZ", 0.3)
        yield ("IYZ", 0.4)
        yield ("XIZ", 0.5)
        yield ("XYI", 0.6)
        yield ("XYZ", 0.7)

def tprint(obj):
    print(type(obj), obj)

def main():
    print("-----------------------------------------------------------------------------")
    qubits = LineQubit.range(32)
    qualtran_x = XGate()
    tprint(qualtran_x)
    qualtran_y = YGate()
    tprint(qualtran_y)

    print("-----------------------------------------------------------------------------")
    cirq_x5 = X(qubits[5])
    tprint(cirq_x5)
    qc_x5 = CirqGateAsBloq(cirq_x5)
    tprint(qc_x5)
    tprint(qc_x5.signature)
    qc_pow_x5 = qc_x5**0.34
    tprint(qc_pow_x5)
    tprint(qc_pow_x5.signature)

    print("-----------------------------------------------------------------------------")
    cirq_y9 = Y(qubits[9])
    tprint(cirq_y9)
    qc_y9 = CirqGateAsBloq(cirq_y9)
    tprint(qc_y9)
    tprint(qc_y9.signature)
    qc_pow_y9 = qc_y9**0.34
    tprint(qc_pow_y9)
    tprint(qc_pow_y9.signature)

    print("-----------------------------------------------------------------------------")
    cirq_xy = cirq_x5 * cirq_y9
    tprint(cirq_xy)
    qc_xy = CirqGateAsBloq(cirq_xy)
    tprint(qc_xy)
    tprint(qc_xy.signature)
    qc_pow_xy = qc_xy**0.34
    tprint(qc_pow_xy)
    tprint(qc_pow_xy.signature)

    print("-----------------------------------------------------------------------------")
    print("DensePauliString")
    # TODO: deal with sparsity dps = DensePauliString.from_dense_pauli_string('IIXIIYIIZII', 0.34)
    dps = DensePauliString.from_dense_pauli_string('XYZ', 0.34)
    print(f"- dps: {dps}")
    tprint(dps.string)
    tprint(dps.coefficient)
    tprint(dps.num_qubits)
    tprint(dps.num_paulis)
    tprint(dps.signature)

    print("-----------------------------------------------------------------------------")
    print("exp(i dps)")
    dpe = exp(1j * dps)
    print(f"- dpe: {dpe}")
    print(t_complexity(dpe))

    print("-----------------------------------------------------------------------------")
    print("controlled exp(i dps)")
    c_dpe = dpe.controlled()
    print(t_complexity(c_dpe))

    print("-----------------------------------------------------------------------------")
    # TODO: Why do these two generate different TComplexity counts?  They are equivalent, so they
    #       should be the same, right?  Is it because one of them has a global phase relative to
    #       the other?  I thought the exponential generated a PauliStringPhasorGate and ended up
    #       with the same result (essentially throwing away the global phase difference between
    #       them)?
    print("cirq's version of exp(i dps)")
    q3 = LineQubit.range(3)
    cirq_dpe1 = CirqGateAsBloq(exp(1j * convert_to_dense_pauli_string(('XYZ', 0.34)).on(*q3)))
    print(cirq_dpe1)
    print(t_complexity(cirq_dpe1))
    from numpy import pi
    cirq_dpe2 = CirqGateAsBloq(convert_to_dense_pauli_string(('XYZ', 1)).on(*q3)**(-0.34*2/pi))
    print(cirq_dpe2)
    print(t_complexity(cirq_dpe2))
    #cirq_c_dpe = cirq_dpe.controlled()
    #print(t_complexity(cirq_c_dpe))


if __name__ == "__main__":
    main()
