# TODO: This file is nowhere near complete testing.  It's at best a quick sanity-check.

from cirq import LineQubit
from dense_pauli_exp import DensePauliString
from numpy import exp
from pyLIQTR.utils.pauli_string_manip import convert_to_dense_pauli_string
from qualtran.bloqs.phase_estimation import TextbookQPE
from qualtran.cirq_interop import CirqGateAsBloq
from qualtran.cirq_interop.t_complexity_protocol import t_complexity, TComplexity
from trotter import build_ramped_trotterized_unitary

class BasicBloqs:
    def n_qubits(self):
        return 3
    def yield_PauliLCU_Info(self, return_as):
        assert(return_as == "strings")
        yield ("IIIIIIXIIIIIIII", 0.1)
        yield ("IIIIIIIYIIIIIII", 0.2)
        yield ("IIIIIIIIZIIIIII", 0.3)
        yield ("IIIIIIIYZIIIIII", 0.4)
        yield ("IIIIIIXIZIIIIII", 0.5)
        yield ("IIIIIIXYIIIIIII", 0.6)
        yield ("IIIIIIXYZIIIIII", 0.7)
        yield ("IIIIIIYZXIIIIII", 0.8)
        yield ("IIIIIIZXYIIIIII", 0.9)

def build_unitary(problem_instance, coefficients, timestep, numsteps):
    rtu = build_ramped_trotterized_unitary(
            problem_instance,
            method=coefficients,
            timestep=timestep,
            numsteps=numsteps)
    Nr = len(coefficients) * numsteps

    print("Complexity of RampedTrotterizedUnitary --",
          f"{len(coefficients)} coefficients x {numsteps} steps:")
    tc_u = t_complexity(rtu)
    print(f'    T-count:   {tc_u.t:g}')
    print(f'    Rotations: {tc_u.rotations:g}')
    print(f'    Cliffords: {tc_u.clifford:g}')

    tc_c = t_complexity(rtu.controlled())
    print("    Complexity of controlled unitary:")
    print(f'        T-count:   {tc_c.t:g} ({Nr} x {tc_c.t//Nr})')
    print(f'        Rotations: {tc_c.rotations:g} ({Nr} x {tc_c.rotations//Nr})')
    print(f'        Cliffords: {tc_c.clifford:g} ({Nr} * {tc_c.clifford//Nr})')

    # TODO: TextbookQPE runs into problems that I haven't yet resolved
    Nq = 7
    qpe = TextbookQPE(rtu, Nq)
    tc_q = t_complexity(qpe)
    print(f"    Complexity of TextbookQPE<RampedTrotterizedUnitary> -- {Nq} qubits:")
    print(f'        T-count:   {tc_q.t:g}')
    print(f'        Rotations: {tc_q.rotations:g}')
    print(f'        Cliffords: {tc_q.clifford:g}')

    print(f"        Remove (2^{Nq}-1) copies of C-U:")
    k = 2**Nq - 1
    print(f'            T-count:   {tc_q.t - k * tc_c.t:g}')
    print(f'            Rotations: {tc_q.rotations - k * tc_c.rotations:g}')
    print(f'            Cliffords: {tc_q.clifford - k * tc_c.clifford:g}')

def main():
    problem_instance = BasicBloqs()

    ramp_cost = TComplexity()
    qubits = LineQubit.range(problem_instance.n_qubits())
    for term in problem_instance.yield_PauliLCU_Info("strings"):
        bloq = exp(1j * DensePauliString.from_dense_pauli_string(*term))
        bloq_cost = t_complexity(bloq)
        ramp_cost += bloq_cost
        print(f"Complexity of exp(i {term[0]}):",
              f"(t: {bloq_cost.t}, rot: {bloq_cost.rotations}, clif: {bloq_cost.clifford})")
    print("Complexity of a single ramp:")
    print(f'    T-count:   {ramp_cost.t:g}')
    print(f'    Rotations: {ramp_cost.rotations:g}')
    print(f'    Cliffords: {ramp_cost.clifford:g}')

    # Resource counts for various ramped trotterized unitaries
    dt = 1.0
    n1 = 50
    coef1 = (1.0,)
    build_unitary(problem_instance, coef1, dt, n1)
    n2 = 10
    coef2 = (0.5, 0.5)
    build_unitary(problem_instance, coef2, dt, n2)
    n3 = 3
    s2 = 1.0 / (2.0 * (4 - 4.0**(1/3)))
    coef3 = (s2,s2,s2,s2,0.5-4*s2,0.5-4*s2,s2,s2,s2,s2)
    build_unitary(problem_instance, coef3, dt, n3)

if __name__ == "__main__":
    main()
