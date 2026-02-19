import cirq
import math
from openfermion import QubitOperator

# =================================================================================================
# This block only for development and testing.  Remove before putting into production.

def outer_term(num_qubits, ket, bra, coef):
    outers = {
        (0,0): lambda j: QubitOperator("", 0.5) + QubitOperator(f"Z{j}", 0.5),
        (0,1): lambda j: QubitOperator(f"X{j}", 0.5) + QubitOperator(f"Y{j}", 0.5j),
        (1,0): lambda j: QubitOperator(f"X{j}", 0.5) - QubitOperator(f"Y{j}", 0.5j),
        (1,1): lambda j: QubitOperator("", 0.5) - QubitOperator(f"Z{j}", 0.5),
    }
    bit = lambda n, j : (n // 2**j) % 2
    outer1 = lambda d_m, d_n, j : outers[(d_m, d_n)](j)
    # TODO: Check if range should be reversed(range)
    return coef * math.prod(outer1(bit(ket, j), bit(bra, j), num_qubits-j-1) for j in range(num_qubits))

def beta_dagger(num_qubits):
    return sum(outer_term(num_qubits, n+1, n, math.sqrt(n+1)) for n in range(2**num_qubits-1))

def beta(num_qubits):
    return sum(outer_term(num_qubits, n-1, n, math.sqrt(n)) for n in range(1, 2**num_qubits))

def qubits_per_bosonic_state(max_num_bosons):
    return math.ceil(math.log2(max_num_bosons + 1))

def shift(op, N):
    return sum(
            QubitOperator(" ".join(f"{op[1]}{op[0]+N}" for op in ops), coef)
            for ops, coef in op.terms.items() if coef != 0)

def creation(max_num_bosons, target_state):
    Q = qubits_per_bosonic_state(max_num_bosons)
    op = beta_dagger(Q)
    return shift(op, target_state * Q)

def annihilation(max_num_bosons, target_state):
    Q = qubits_per_bosonic_state(max_num_bosons)
    op = beta(Q)
    return shift(op, target_state * Q)

def matrix_from_Pauli_strings(pauli_strings, num_qubits):
    q = cirq.LineQubit.range(num_qubits)
    op = {'I': cirq.I, 'X': cirq.X, 'Y': cirq.Y, 'Z': cirq.Z}
    hamiltonian_gate = sum(cirq.PauliString((op[p](q[i]) for i, p in s), coefficient=c)
                           for s, c in pauli_strings.items())
    return hamiltonian_gate.matrix()

# =================================================================================================
# This block is the production implementation

class BosonicBinaryEncoding:
    def __init__(self, max_bosons_per_state_user):
        # number of qubits per state
        self._Q = math.ceil(math.log2(max_bosons_per_state_user + 1))
        # maximum number of bosons per state (greater than or equal to user-provided value)
        self._B = 2**self._Q - 1
    def _outer_term(self, ket, bra, coef):
        I = lambda j, coef : QubitOperator("", coef)
        X = lambda j, coef : QubitOperator(f"X{j}", coef)
        Y = lambda j, coef : QubitOperator(f"Y{j}", coef)
        Z = lambda j, coef : QubitOperator(f"Z{j}", coef)
        outers = {
            (0,0): lambda j: I(j, 0.5) + Z(j, 0.5),
            (0,1): lambda j: X(j, 0.5) + Y(j, 0.5j),
            (1,0): lambda j: X(j, 0.5) - Y(j, 0.5j),
            (1,1): lambda j: I(j, 0.5) - Z(j, 0.5),
        }
        def outer1(ket, bra, j):
            bit = lambda n, j : (n // 2**j) % 2
            d_m = bit(ket, j)     # bit j of number in |ket>
            d_n = bit(bra, j)     # bit j of number in <bra|
            key = (d_m, d_n)      # pair bits together as a key
            q = self._Q - j - 1   # put on qubit (order reversed)
            return outers[key](q) # build the appropriate QubitOperator
        return coef * math.prod(outer1(ket, bra, j) for j in range(self._Q))
    def _beta_dagger(self):
        """The bosonic creation operator on a single state."""
        return sum(self._outer_term(n+1, n, math.sqrt(n+1)) for n in range(self._B))
    def _beta(self):
        """The bosonic annihilation operator on a single state."""
        return sum(self._outer_term(n-1, n, math.sqrt(n)) for n in range(1, self._B + 1))
    def maximum_bosons_per_state(self):
        """The maximum number of bosons that can be in a single state.  Greater than or equal to
        the user-requested value."""
        return self._B
    def number_of_qubits_per_state(self):
        """The number of qubits used to encode the count in a single bosonic state."""
        return self._Q
    def creation(self, target_state):
        """The bosonic creation operator for multiple states, incrementing the boson count in
        target_state.  Incrementing the boson count above the maximum number of bosons permitted by
        this encoding results in a null state."""
        return shift(self._beta_dagger(), target_state * self._Q)
    def annihilation(self, target_state):
        """The bosonic annihilation operator for multiple states, decrementing the boson count in
        target_state.  Decrementing the boson count below zero results in a null state."""
        return shift(self._beta(), target_state * self._Q)
    def number(self, target_state):
        """The number operator for multiple states, querying the boson count in target_state."""
        return self.creation(target_state) * self.annihilation(target_state)

# =================================================================================================
# This block is only for testing

def main():
    number_of_bosonic_states = 3
    maximum_number_of_bosons_per_state = 3
    total_qubits = number_of_bosonic_states * \
            qubits_per_bosonic_state(maximum_number_of_bosons_per_state)
    nn = dict()
    for i in range(number_of_bosonic_states):
        c = creation(maximum_number_of_bosons_per_state, i)
        a = annihilation(maximum_number_of_bosons_per_state, i)
        n = c * a
        nn[i] = n
    n = sum(nn.values())
    N1 = matrix_from_Pauli_strings(n.terms, total_qubits)
    print(N1.shape)
    for i in range(N1.shape[0]):
        print("[", end='')
        for j in range(N1.shape[1]):
            x = N1[i,j]
            assert x.imag == 0
            x = x.real
            if abs(x) <= 1.0e-9:
                print(f" {'':4}", end='')
            else:
                print(f" {x:4.1f}", end='')
        print(" ]")

    encoding = BosonicBinaryEncoding(maximum_number_of_bosons_per_state)
    n = sum(encoding.number(i) for i in range(number_of_bosonic_states))
    N2 = matrix_from_Pauli_strings(n.terms, total_qubits)

    for i in range(number_of_bosonic_states):
        assert encoding.creation(i) == creation(maximum_number_of_bosons_per_state, i), \
                f"Creation operator {i} DIFFERS."
    print("Creation operators all match.")

    for i in range(number_of_bosonic_states):
        assert encoding.annihilation(i) == annihilation(maximum_number_of_bosons_per_state, i), \
                f"Annihilation operator {i} DIFFERS."
    print("Annihilation operators all match.")

    assert (N1 == N2).all(), "Number-sum matrices DIFFER."
    print("Number-sum matrices match.")

if __name__ == "__main__":
    main()
