from openfermion import InteractionOperator, QubitOperator

def shift_qubit_operator(op, N):
    return sum(
            QubitOperator(" ".join(f"{op[1]}{op[0]+N}" for op in ops), coef)
            for ops, coef in op.terms.items() if coef != 0)

class ZeroAccessor:
    def __getitem__(self, key):
        return 0.0

class MixedFermionBosonOperator:
    def __init__(self,
                 constant, fermionic_one_body_tensor, fermionic_two_body_tensor,
                 bosonic_number_coefficient,
                 interaction_tensor):
        # verify tensor shapes and sizes
        # -- store number of fermionic and bosonic states
        #    -- fermions are one qubit per state
        #    -- bosons may be more than one qubit per state (see "boson counts" below)
        self._Nf = self._verify_fermionic(fermionic_one_body_tensor, fermionic_two_body_tensor)
        self._Nb = self._verify_bosonic(interaction_tensor, self._Nf)
        # store data
        self._fermionic = InteractionOperator(
                constant, fermionic_one_body_tensor, fermionic_two_body_tensor)
        self._bosonic = bosonic_number_coefficient
        self._interaction = interaction_tensor
    def _verify_fermionic(self, tensor1b, tensor2b):
        assert len(tensor1b.shape) == 2
        N = tensor1b.shape[0]
        for L in tensor1b.shape:
            assert L == N
        assert len(tensor2b.shape) == 4
        for L in tensor2b.shape:
            assert L == N
        return N
    def _verify_bosonic(self, tensor, Nf):
        assert len(tensor.shape) == 3
        assert tensor.shape[0] == Nf
        assert tensor.shape[1] == Nf
        Nb = tensor.shape[2]
        return Nb
    def num_fermionic_states(self):
        return self._Nf
    def num_bosonic_states(self):
        return self._Nb
    def max_bosons_per_state(self):
        return self._bencode.maximum_bosons_per_state()
    def qubits_per_bosonic_state(self):
        return self._bencode.number_of_qubits_per_state()
    def total_bosonic_qubits(self):
        return self.qubits_per_bosonic_state() * self.num_bosonic_states()
    def num_qubits(self):
        return self.total_bosonic_qubits() + self.num_fermionic_states()
    def set_fermionic_encoding(self, fermionic_encoder):
        self._fencode = fermionic_encoder
    def set_bosonic_encoding(self, bosonic_encoder):
        self._bencode = bosonic_encoder
    def energy_shift(self, dE):
        self._fermionic = InteractionOperator(
                self._fermionic.constant + dE,
                self._fermionic.one_body_tensor,
                self._fermionic.two_body_tensor)
    def generate_qubit_operator(self):
        # purely-fermionic component
        qo = self._fencode(self._fermionic)
        # interaction component
        for n in range(self.num_bosonic_states()):
           # fermionic piece
           qo_if = self._fencode(InteractionOperator(0, self._interaction[:,:,n], ZeroAccessor()))
           # bosonic piece
           c = self._bencode.creation(n)
           a = self._bencode.annihilation(n)
           qo_ib = shift_qubit_operator(c + a, self.num_fermionic_states())
           # accumulate
           qo += qo_if * qo_ib
        # purely-bosonic component
        for n in range(self.num_bosonic_states()):
           c = self._bencode.creation(n)
           a = self._bencode.annihilation(n)
           qo += self._bosonic * shift_qubit_operator(c * a, self.num_fermionic_states())
        # Eliminate "+0j" (also trims "small" imaginary and real components)
        qo.compress()
        # Return full operator
        return qo

def main():
    from common.bosons_binary import BosonicBinaryEncoding
    import numpy as np
    from openfermion import jordan_wigner
    Nf = 3
    Nb = 2
    fermionic = True
    bosonic = True
    interaction = True
    fermion_constant = 1.0 if fermionic else 0.0
    fermion_one_body = np.ones((Nf,Nf)) if fermionic else np.zeros((Nf,Nf))
    fermion_two_body = np.ones((Nf,Nf,Nf,Nf)) if fermionic else np.zeros((Nf,Nf,Nf,Nf))
    boson_scalar = 1.0 if bosonic else 0.0
    interaction_tensor = np.ones((Nf,Nf,Nb)) if interaction else np.zeros((Nf,Nf,Nb))
    max_bosons_per_state = 2
    op = MixedFermionBosonOperator(
            fermion_constant, fermion_one_body, fermion_two_body,
            boson_scalar, 
            interaction_tensor)
    op.set_fermionic_encoding(jordan_wigner)
    op.set_bosonic_encoding(BosonicBinaryEncoding(max_bosons_per_state))
    print(f"op.num_fermionic_states()       = {op.num_fermionic_states()}")
    assert op.num_fermionic_states() == Nf
    print(f"op.num_bosonic_states()         = {op.num_bosonic_states()}")
    assert op.num_bosonic_states() == Nb
    print(f"op.max_bosons_per_state()       = {op.max_bosons_per_state()}")
    assert op.max_bosons_per_state() == 3
    print(f"op.qubits_per_bosonic_state()   = {op.qubits_per_bosonic_state()}")
    assert op.qubits_per_bosonic_state() == 2
    print(f"op.total_bosonic_qubits()       = {op.total_bosonic_qubits()}")
    assert op.total_bosonic_qubits() == 2 * Nb
    print(f"op.num_qubits()                 = {op.num_qubits()}")
    assert op.num_qubits() == Nf + 2 * Nb
    qo = op.generate_qubit_operator()
    print(f"len(qo.terms)                   = {len(qo.terms)}")
    print(qo)
    print("SUCCESS")

if __name__ == "__main__":
    main()
