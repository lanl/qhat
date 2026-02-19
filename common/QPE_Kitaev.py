"""
An implementation of the "standard" quantum phase estimation circuit, originally introduced by
Alexei Kitaev in https://arxiv.org/abs/quant-ph/9511026.  Also discussed and illustrated at
https://en.wikipedia.org/wiki/Quantum_phase_estimation_algorithm#cite_note-kitaev-1.
"""

# TODO: Compare qualtran.bloqs.phase_estimation.TextbookQPE.  It may be that this is an interesting
#       exercise to learn about building things in Qualtran, but not something we'll actually use.
#       -- Note that in the same directory Qualtran has a variety of phase estimation circuits that
#          we could use in addition to / instead of what pyLIQTR provides.

from attrs import frozen
from qualtran import Bloq, BloqBuilder, QBit, Register, Signature, SoquetT
from qualtran.bloqs.basic_gates import Hadamard
from qualtran.bloqs.qft.qft_text_book import QFTTextBook

@frozen
class QPE_Kitaev(Bloq):
    U: Bloq
    num_phase_qubits: int
    # TODO: Can we automatically read this from the unitary?
    num_state_qubits: int

    @property
    def signature(self):
        return Signature.build(phase=self.num_phase_qubits, state=self.num_state_qubits)

    def build_composite_bloq(
            self, bb: BloqBuilder, *, phase: SoquetT, state: SoquetT) -> dict[str, SoquetT]:

        # Split the phase register to access individual qubits
        phase_s = bb.split(phase)

        # Wall of Hadamard gates on the phase register
        for i in range(self.num_phase_qubits):
            phase_s[i] = bb.add(Hadamard(), q=phase_s[i])

        # Controlled-unitary powers on state register, controlled by phase register bits
        for i in range(self.num_phase_qubits):
            # TODO: How do we know that 'x' will be the correct name?
            #       -- See unitary.py and its use of `add_and_partition`?
            phase_s[i], state = bb.add(cirq.pow(U, 2**i).controlled(), ctrl=phase_s[i], x=state)

        # Merge the individual phase qubits back into a single register
        phase = bb.join(phase_s)

        # Inverse Fourier transform on the phase register
        phase = bb.add(QFTTextBook(bitsize=self.num_phase_qubits).adjoint(), q=phase)

        return {
                'phase': phase,
                'state': state,
               }

if __name__ == "__main__":

    # I just need some unitary for testing
    from qualtran.bloqs.rotations import QvrZPow
    U = QvrZPow.from_bitsize(5, gamma=0.1, eps=1e-2)

    mybloq0 = QPE_Kitaev(U, 3, 5)
    mybloq1 = mybloq0.decompose_bloq()

    # Note that the below appears to only work in a Jupyter notebook
    draw_score = True
    if draw_score:
        from qualtran.drawing import get_musical_score_data, draw_musical_score 
        msd = get_musical_score_data(mybloq1)
        fig, ax = draw_musical_score(msd)
        fig.tight_layout()
    else:
        from qualtran.drawing import show_bloq
        show_bloq(mybloq1)

