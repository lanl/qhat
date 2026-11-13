from abc import ABC, abstractmethod

class BaseTerm(ABC):
    """The abstract base class that defines the generic interface for Hamiltonian terms."""

    @abstractmethod
    def num_qubits(self) -> int:
        pass

    @abstractmethod
    def generate_evolution_operator(self, t: float, hbar: float = 1.0, *, scale: float = 1.0,
                                    **evolution_options):
        """The bloq for exp(-i * scale * term * t / hbar).

        `scale` multiplies the term (e.g. a Trotter ramp coefficient).  `evolution_options` are
        the adjoint and control options of the evolution bloq (`is_adjoint`, `ctrl_values`,
        `multi_control_method`, `controlled_rotation_method`, `omit_control_phase`).
        """

    @abstractmethod
    def has_rotations(self) -> bool:
        """False only if the term is a multiple of the identity (it needs no rotation)."""

    @abstractmethod
    def sign_controlled(self) -> 'BaseTerm':
        """The term on one more qubit, `q[0]`, such that exp(-i * term' * t) is exp(-i * term * t)
        with the opposite sign of t on |0> of `q[0]` and the same on |1> (`-term` times Z on the
        new qubit)."""

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

