from abc import ABC, abstractmethod

class BaseTerm(ABC):
    """The abstract base class that defines the generic interface for Hamiltonian terms."""

    @abstractmethod
    def num_qubits(self) -> int:
        pass

    @abstractmethod
    def generate_evolution_operator(self, t: float, hbar: float = 1.0):
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

