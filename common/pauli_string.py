"""Representation-independent Pauli string value type."""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from numbers import Integral
from typing import TypeAlias


SparsePauli: TypeAlias = tuple[tuple[int, str], ...]
SparsePauliInput: TypeAlias = Mapping[int, str] | Iterable[tuple[int, str]]

_VALID_OPERATORS = frozenset("IXYZ")


def _validate_num_qubits(num_qubits: int) -> int:
    if isinstance(num_qubits, bool) or not isinstance(num_qubits, Integral):
        raise TypeError("num_qubits must be an integer.")
    if num_qubits < 0:
        raise ValueError("num_qubits must be non-negative.")
    return int(num_qubits)


def _normalize_sparse(operators: SparsePauliInput, num_qubits: int) -> SparsePauli:
    if isinstance(operators, Mapping):
        items = operators.items()
    elif isinstance(operators, (str, bytes)) or not isinstance(operators, Iterable):
        raise TypeError("Sparse Pauli input must be a mapping or an iterable of pairs.")
    else:
        items = operators

    normalized = []
    seen_indices = set()
    for item in items:
        try:
            index, operator = item
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "Each sparse Pauli operator must be an (index, operator) pair."
            ) from exc

        if isinstance(index, bool) or not isinstance(index, Integral):
            raise TypeError("Pauli operator indices must be integers.")
        index = int(index)
        if index < 0 or index >= num_qubits:
            raise ValueError(
                f"Pauli operator index {index} is outside the range "
                f"[0, {num_qubits})."
            )
        if index in seen_indices:
            raise ValueError(f"Duplicate Pauli operator index: {index}.")
        seen_indices.add(index)

        if not isinstance(operator, str) or len(operator) != 1:
            raise TypeError("Pauli operators must be single-character strings.")
        if operator not in _VALID_OPERATORS:
            raise ValueError(
                f"Invalid Pauli operator {operator!r}; expected one of I, X, Y, or Z."
            )
        if operator != "I":
            normalized.append((index, operator))

    return tuple(sorted(normalized))


@dataclass(frozen=True, slots=True)
class PauliString:
    """An immutable Pauli string with dense and sparse representations.

    Coefficients are intentionally kept outside this class so instances can be
    used as keys in Hamiltonian dictionaries.
    """

    num_qubits: int
    _operators: SparsePauli

    def __post_init__(self) -> None:
        num_qubits = _validate_num_qubits(self.num_qubits)
        operators = _normalize_sparse(self._operators, num_qubits)
        object.__setattr__(self, "num_qubits", num_qubits)
        object.__setattr__(self, "_operators", operators)

    @classmethod
    def from_dense(cls, value: str) -> "PauliString":
        """Construct a Pauli string from a dense value such as ``"IXYZ"``."""
        if not isinstance(value, str):
            raise TypeError("Dense Pauli input must be a string.")
        invalid = next(
            (operator for operator in value if operator not in _VALID_OPERATORS),
            None,
        )
        if invalid is not None:
            raise ValueError(
                f"Invalid Pauli operator {invalid!r}; expected one of I, X, Y, or Z."
            )

        operators = tuple(
            (index, operator)
            for index, operator in enumerate(value)
            if operator != "I"
        )
        return cls(num_qubits=len(value), _operators=operators)

    @classmethod
    def from_sparse(
        cls,
        operators: SparsePauliInput,
        num_qubits: int,
    ) -> "PauliString":
        """Construct from sparse ``(index, operator)`` pairs or a mapping."""
        num_qubits = _validate_num_qubits(num_qubits)
        return cls(
            num_qubits=num_qubits,
            _operators=_normalize_sparse(operators, num_qubits),
        )

    def to_dense(self) -> str:
        """Return the full dense string, including identity operators."""
        dense = ["I"] * self.num_qubits
        for index, operator in self._operators:
            dense[index] = operator
        return "".join(dense)

    def to_sparse(self) -> SparsePauli:
        """Return the canonical sparse tuple used by existing QHAT APIs."""
        return self._operators

    def to_sparse_dict(self) -> dict[int, str]:
        """Return the sparse dictionary representation described in issue #36."""
        return dict(self._operators)

    def __len__(self) -> int:
        return self.num_qubits

    def __str__(self) -> str:
        return self.to_dense()

    def __repr__(self) -> str:
        return f"PauliString.from_dense({self.to_dense()!r})"
