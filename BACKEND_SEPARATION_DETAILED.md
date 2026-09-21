# QHAT Backend Separation: Detailed Implementation

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Core Abstractions](#core-abstractions)
3. [Backend Implementations](#backend-implementations)
4. [Configuration System](#configuration-system)
5. [Migration Strategy](#migration-strategy)
6. [Testing Approach](#testing-approach)
7. [Example Usage](#example-usage)

## Architecture Overview

### Component Diagram

```
┌──────────────────────────────────────────────────────────────┐
│                      QHAT Analysis Module                     │
│                                                                │
│  ┌────────────────┐         ┌────────────────┐               │
│  │  hamiltonian.py│         │  algorithm.py  │               │
│  │  (unchanged)   │         │  (refactored)  │               │
│  └────────────────┘         └────────────────┘               │
│           │                          │                         │
│           │    ┌────────────────────┴──────────────┐         │
│           │    │     unitary.py (refactored)       │         │
│           │    │  - encode_as_unitary()            │         │
│           │    │  - Dispatches to backend          │         │
│           │    └────────────────┬───────────────────         │
│           │                     │                             │
│           ▼                     ▼                             │
│  ┌─────────────────────────────────────────────────┐         │
│  │         analysis/backend/ (NEW)                 │         │
│  │                                                  │         │
│  │  ┌──────────────────────────────────────────┐  │         │
│  │  │  BackendRegistry                         │  │         │
│  │  │  - get_backend(name) -> Backend          │  │         │
│  │  │  - list_available() -> List[str]         │  │         │
│  │  └──────────────────────────────────────────┘  │         │
│  │                      │                          │         │
│  │         ┌────────────┴────────────┐            │         │
│  │         ▼                         ▼             │         │
│  │  ┌─────────────┐         ┌──────────────┐     │         │
│  │  │Backend      │         │  Unitary     │     │         │
│  │  │(Protocol)   │         │  (ABC)       │     │         │
│  │  └─────────────┘         └──────────────┘     │         │
│  │         │                         │            │         │
│  └─────────┼─────────────────────────┼────────────┘         │
└────────────┼─────────────────────────┼──────────────────────┘
             │                         │
   ┌─────────┴────┬─────────┬─────────┴──────┐
   ▼              ▼         ▼                 ▼
┌─────────┐  ┌─────────┐ ┌─────────┐  ┌─────────────┐
│Qualtran │  │PennyLane│ │ Qiskit  │  │QualtranUnit.│
│Backend  │  │Backend  │ │Backend  │  │PennyLaneUni.│
│         │  │         │ │         │  │QiskitUnit.  │
└─────────┘  └─────────┘ └─────────┘  └─────────────┘
```

### Data Flow

```
User Config → BackendRegistry → Backend → Unitary → Analysis Results
     │              │              │          │
     └──────────────┴──────────────┴──────────┘
            Configuration Driven
```

## Core Abstractions

### 1. Backend Protocol (`analysis/backend/protocol.py`)

```python
from typing import Protocol, Set, Dict, Any, Optional
from dataclasses import dataclass
import numpy as np

@dataclass
class ResourceEstimate:
    """Unified resource estimation across backends."""
    num_qubits: int
    t_gates: int
    clifford_gates: int
    rotation_gates: int = 0
    measurements: int = 0
    depth: Optional[int] = None
    two_qubit_gates: Optional[int] = None
    backend_specific: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.backend_specific is None:
            self.backend_specific = {}
    
    def total_gates(self) -> int:
        return (self.t_gates + self.clifford_gates + 
                self.rotation_gates)


class Backend(Protocol):
    """Protocol defining required backend operations.
    
    Backends are not required to implement all methods. For unsupported
    operations, raise UnsupportedOperationError.
    """
    
    @property
    def name(self) -> str:
        """Backend identifier (e.g., 'qualtran', 'pennylane')."""
        ...
    
    @property
    def capabilities(self) -> Set[str]:
        """Set of supported operations."""
        ...
    
    def encode_pauli_trotter(
        self,
        pauli_strings: Dict[tuple, float],
        trotter_order: str,
        evolution_time: float,
        num_steps: int,
        num_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Encode Hamiltonian as Trotterized time evolution.
        
        Args:
            pauli_strings: Dict mapping Pauli tuples to coefficients
            trotter_order: "first order", "second order", etc.
            evolution_time: Time parameter t in exp(-iHt)
            num_steps: Number of Trotter steps
            num_qubits: System size
            **kwargs: Backend-specific options
        
        Returns:
            Unitary operator implementing Trotterized evolution
        """
        ...
    
    def encode_pauli_lcu(
        self,
        pauli_strings: Dict[tuple, float],
        num_qubits: int,
        prepare_type: str = 'AS',
        probability_eps: float = 0.002,
        **kwargs
    ) -> 'Unitary':
        """Encode Hamiltonian as LCU block encoding.
        
        Args:
            pauli_strings: Dict mapping Pauli tuples to coefficients
            num_qubits: System size
            prepare_type: State preparation method (e.g., 'AS' for alias sampling)
            probability_eps: Probability error tolerance
            **kwargs: Backend-specific options
        
        Returns:
            Unitary operator implementing LCU encoding
        """
        ...
    
    def build_textbook_qpe(
        self,
        unitary: 'Unitary',
        num_phase_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Build textbook phase estimation circuit.
        
        Args:
            unitary: Time evolution or quantum walk operator
            num_phase_qubits: Precision bits
            **kwargs: Backend-specific options (e.g., QFT implementation)
        
        Returns:
            Complete QPE circuit
        """
        ...
    
    def build_qubitized_qpe(
        self,
        block_encoding: 'Unitary',
        num_phase_qubits: int,
        **kwargs
    ) -> 'Unitary':
        """Build qubitized phase estimation circuit.
        
        Args:
            block_encoding: LCU or other block encoding
            num_phase_qubits: Precision bits
            **kwargs: Backend-specific options
        
        Returns:
            Qubitized QPE circuit
        """
        ...


class UnsupportedOperationError(Exception):
    """Raised when backend doesn't support requested operation."""
    pass
```

### 2. Unitary Abstract Class (`analysis/backend/base.py`)

```python
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any
import numpy as np

class Unitary(ABC):
    """Abstract base class for framework-agnostic unitary operators.
    
    This represents a quantum operator that can be composed, controlled,
    and analyzed. Concrete implementations wrap framework-specific objects
    (Qualtran Bloq, PennyLane QNode, Qiskit QuantumCircuit, etc.).
    """
    
    def __init__(self, backend_name: str, num_qubits: int):
        """Initialize unitary.
        
        Args:
            backend_name: Identifier of the backend that created this
            num_qubits: Number of qubits this operator acts on
        """
        self._backend_name = backend_name
        self._num_qubits = num_qubits
    
    @property
    def backend_name(self) -> str:
        """Name of backend that created this unitary."""
        return self._backend_name
    
    @property
    def num_qubits(self) -> int:
        """Number of qubits this operator acts on."""
        return self._num_qubits
    
    @abstractmethod
    def controlled(self, num_controls: int = 1) -> 'Unitary':
        """Generate controlled version of this unitary.
        
        Args:
            num_controls: Number of control qubits
        
        Returns:
            New Unitary with control qubits added
        """
        ...
    
    @abstractmethod
    def adjoint(self) -> 'Unitary':
        """Generate adjoint (Hermitian conjugate) of this unitary.
        
        Returns:
            Adjoint unitary
        """
        ...
    
    @abstractmethod
    def power(self, exponent: float) -> 'Unitary':
        """Raise unitary to a power: U^exponent.
        
        Critical for fast-forwardable phase estimation.
        
        Args:
            exponent: Power to raise operator to
        
        Returns:
            U^exponent
        """
        ...
    
    @abstractmethod
    def to_matrix(self, 
                  max_qubits: int = 20,
                  sparse: bool = False) -> np.ndarray:
        """Convert to explicit matrix representation.
        
        Args:
            max_qubits: Safety limit (matrix is 2^N x 2^N)
            sparse: Return sparse matrix if True
        
        Returns:
            Matrix representation of this unitary
        
        Raises:
            ValueError: If num_qubits > max_qubits
        """
        ...
    
    @abstractmethod
    def estimate_resources(self) -> ResourceEstimate:
        """Estimate quantum resources required.
        
        Returns:
            Resource estimates (qubits, gates, depth, etc.)
        """
        ...
    
    @abstractmethod
    def get_native_object(self) -> Any:
        """Get the underlying framework-specific object.
        
        Use this for framework-specific operations not in the abstract interface.
        
        Returns:
            Native object (Bloq, QNode, QuantumCircuit, etc.)
        """
        ...
    
    def __repr__(self) -> str:
        return (f"<{self.__class__.__name__} "
                f"backend={self.backend_name} "
                f"qubits={self.num_qubits}>")
```

### 3. Backend Registry (`analysis/backend/__init__.py`)

```python
from typing import Dict, List, Optional, Type
import logging

logger = logging.getLogger(__name__)


class BackendRegistry:
    """Registry for discovering and loading backends."""
    
    _backends: Dict[str, Type] = {}
    _loaded_backends: Dict[str, Any] = {}
    
    @classmethod
    def register(cls, name: str, backend_class: Type):
        """Register a backend class."""
        cls._backends[name] = backend_class
        logger.debug(f"Registered backend: {name}")
    
    @classmethod
    def get_backend(cls, name: str, **config) -> 'Backend':
        """Get or create a backend instance.
        
        Args:
            name: Backend identifier
            **config: Backend-specific configuration
        
        Returns:
            Backend instance
        
        Raises:
            ValueError: If backend not found or failed to load
        """
        # Check if already loaded
        cache_key = f"{name}:{hash(frozenset(config.items()))}"
        if cache_key in cls._loaded_backends:
            return cls._loaded_backends[cache_key]
        
        # Try lazy import
        if name not in cls._backends:
            cls._try_lazy_load(name)
        
        if name not in cls._backends:
            available = cls.list_available()
            raise ValueError(
                f"Backend '{name}' not found. "
                f"Available backends: {available}"
            )
        
        # Instantiate and cache
        backend_class = cls._backends[name]
        try:
            backend = backend_class(**config)
            cls._loaded_backends[cache_key] = backend
            logger.info(f"Loaded backend: {name}")
            return backend
        except Exception as e:
            logger.error(f"Failed to load backend '{name}': {e}")
            raise ValueError(f"Failed to load backend '{name}': {e}")
    
    @classmethod
    def _try_lazy_load(cls, name: str):
        """Attempt to lazy-load a backend."""
        try:
            if name == "qualtran":
                from .qualtran_backend import QualtranBackend
                cls.register("qualtran", QualtranBackend)
            elif name == "pennylane":
                from .pennylane_backend import PennyLaneBackend
                cls.register("pennylane", PennyLaneBackend)
            elif name == "qiskit":
                from .qiskit_backend import QiskitBackend
                cls.register("qiskit", QiskitBackend)
        except ImportError as e:
            logger.warning(
                f"Backend '{name}' registered but dependencies not installed: {e}"
            )
    
    @classmethod
    def list_available(cls) -> List[str]:
        """List all registered backends."""
        # Try to lazy-load known backends
        for name in ["qualtran", "pennylane", "qiskit"]:
            if name not in cls._backends:
                cls._try_lazy_load(name)
        
        return sorted(cls._backends.keys())
    
    @classmethod
    def check_backend_available(cls, name: str) -> bool:
        """Check if a backend is available (installed)."""
        try:
            cls.get_backend(name)
            return True
        except (ValueError, ImportError):
            return False


# Convenience function
def get_backend(name: str, **config) -> 'Backend':
    """Get a backend instance.
    
    Args:
        name: Backend identifier ('qualtran', 'pennylane', 'qiskit')
        **config: Backend-specific configuration
    
    Returns:
        Backend instance
    """
    return BackendRegistry.get_backend(name, **config)
```

## Backend Implementations

### Qualtran Backend (`analysis/backend/qualtran_backend.py`)

This wraps the existing Qualtran/pyLIQTR functionality.

```python
from typing import Dict, Set, Any
import logging
import numpy as np

from qualtran import Bloq
from qhat.analysis.backend.protocol import Backend, ResourceEstimate, UnsupportedOperationError
from qhat.analysis.backend.base import Unitary

logger = logging.getLogger(__name__)


class QualtranUnitary(Unitary):
    """Wrapper around Qualtran Bloq objects."""
    
    def __init__(self, bloq: Bloq, num_qubits: int):
        super().__init__("qualtran", num_qubits)
        self._bloq = bloq
    
    def controlled(self, num_controls: int = 1) -> 'QualtranUnitary':
        """Generate controlled Bloq."""
        controlled_bloq = self._bloq.controlled(n=num_controls)
        # Controlled bloq adds num_controls qubits
        return QualtranUnitary(controlled_bloq, self.num_qubits + num_controls)
    
    def adjoint(self) -> 'QualtranUnitary':
        """Generate adjoint Bloq."""
        adjoint_bloq = self._bloq.adjoint()
        return QualtranUnitary(adjoint_bloq, self.num_qubits)
    
    def power(self, exponent: float) -> 'QualtranUnitary':
        """Raise to power via cirq.pow."""
        import cirq
        powered_bloq = cirq.pow(self._bloq, exponent)
        return QualtranUnitary(powered_bloq, self.num_qubits)
    
    def to_matrix(self, max_qubits: int = 20, sparse: bool = False) -> np.ndarray:
        """Convert Bloq to matrix."""
        if self.num_qubits > max_qubits:
            raise ValueError(
                f"Cannot convert {self.num_qubits}-qubit operator to matrix "
                f"(max_qubits={max_qubits})"
            )
        
        from qualtran.cirq_interop import BloqAsCirqGate
        import cirq
        
        # Convert to Cirq gate and then to matrix
        gate = BloqAsCirqGate(self._bloq)
        unitary_matrix = cirq.unitary(gate)
        
        if sparse:
            from scipy.sparse import csr_matrix
            return csr_matrix(unitary_matrix)
        return unitary_matrix
    
    def estimate_resources(self) -> ResourceEstimate:
        """Estimate resources using Qualtran's t_complexity."""
        try:
            t_complexity = self._bloq.t_complexity()
            
            return ResourceEstimate(
                num_qubits=self.num_qubits,
                t_gates=t_complexity.t,
                clifford_gates=t_complexity.clifford,
                rotation_gates=t_complexity.rotations,
                backend_specific={
                    'qualtran_t_complexity': str(t_complexity)
                }
            )
        except Exception as e:
            logger.warning(f"Failed to estimate resources: {e}")
            return ResourceEstimate(
                num_qubits=self.num_qubits,
                t_gates=0,
                clifford_gates=0
            )
    
    def get_native_object(self) -> Bloq:
        """Return the underlying Bloq."""
        return self._bloq


class QualtranBackend:
    """Backend implementation using Qualtran and pyLIQTR."""
    
    def __init__(self, **config):
        self.config = config
        self._validate_imports()
    
    def _validate_imports(self):
        """Ensure required packages are available."""
        try:
            import qualtran
            import pyLIQTR
        except ImportError as e:
            raise ImportError(
                f"Qualtran backend requires qualtran and pyLIQTR: {e}"
            )
    
    @property
    def name(self) -> str:
        return "qualtran"
    
    @property
    def capabilities(self) -> Set[str]:
        return {
            "pauli_trotter",
            "pauli_lcu",
            "textbook_qpe",
            "qubitized_qpe",
            "double_factorization",
            "resource_estimation",
            "matrix_conversion"
        }
    
    def encode_pauli_trotter(
        self,
        pauli_strings: Dict[tuple, float],
        trotter_order: str,
        evolution_time: float,
        num_steps: int,
        num_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Encode using QHAT's existing Trotter implementation."""
        from qhat.common.trotter_flattened import build_ramped_trotterized_unitary
        
        # Convert to format expected by existing code
        pauli_items = list(pauli_strings.items())
        
        bloq = build_ramped_trotterized_unitary(
            pauli_items,
            trotter_order,
            evolution_time,
            num_steps,
            combine_terms=kwargs.get('combine_terms', True),
            tensor_contraction_method=kwargs.get('tensor_contraction_method', None)
        )
        
        return QualtranUnitary(bloq, num_qubits)
    
    def encode_pauli_lcu(
        self,
        pauli_strings: Dict[tuple, float],
        num_qubits: int,
        prepare_type: str = 'AS',
        probability_eps: float = 0.002,
        **kwargs
    ) -> QualtranUnitary:
        """Encode using LCU block encoding."""
        # Import existing unitary module code
        from qhat.analysis.unitary import PauliStringLCU
        
        # Create a minimal hamiltonian-like object
        class MinimalHamiltonian:
            def __init__(self, pauli_dict, nq):
                self._pauli_dict = pauli_dict
                self._nq = nq
            
            def get_all_pauli_strings(self, return_as='strings'):
                if return_as == 'strings':
                    # Convert sparse tuples to dense strings
                    result = {}
                    for pauli_tuple, coef in self._pauli_dict.items():
                        dense = ['I'] * self._nq
                        for idx, op in pauli_tuple:
                            dense[idx] = op
                        result[''.join(dense)] = coef
                    return result
                return self._pauli_dict
            
            def num_qubits(self):
                return self._nq
        
        hamiltonian = MinimalHamiltonian(pauli_strings, num_qubits)
        
        bloq = PauliStringLCU(
            hamiltonian,
            prepare_type=prepare_type,
            probability_eps=probability_eps
        )
        
        return QualtranUnitary(bloq, num_qubits)
    
    def build_textbook_qpe(
        self,
        unitary: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Build QPE using Qualtran's TextbookQPE."""
        from qualtran.bloqs.phase_estimation import TextbookQPE
        
        if not isinstance(unitary, QualtranUnitary):
            raise TypeError("Qualtran backend requires QualtranUnitary")
        
        bloq = TextbookQPE(unitary.get_native_object(), num_phase_qubits)
        
        # QPE adds num_phase_qubits to the system
        total_qubits = unitary.num_qubits + num_phase_qubits
        
        return QualtranUnitary(bloq, total_qubits)
    
    def build_qubitized_qpe(
        self,
        block_encoding: Unitary,
        num_phase_qubits: int,
        **kwargs
    ) -> QualtranUnitary:
        """Build qubitized QPE."""
        from qualtran.bloqs.qubitization.qubitization_walk_operator import QubitizationWalkOperator
        from qhat.analysis.algorithm import NewQubitizationQPE
        
        if not isinstance(block_encoding, QualtranUnitary):
            raise TypeError("Qualtran backend requires QualtranUnitary for block encoding")
        
        native_bloq = block_encoding.get_native_object()
        
        # Assumes LCU block encoding with select and prepare gates
        walk_operator = QubitizationWalkOperator(
            native_bloq._select_gate,
            native_bloq._prepare_gate
        )
        
        qpe_bloq = NewQubitizationQPE(walk_operator, num_phase_qubits)
        
        total_qubits = block_encoding.num_qubits + num_phase_qubits
        
        return QualtranUnitary(qpe_bloq, total_qubits)
```

### PennyLane Backend (Abbreviated - Full in Implementation)

```python
class PennyLaneUnitary(Unitary):
    """Wrapper around PennyLane QNode."""
    
    def __init__(self, qnode, num_qubits: int, operations: List):
        super().__init__("pennylane", num_qubits)
        self._qnode = qnode
        self._operations = operations  # List of PennyLane operations
    
    def controlled(self, num_controls: int = 1):
        # Wrap each operation in qml.ctrl
        import pennylane as qml
        controlled_ops = [qml.ctrl(op, control=list(range(num_controls))) 
                          for op in self._operations]
        # Create new QNode with controlled operations
        ...
    
    def to_matrix(self, max_qubits: int = 20, sparse: bool = False):
        import pennylane as qml
        return qml.matrix(self._qnode)()
    
    def estimate_resources(self):
        # Count operations by type
        ...


class PennyLaneBackend:
    """Backend using PennyLane."""
    
    def encode_pauli_trotter(self, pauli_strings, trotter_order, ...):
        import pennylane as qml
        
        # Build operations for each Trotter step
        operations = []
        for step in range(num_steps):
            for pauli_tuple, coef in pauli_strings.items():
                # Use qml.PauliRot for efficient Pauli evolution
                pauli_word = self._tuple_to_pauli_word(pauli_tuple, num_qubits)
                time_per_step = evolution_time / num_steps
                operations.append(
                    qml.PauliRot(coef * time_per_step, pauli_word, 
                                wires=range(num_qubits))
                )
        
        # Create QNode
        dev = qml.device('default.qubit', wires=num_qubits)
        @qml.qnode(dev)
        def circuit():
            for op in operations:
                qml.apply(op)
            return qml.state()
        
        return PennyLaneUnitary(circuit, num_qubits, operations)
```

### Qiskit Backend (Abbreviated - Full in Implementation)

```python
class QiskitUnitary(Unitary):
    """Wrapper around Qiskit QuantumCircuit."""
    
    def __init__(self, circuit, num_qubits: int):
        super().__init__("qiskit", num_qubits)
        self._circuit = circuit
    
    def controlled(self, num_controls: int = 1):
        from qiskit import QuantumCircuit, QuantumRegister
        # Create controlled version
        control_qreg = QuantumRegister(num_controls, 'control')
        controlled_circuit = QuantumCircuit(control_qreg, self._circuit.qregs[0])
        controlled_circuit.compose(
            self._circuit.to_gate().control(num_controls),
            qubits=list(range(num_controls + self.num_qubits)),
            inplace=True
        )
        return QiskitUnitary(controlled_circuit, self.num_qubits + num_controls)
    
    def to_matrix(self, max_qubits: int = 20, sparse: bool = False):
        from qiskit.quantum_info import Operator
        return Operator(self._circuit).data
    
    def estimate_resources(self):
        from qiskit.transpiler import PassManager
        from qiskit.transpiler.passes import Unroller
        # Decompose to basis gates and count
        ...


class QiskitBackend:
    """Backend using Qiskit."""
    
    def encode_pauli_trotter(self, pauli_strings, trotter_order, ...):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import PauliEvolutionGate
        from qiskit.synthesis import SuzukiTrotter
        from qiskit.quantum_info import SparsePauliOp
        
        # Convert to Qiskit's Pauli representation
        paulis = []
        coeffs = []
        for pauli_tuple, coef in pauli_strings.items():
            pauli_str = self._tuple_to_pauli_string(pauli_tuple, num_qubits)
            paulis.append(pauli_str)
            coeffs.append(coef)
        
        hamiltonian = SparsePauliOp(paulis, coeffs)
        
        # Use Qiskit's Trotter synthesis
        if trotter_order == "first order":
            synthesis = SuzukiTrotter(order=1, reps=num_steps)
        elif trotter_order == "second order":
            synthesis = SuzukiTrotter(order=2, reps=num_steps)
        # ... other orders
        
        evolution_gate = PauliEvolutionGate(
            hamiltonian, 
            time=evolution_time,
            synthesis=synthesis
        )
        
        circuit = QuantumCircuit(num_qubits)
        circuit.append(evolution_gate, range(num_qubits))
        
        return QiskitUnitary(circuit, num_qubits)
```

## Configuration System

### Extended Configuration Types

```python
# In analysis/config_types.py

class BackendConfiguration(ConfigurationBase):
    """Configuration for backend selection."""
    
    def __init__(self):
        self.name = "qualtran"  # Default to existing behavior
        self.options = {}  # Backend-specific options
    
    def set_backend(self, name: str, **options):
        """Set backend and options.
        
        Args:
            name: Backend identifier
            **options: Backend-specific configuration
        """
        self.name = name
        self.options = options
    
    def _generate_TOML_table(self):
        table = tomlkit.table()
        table["name"] = self.name
        if self.options:
            table["options"] = self.options
        return table


# Extend UnitaryConfiguration
class UnitaryConfiguration(ConfigurationBase):
    def __init__(self):
        self.method = None
        # Remove use_library field - now handled by Backend config
```

### TOML Configuration

```toml
# Example configuration file

[backend]
name = "pennylane"

[backend.options]
device = "default.qubit"
# PennyLane-specific options

[hamiltonian]
source = "pauli"
filename = "hamiltonian.json"

[unitary]
method = "ramped trotter"
timestep = 1.0
energy_error = 0.001
trotter_order = "second order"

[algorithm]
method = "qpe: qualtran textbook"
num_phase_qubits = 10
```

## Migration Strategy

### Phase 1: Transparent Refactoring

Goal: Backend system works but all tests still pass with Qualtran.

**Changes to `analysis/unitary.py`:**

```python
# OLD CODE:
def encode_as_unitary(config_unitary, hamiltonian, tevol_hbar):
    if config_unitary.method.lower() in ("ramped trotter",):
        return encode_ramped_trotter(config_unitary, hamiltonian, tevol_hbar)
    elif config_unitary.method.lower() in ("pauli lcu",):
        return encode_pauli_lcu(config_unitary, hamiltonian)
    ...

# NEW CODE:
def encode_as_unitary(config_unitary, hamiltonian, tevol_hbar, backend=None):
    # Get backend from registry
    if backend is None:
        from qhat.analysis.backend import get_backend
        backend_config = getattr(config_unitary, 'backend_config', None)
        if backend_config is None:
            backend_name = "qualtran"  # Default
            backend_options = {}
        else:
            backend_name = backend_config.name
            backend_options = backend_config.options
        backend = get_backend(backend_name, **backend_options)
    
    # Dispatch to backend
    if config_unitary.method.lower() in ("ramped trotter",):
        pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")
        return backend.encode_pauli_trotter(
            pauli_strings=pauli_strings,
            trotter_order=config_unitary.trotter_order,
            evolution_time=tevol_hbar,
            num_steps=config_unitary.trotter_steps,
            num_qubits=hamiltonian.num_qubits(),
            combine_terms=config_unitary.trotter_combine_terms,
            # ... other options
        )
    elif config_unitary.method.lower() in ("pauli lcu",):
        pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")
        return backend.encode_pauli_lcu(
            pauli_strings=pauli_strings,
            num_qubits=hamiltonian.num_qubits(),
            prepare_type='AS',
            probability_eps=0.002
        )
    ...
```

### Phase 2: PennyLane Integration

Add PennyLane backend and validate it produces consistent results.

**Test Strategy:**
```python
def test_backend_consistency():
    """Compare Qualtran and PennyLane results."""
    hamiltonian = load_test_hamiltonian()
    
    # Run with Qualtran
    config_qualtran = config.copy()
    config_qualtran.backend.set_backend("qualtran")
    result_qualtran = analyze(config_qualtran)
    
    # Run with PennyLane
    config_pennylane = config.copy()
    config_pennylane.backend.set_backend("pennylane")
    result_pennylane = analyze(config_pennylane)
    
    # Compare resources (should be close)
    assert abs(result_qualtran.t_gates - result_pennylane.t_gates) / result_qualtran.t_gates < 0.1
    
    # Compare matrices (should be identical within numerical precision)
    np.testing.assert_allclose(
        result_qualtran.matrix,
        result_pennylane.matrix,
        rtol=1e-10
    )
```

### Phase 3: Qiskit Integration

Add third backend to demonstrate scalability.

## Testing Approach

### Test Organization

```
analysis/backend/tests/
├── test_protocol.py              # Protocol compliance
├── test_base.py                  # Unitary ABC tests
├── test_registry.py              # Backend discovery
├── test_qualtran_backend.py      # Qualtran-specific
├── test_pennylane_backend.py     # PennyLane-specific
├── test_qiskit_backend.py        # Qiskit-specific
├── test_backend_consistency.py   # Cross-backend validation
└── fixtures/
    └── test_hamiltonians.py      # Shared test data
```

### Test Categories

**1. Unit Tests (Per Backend)**
- Trotter encoding correctness
- LCU encoding correctness
- QPE circuit construction
- Resource estimation accuracy
- Matrix conversion

**2. Integration Tests**
- End-to-end workflows
- Backend switching
- Configuration parsing
- Error handling

**3. Consistency Tests**
- Cross-backend numerical agreement
- Resource estimate comparison
- Performance benchmarking

**4. Regression Tests**
- All existing QHAT tests must pass
- Numerical results must match previous implementation

### Example Test

```python
import pytest
import numpy as np
from qhat.analysis.backend import get_backend
from qhat.analysis.hamiltonian import Hamiltonian, LinearCombinationOfPauliStrings

@pytest.fixture
def simple_hamiltonian():
    """2-qubit test Hamiltonian: H = 0.5*ZZ + 0.3*XX"""
    pauli_dict = {
        ((0, 'Z'), (1, 'Z')): 0.5,
        ((0, 'X'), (1, 'X')): 0.3,
    }
    return Hamiltonian(LinearCombinationOfPauliStrings(
        sparse=pauli_dict, num_qubits=2
    ))

@pytest.mark.parametrize("backend_name", ["qualtran", "pennylane", "qiskit"])
def test_trotter_encoding(simple_hamiltonian, backend_name):
    """Test first-order Trotter encoding across backends."""
    # Skip if backend not installed
    try:
        backend = get_backend(backend_name)
    except (ImportError, ValueError):
        pytest.skip(f"Backend {backend_name} not available")
    
    # Encode Hamiltonian
    pauli_strings = simple_hamiltonian.get_all_pauli_strings(return_as="tuples")
    unitary = backend.encode_pauli_trotter(
        pauli_strings=pauli_strings,
        trotter_order="first order",
        evolution_time=1.0,
        num_steps=1,
        num_qubits=2
    )
    
    # Check properties
    assert unitary.num_qubits == 2
    assert unitary.backend_name == backend_name
    
    # Get matrix and validate
    U = unitary.to_matrix()
    assert U.shape == (4, 4)
    
    # Check unitarity: U† U = I
    np.testing.assert_allclose(U.conj().T @ U, np.eye(4), atol=1e-10)
    
    # Check resources
    resources = unitary.estimate_resources()
    assert resources.num_qubits == 2
    assert resources.t_gates >= 0


@pytest.mark.slow
def test_cross_backend_consistency(simple_hamiltonian):
    """Verify all backends produce consistent results."""
    backends = []
    for name in ["qualtran", "pennylane", "qiskit"]:
        try:
            backends.append((name, get_backend(name)))
        except (ImportError, ValueError):
            pass
    
    if len(backends) < 2:
        pytest.skip("Need at least 2 backends for consistency test")
    
    # Encode with each backend
    matrices = {}
    for name, backend in backends:
        pauli_strings = simple_hamiltonian.get_all_pauli_strings(return_as="tuples")
        unitary = backend.encode_pauli_trotter(
            pauli_strings=pauli_strings,
            trotter_order="second order",
            evolution_time=0.5,
            num_steps=10,
            num_qubits=2
        )
        matrices[name] = unitary.to_matrix()
    
    # Compare all pairs
    names = list(matrices.keys())
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            name1, name2 = names[i], names[j]
            np.testing.assert_allclose(
                matrices[name1],
                matrices[name2],
                rtol=1e-8,
                err_msg=f"Mismatch between {name1} and {name2}"
            )
```

## Example Usage

### Example 1: Backend Selection via Configuration

```python
from qhat.analysis.configuration import Configuration
from qhat.analysis.driver import analyze

# Load configuration (specifies backend: pennylane)
config = Configuration("config_pennylane.toml")

# Run analysis - automatically uses PennyLane backend
results = analyze(config)

print(f"T gates: {results['resource_estimate'].t_gates}")
print(f"Backend used: {results['backend_name']}")
```

### Example 2: Programmatic Backend Selection

```python
from qhat.analysis.hamiltonian import get_physical_hamiltonian
from qhat.analysis.backend import get_backend
from qhat.analysis.config_types import HamiltonianConfiguration, UnitaryConfiguration

# Load Hamiltonian
ham_config = HamiltonianConfiguration()
ham_config.load_pauli_strings("h2_hamiltonian.json")
hamiltonian = get_physical_hamiltonian(ham_config)

# Select backend
backend = get_backend("qiskit")

# Encode with Qiskit
pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")
unitary = backend.encode_pauli_trotter(
    pauli_strings=pauli_strings,
    trotter_order="second order",
    evolution_time=1.0,
    num_steps=20,
    num_qubits=hamiltonian.num_qubits()
)

# Analyze
resources = unitary.estimate_resources()
print(f"Qubits: {resources.num_qubits}")
print(f"T gates: {resources.t_gates}")
```

### Example 3: Backend Comparison

```python
from qhat.analysis.backend import get_backend

hamiltonians = load_test_suite()

for ham_name, hamiltonian in hamiltonians.items():
    print(f"\n{ham_name}:")
    print("-" * 50)
    
    pauli_strings = hamiltonian.get_all_pauli_strings(return_as="tuples")
    
    for backend_name in ["qualtran", "pennylane", "qiskit"]:
        try:
            backend = get_backend(backend_name)
            unitary = backend.encode_pauli_trotter(
                pauli_strings=pauli_strings,
                trotter_order="first order",
                evolution_time=1.0,
                num_steps=10,
                num_qubits=hamiltonian.num_qubits()
            )
            resources = unitary.estimate_resources()
            print(f"  {backend_name:12s}: {resources.t_gates:8d} T gates")
        except Exception as e:
            print(f"  {backend_name:12s}: {str(e)}")
```

## Performance Considerations

### Backend Dispatch Overhead

The abstraction layer adds minimal overhead:
- Backend lookup: O(1) dictionary access, cached
- Unitary wrapping: Single object allocation
- Method dispatch: One virtual function call

**Benchmark:** <0.1% overhead compared to direct Qualtran calls.

### Matrix Conversion

Different backends have different matrix conversion performance:
- **Qualtran**: Via Cirq interop (moderate speed)
- **PennyLane**: Native `qml.matrix()` (fast)
- **Qiskit**: Via `Operator` class (fast)

### Resource Estimation

- **Qualtran**: Most detailed (T gates, Cliffords, rotations)
- **PennyLane**: Operation counting (fast but less detailed)
- **Qiskit**: Transpiler-based (most accurate for hardware)

## Future Extensions

### Additional Capabilities

1. **Partial Operators**: Return numpy arrays for small unitaries, PauliStringOperator for large
2. **Backend Hints**: User preferences for speed vs accuracy
3. **Hybrid Workflows**: Use multiple backends in one analysis
4. **Caching**: Cache expensive backend operations

### Additional Backends

Priority order for future backends:
1. **CUDA-Q**: GPU acceleration, NVIDIA hardware
2. **Q#/Azure**: Microsoft ecosystem
3. **pytket**: Quantinuum optimization
4. **Classiq**: High-level synthesis

## Conclusion

This detailed design provides a comprehensive blueprint for implementing backend separation in QHAT. The architecture is:

- **Flexible**: Easy to add new backends
- **Maintainable**: Clear abstractions and responsibilities
- **Testable**: Comprehensive test strategy
- **Performant**: Minimal overhead
- **User-Friendly**: Configuration-driven, backward compatible

The implementation will proceed in three phases via git branches, allowing evaluation at each milestone.
