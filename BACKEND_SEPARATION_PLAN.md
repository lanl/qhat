# QHAT Backend Separation: High-Level Plan

## Executive Summary

This document outlines the implementation of QHAT's separation into a **front end** (user-facing logic and unified framework) and a **back end** (dispatching to quantum computing frameworks). This branch will fully implement the abstraction layer along with three concrete backends: **Qualtran** (existing), **PennyLane** (new), and **Qiskit** (new).

## Scope of This Implementation

### What Will Be Delivered

1. **Backend Abstraction Layer**: Complete protocol and base classes
2. **Three Working Backends**:
   - Qualtran backend (wrapping existing functionality)
   - PennyLane backend (new implementation)
   - Qiskit backend (new implementation)
3. **Configuration System**: Backend selection and configuration
4. **Comprehensive Tests**: Validation across all backends
5. **Documentation**: Architecture guide and usage examples

### Git Branch Strategy

This work will be organized into logical branches for evaluation:

```
main
 └─> bkk_backend (current)
      └─> bkk_backend_pennylane
           └─> bkk_backend_qiskit
```

**Branch Breakdown:**
- `bkk_backend`: Backend abstraction + Qualtran backend (validates existing functionality)
- `bkk_backend_pennylane`: Adds PennyLane backend (validates abstraction design)
- `bkk_backend_qiskit`: Adds Qiskit backend (demonstrates scalability)

Each branch represents a stable checkpoint where work can be evaluated.

## Current State Analysis

### Hard-Coded Dependencies

**analysis/unitary.py:**
- Direct imports: `qualtran`, `pyLIQTR`, `cirq`
- Functions hard-code framework: `encode_pauli_lcu_qualtran()`, `encode_pauli_lcu_pyliqtr()`
- Returns Qualtran `Bloq` objects directly

**analysis/algorithm.py:**
- Imports: `qualtran.bloqs.phase_estimation`, `pyLIQTR.qubitization`
- Functions: `build_qpe_qualtran_textbook()`, `build_qpe_pyliqtr_qubitized()`

**analysis/operators.py:**
- Uses numpy arrays and scipy for matrix operations
- Could be framework-agnostic with proper abstraction

### What Stays in Front End

- `Hamiltonian` class and file readers (`analysis/hamiltonian.py`)
- Configuration system (`analysis/config_types.py`)
- Analysis orchestration (`analysis/analysis.py`)
- Matrix operations and error analysis
- Pauli string utilities (`common/pauli_string.py`)

## Proposed Architecture

### Core Abstractions

#### 1. Backend Protocol

```python
class Backend(Protocol):
    """Protocol defining required backend operations."""
    
    def encode_pauli_trotter(...) -> Unitary: ...
    def encode_pauli_lcu(...) -> Unitary: ...
    def build_phase_estimation(...) -> Unitary: ...
    def estimate_resources(...) -> ResourceEstimate: ...
    def to_matrix(...) -> np.ndarray: ...
```

#### 2. Unitary Abstract Class

```python
class Unitary(ABC):
    """Framework-agnostic unitary operator."""
    
    @abstractmethod
    def controlled(self, num_controls=1) -> 'Unitary': ...
    
    @abstractmethod
    def adjoint(self) -> 'Unitary': ...
    
    @abstractmethod
    def to_matrix(self) -> np.ndarray: ...
    
    @abstractmethod
    def estimate_resources(self) -> ResourceEstimate: ...
```

#### 3. Backend-Specific Unitary Implementations

```python
class QualtranUnitary(Unitary):
    """Wraps Qualtran Bloq objects."""
    def __init__(self, bloq: Bloq): ...

class PennyLaneUnitary(Unitary):
    """Wraps PennyLane QNode objects."""
    def __init__(self, qnode: qml.QNode): ...

class QiskitUnitary(Unitary):
    """Wraps Qiskit QuantumCircuit objects."""
    def __init__(self, circuit: QuantumCircuit): ...
```

### Module Structure

```
analysis/
├── backend/
│   ├── __init__.py          # Backend registry and loading
│   ├── protocol.py          # Backend protocol definition
│   ├── base.py              # Base Unitary class
│   ├── qualtran_backend.py  # Qualtran implementation
│   ├── pennylane_backend.py # PennyLane implementation
│   ├── qiskit_backend.py    # Qiskit implementation
│   └── types.py             # Shared types (ResourceEstimate, etc.)
├── unitary.py               # Refactored to use backends
├── algorithm.py             # Refactored to use backends
└── ... (existing files)
```

## Implementation Plan

### Phase 1: Abstraction + Qualtran Backend (Branch: bkk_backend)

**Goal:** Establish abstraction without breaking existing functionality

**Tasks:**

1. **Create Backend Infrastructure** (`analysis/backend/`)
   - Define `Backend` protocol
   - Implement `Unitary` abstract class
   - Create `BackendRegistry` for discovery and loading
   - Define `ResourceEstimate` and other shared types

2. **Implement Qualtran Backend**
   - Extract existing Qualtran/pyLIQTR code from `unitary.py`
   - Wrap in `QualtranBackend` class
   - Implement `QualtranUnitary` wrapper around `Bloq`
   - Ensure all existing operations work

3. **Refactor Core Modules**
   - Update `unitary.py` to dispatch through backend
   - Update `algorithm.py` to use `Unitary` abstraction
   - Modify configuration to support backend selection

4. **Testing**
   - All existing tests must pass
   - Add backend-specific unit tests
   - Add integration tests for backend switching

**Success Criteria:**
- ✓ All existing tests pass
- ✓ Can select backend via configuration
- ✓ Zero functional changes to analysis results
- ✓ Code coverage ≥ previous level

### Phase 2: PennyLane Backend (Branch: bkk_backend_pennylane)

**Goal:** Validate abstraction with a second framework

**Tasks:**

1. **Implement PennyLane Backend** (`analysis/backend/pennylane_backend.py`)
   - Pauli string time evolution using `qml.PauliRot`
   - Trotterization using sequential PennyLane operations
   - Resource estimation via operation counting
   - Matrix representation via `qml.matrix()`

2. **PennyLane-Specific Features**
   - Leverage PennyLane's automatic differentiation (future)
   - Support PennyLane's built-in optimizers
   - Interface with PennyLane's hardware plugins

3. **Testing**
   - Compare PennyLane results with Qualtran
   - Validate resource estimates match
   - Test backend switching mid-workflow

4. **Documentation**
   - PennyLane backend capabilities and limitations
   - Installation and configuration guide
   - Example workflows

**Success Criteria:**
- ✓ PennyLane backend passes core test suite
- ✓ Resource estimates within 5% of Qualtran
- ✓ Can run complete analysis with PennyLane only
- ✓ Documentation includes PennyLane examples

### Phase 3: Qiskit Backend (Branch: bkk_backend_qiskit)

**Goal:** Demonstrate scalability to third framework

**Tasks:**

1. **Implement Qiskit Backend** (`analysis/backend/qiskit_backend.py`)
   - Pauli string evolution using Qiskit's `PauliEvolutionGate`
   - Trotterization using `TrotterQRTE`
   - Resource estimation via Qiskit's transpiler
   - Matrix representation via `Operator(circuit).data`

2. **Qiskit-Specific Features**
   - Leverage Qiskit's hardware optimization
   - Support Qiskit's noise models
   - Interface with IBM hardware backends (future)

3. **Testing**
   - Three-way comparison: Qualtran vs PennyLane vs Qiskit
   - Validate consistency across backends
   - Performance benchmarking

4. **Documentation**
   - Complete backend comparison guide
   - When to use which backend
   - Migration guide from old API

**Success Criteria:**
- ✓ Qiskit backend passes core test suite
- ✓ Three backends produce consistent results
- ✓ Documentation covers all backends
- ✓ Example demonstrating backend selection

## Design Decisions

### 1. Protocol vs ABC for Backend

**Decision:** Use `Protocol` for Backend, `ABC` for Unitary

**Rationale:**
- Backend uses Protocol: allows wrapping external classes without inheritance
- Unitary uses ABC: we control these classes, want enforcement

### 2. Lazy Backend Imports

**Decision:** Backends imported only when selected

```python
def get_backend(name: str) -> Backend:
    if name == "qualtran":
        from .qualtran_backend import QualtranBackend
        return QualtranBackend()
    elif name == "pennylane":
        from .pennylane_backend import PennyLaneBackend
        return PennyLaneBackend()
    # ...
```

**Rationale:**
- Users don't need all backends installed
- Reduces startup time
- Avoids dependency conflicts

### 3. Capability Discovery

**Decision:** Backends declare capabilities explicitly

```python
class Backend(Protocol):
    @property
    def capabilities(self) -> Set[str]: ...

# Usage
if "qubitized_phase_estimation" in backend.capabilities:
    result = backend.build_qubitized_qpe(...)
else:
    logger.warning(f"{backend.name} doesn't support qubitized QPE")
```

### 4. Resource Estimation Normalization

**Decision:** Unified `ResourceEstimate` dataclass

```python
@dataclass
class ResourceEstimate:
    num_qubits: int
    t_gates: int
    clifford_gates: int
    rotation_gates: int = 0
    measurements: int = 0
    depth: Optional[int] = None
    backend_specific: Dict[str, Any] = field(default_factory=dict)
```

### 5. Error Handling Strategy

**Decision:** Explicit exceptions for unsupported operations

```python
class UnsupportedOperationError(Exception):
    """Raised when backend doesn't support requested operation."""
    pass

# In backend implementation
def encode_double_factorization(self, ...):
    raise UnsupportedOperationError(
        f"{self.name} doesn't support double factorization encoding"
    )
```

## Backend Comparison

### Qualtran Backend

**Strengths:**
- Existing QHAT integration
- Comprehensive Bloq library
- Strong pyLIQTR integration
- Well-tested in QHAT context

**Limitations:**
- Pre-1.0 release (API instability)
- Limited active development visibility
- Installation can be challenging
- Pinned to specific versions by pyLIQTR

**Best For:**
- Block encoding methods
- LCU-based approaches
- Existing QHAT workflows

### PennyLane Backend

**Strengths:**
- Active development and community
- Excellent documentation
- Automatic differentiation
- Hardware integration (various providers)
- Quantum ML capabilities

**Limitations:**
- Less focus on fault-tolerant computing
- Resource estimation less mature
- May not support all QHAT encoding methods

**Best For:**
- Variational algorithms (future)
- Hardware execution
- Research and prototyping
- Gradient-based optimization

### Qiskit Backend

**Strengths:**
- Industry standard
- Extensive tooling
- IBM hardware access
- Mature resource estimation
- Large community

**Limitations:**
- Complex API surface
- Some features require Qiskit Terra/Aer splits
- Hardware-focused (not always ideal for simulation)

**Best For:**
- Industry collaboration
- IBM hardware targeting
- Standard compliance
- Resource estimation for real hardware

## Testing Strategy

### Unit Tests

Each backend has isolated unit tests:
- `test_qualtran_backend.py`
- `test_pennylane_backend.py`
- `test_qiskit_backend.py`

Test coverage:
- Pauli string time evolution
- Trotterization (all orders)
- Block encodings (where supported)
- Phase estimation
- Resource estimation
- Matrix representation

### Integration Tests

Cross-backend validation:
- `test_backend_consistency.py`: Compare results across backends
- `test_backend_switching.py`: Switch backends mid-workflow
- `test_configuration.py`: Backend selection from config files

### Regression Tests

Ensure new architecture doesn't break existing functionality:
- All existing test suite must pass with Qualtran backend
- Numerical results must match previous implementation
- Performance must be within 10% of previous implementation

### Performance Benchmarks

Compare backend performance:
- Resource estimation speed
- Matrix generation speed
- Memory usage
- Scaling with system size

## Migration Path for Users

### Configuration File Changes

**Old (implicit Qualtran):**
```toml
[unitary]
method = "ramped trotter"
use_library = "qualtran"
```

**New (explicit backend selection):**
```toml
[backend]
name = "qualtran"  # or "pennylane", "qiskit"

[unitary]
method = "ramped trotter"
```

### Python API (No Changes Required)

```python
# User code remains unchanged
config = Configuration("config.toml")
results = analyze(config)
```

Backend selection is purely configuration-driven.

## Risks and Mitigation

### Risk 1: Backends Produce Different Results

**Mitigation:**
- Comprehensive cross-backend validation tests
- Document expected deviations (e.g., different gate sets)
- Provide comparison tools for users

### Risk 2: Performance Regression

**Mitigation:**
- Benchmark suite comparing old vs new
- Profile hot paths
- Optimize backend dispatch if needed

### Risk 3: Installation Complexity

**Mitigation:**
- Make backends optional dependencies
- Clear installation documentation
- Conda environment files for each backend

### Risk 4: Maintenance Burden

**Mitigation:**
- Automated testing across all backends
- Clear backend responsibilities
- Community contributions encouraged

## Success Metrics

### Immediate (This Branch Series)

- [ ] All three backends implemented and tested
- [ ] 100% existing test pass rate
- [ ] <5% performance overhead
- [ ] Documentation complete for all backends

### Medium-Term (6 Months)

- [ ] At least one user adopts non-Qualtran backend
- [ ] Community contribution of a 4th backend
- [ ] Zero critical bugs in backend system

### Long-Term (1+ Years)

- [ ] Default backend can change without code rewrites
- [ ] Backend ecosystem enables rapid tool evaluation
- [ ] QHAT resilient to upstream framework changes

## Future Extensions

### Additional Backends (Priority Order)

1. **CUDA-Q**: NVIDIA hardware optimization, new Logical toolkit
2. **Q# / Azure**: Microsoft ecosystem, chemistry tools
3. **pytket**: Quantinuum hardware optimization
4. **Classiq**: High-level synthesis

### Enhanced Capabilities

- **Partial Backend Support**: Hybrid workflows using multiple backends
- **Backend Hints**: User preferences for speed vs accuracy
- **Caching Layer**: Reuse expensive computations across backends
- **Hardware Execution**: Direct submission to quantum hardware

## Deliverables Checklist

### Documentation

- [ ] `BACKEND_SEPARATION_PLAN.md` (this file)
- [ ] `BACKEND_SEPARATION_DETAILED.md` (detailed design)
- [ ] `docs/BACKEND_GUIDE.md` (user guide)
- [ ] `docs/BACKEND_DEVELOPER.md` (backend developer guide)
- [ ] Updated `README.md` with backend information

### Code

- [ ] `analysis/backend/` module complete
- [ ] `QualtranBackend` implementation
- [ ] `PennyLaneBackend` implementation
- [ ] `QiskitBackend` implementation
- [ ] Refactored `unitary.py` and `algorithm.py`
- [ ] Updated configuration system

### Tests

- [ ] Backend unit tests (3 files)
- [ ] Integration tests
- [ ] Regression tests
- [ ] Performance benchmarks

### Examples

- [ ] Example using Qualtran backend
- [ ] Example using PennyLane backend
- [ ] Example using Qiskit backend
- [ ] Example comparing backends

## Timeline Estimate

Given LLM-assisted development:

- **Phase 1 (Qualtran)**: ~4-6 hours
  - Abstraction design: 1 hour
  - Implementation: 2-3 hours
  - Testing: 1-2 hours

- **Phase 2 (PennyLane)**: ~3-4 hours
  - Backend implementation: 2 hours
  - Testing and validation: 1-2 hours

- **Phase 3 (Qiskit)**: ~3-4 hours
  - Backend implementation: 2 hours
  - Testing and comparison: 1-2 hours

- **Documentation**: ~2 hours throughout

**Total: ~12-16 hours of focused work**

## Conclusion

This plan provides a concrete roadmap for implementing backend separation in QHAT with three working backends. The phased approach via git branches allows evaluation at each milestone while building toward a comprehensive solution that future-proofs QHAT against ecosystem changes.

The implementation will demonstrate:
1. **Feasibility**: Backend abstraction works in practice
2. **Flexibility**: Multiple frameworks can coexist
3. **Maintainability**: Clear architecture and testing
4. **Value**: Real benefits from framework diversity

This is an ambitious but achievable goal that will significantly enhance QHAT's capabilities and longevity.
