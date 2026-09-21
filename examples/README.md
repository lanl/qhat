# QHAT Backend Examples

This directory contains examples demonstrating how to use QHAT with different quantum computing backends.

## Available Examples

### 1. Backend Comparison (`backend_comparison_example.py`)

Compares all three backends (Qualtran, PennyLane, Qiskit) side-by-side:
- Creates a simple 3-qubit Hamiltonian
- Encodes it using Trotterization with each backend
- Compares resource estimates across backends
- Validates numerical consistency

**Run:**
```bash
python3.11 backend_comparison_example.py
```

**What it demonstrates:**
- How to use `get_backend()` to load different backends
- Resource estimation differences between backends
- Matrix conversion and unitarity checks
- Cross-backend consistency validation

### 2. PennyLane Backend (`pennylane_backend_example.py`)

Focused example using PennyLane backend:
- First and second-order Trotterization
- Resource comparison between Trotter orders
- Controlled and adjoint operations
- Matrix-based verification

**Run:**
```bash
python3.11 pennylane_backend_example.py
```

**What it demonstrates:**
- PennyLane-specific features
- Comparing different Trotter formulas
- Testing controlled() and adjoint() operations
- Verifying mathematical properties (U† U = I)

## Requirements

Each example requires different dependencies:

**All examples:**
```bash
# QHAT with backend system (already in this branch)
```

**Qualtran backend:**
```bash
pip install qualtran pyLIQTR cirq
```

**PennyLane backend:**
```bash
pip install pennylane
```

**Qiskit backend:**
```bash
pip install qiskit qiskit-aer
```

Examples will skip backends that aren't installed.

## Usage Patterns

### Basic Backend Usage

```python
from qhat.analysis.backend import get_backend

# Load a backend
backend = get_backend("pennylane")

# Check capabilities
print(backend.capabilities)

# Encode Hamiltonian
unitary = backend.encode_pauli_trotter(
    pauli_strings=pauli_dict,
    trotter_order="second order",
    evolution_time=1.0,
    num_steps=10,
    num_qubits=3
)

# Get resources
resources = unitary.estimate_resources()
print(f"T gates: {resources.t_gates}")
```

### Using Backends with QHAT Configuration

```python
from qhat.analysis.config_types import BackendConfiguration, UnitaryConfiguration
from qhat.analysis.unitary import encode_as_unitary
from qhat.analysis.hamiltonian import get_physical_hamiltonian

# Configure backend
backend_config = BackendConfiguration()
backend_config.set_backend("qiskit")

# Configure unitary encoding
unitary_config = UnitaryConfiguration()
unitary_config.encode_ramped_trotter(
    timestep=1.0,
    energy_error=0.001,
    trotter_order="second order"
)
unitary_config.backend_name = backend_config.name
unitary_config.backend_options = backend_config.options

# Load Hamiltonian
ham_config = HamiltonianConfiguration()
ham_config.load_pauli_strings("hamiltonian.json")
hamiltonian = get_physical_hamiltonian(ham_config)

# Encode with specified backend
from qhat.analysis.backend import get_backend
backend = get_backend(backend_config.name, **backend_config.options)
unitary = encode_as_unitary(unitary_config, hamiltonian, 1.0, backend=backend)
```

### Backend Selection via TOML Configuration

```toml
[backend]
name = "pennylane"

[backend.options]
device = "default.qubit"

[unitary]
method = "ramped trotter"
trotter_order = "second order"
energy_error = 0.001
```

## Tips for Writing Your Own Examples

1. **Always check if backend is available:**
   ```python
   try:
       backend = get_backend("backend_name")
   except ImportError:
       print("Backend not available")
       return
   ```

2. **Use small systems for matrix validation:**
   - Matrices grow as 2^N × 2^N
   - Keep N ≤ 10 for matrix operations
   - Use resource estimation for larger systems

3. **Compare backends on same problem:**
   - Use identical Hamiltonians and parameters
   - Check both resources and numerical results
   - Document any differences

4. **Handle backend-specific features gracefully:**
   ```python
   if "qubitized_qpe" in backend.capabilities:
       result = backend.build_qubitized_qpe(...)
   else:
       print(f"{backend.name} doesn't support qubitized QPE")
   ```

## Expected Output

Examples will show:
- ✓ Successfully loaded backend
- Resource estimates (qubits, gates, depth)
- Matrix properties (norm, unitarity)
- Comparisons between methods/backends
- ✗ Skip messages for unavailable backends

## Troubleshooting

**"Backend not found":**
- Check the backend name ("qualtran", "pennylane", or "qiskit")
- Ensure dependencies are installed

**"UnsupportedOperationError":**
- The backend doesn't support that operation
- Try a different backend or different method
- Check `backend.capabilities` for supported operations

**Import errors:**
- Install the required backend package
- Check Python version compatibility
- Verify virtual environment is activated

**Numerical differences:**
- Small differences (< 1e-6) are normal due to floating point
- Large differences may indicate bugs or different algorithms
- Check Trotter parameters (order, steps) are identical

## Next Steps

After running these examples:
1. Try with your own Hamiltonians
2. Compare performance on larger systems
3. Explore backend-specific optimizations
4. Integrate backends into your QHAT workflows

## Contributing

To add a new example:
1. Create `<topic>_example.py`
2. Follow the pattern of existing examples
3. Add entry to this README
4. Include error handling for missing dependencies
5. Add clear documentation and print statements
