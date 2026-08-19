# Command-Line Parameters for Configuration Files

This feature allows you to pass parameters from the command line into your configuration files, enabling easy parametric studies and configuration reuse without modifying the config file itself.

## Usage

Pass parameters using the `-p` or `--param` flag with `KEY=VALUE` format:

```bash
python3.11 driver.py my_config.py -p distance=1.5 -p basis=sto-6g
```

Multiple parameters can be specified:

```bash
python3.11 driver.py my_config.py -p distance=1.5 -p basis=sto-6g -p num_steps=100
```

## Accessing Parameters in Configuration Files

Parameters are available in two ways within your configuration script:

### Method 1: Direct access (if parameter is guaranteed to exist)

```python
# Access parameter directly as a variable
hamiltonian.add_atom("H", distance, 0, 0)
hamiltonian.basis = basis
```

**Note:** This will raise a `NameError` if the parameter wasn't provided on the command line.

### Method 2: Using the `params` dictionary with defaults (recommended)

```python
# Access with a default value
distance = params.get('distance', 1.30)  # Defaults to 1.30 if not provided
basis_set = params.get('basis', 'sto-3g')  # Defaults to 'sto-3g' if not provided

hamiltonian.add_atom("H", distance, 0, 0)
hamiltonian.basis = basis_set
```

## Type Conversion

Values are automatically evaluated as Python expressions when possible:

| Command Line | Type | Value |
|--------------|------|-------|
| `-p x=1.5` | float | `1.5` |
| `-p n=42` | int | `42` |
| `-p items=[1,2,3]` | list | `[1, 2, 3]` |
| `-p flag=True` | bool | `True` |
| `-p name=hello` | string | `'hello'` (eval fails, treated as string) |

## Examples

### Example 1: Parametric distance study

**Configuration file** (`distance_study_config.py`):
```python
# Use command-line parameter for distance
distance = params.get('distance', 1.30)

general.print_verbose()
general.file_stub = f"Be-H_{distance:.2f}_sto-6g"
general.file_format = "default"

hamiltonian.add_atom("Be", 0.00, 0, 0)
hamiltonian.add_atom("H", distance, 0, 0)
hamiltonian.basis = "sto-6g"
hamiltonian.num_active_occupied = 3
hamiltonian.num_active_vacant = 3
```

**Run with different distances**:
```bash
python3.11 driver.py distance_study_config.py -p distance=1.2
python3.11 driver.py distance_study_config.py -p distance=1.5
python3.11 driver.py distance_study_config.py -p distance=2.0
```

### Example 2: Basis set comparison

```bash
# Same config, different basis sets
python3.11 driver.py config.py -p basis=sto-3g
python3.11 driver.py config.py -p basis=sto-6g
python3.11 driver.py config.py -p basis=6-31g
```

### Example 3: Multiple parameters

**Configuration file**:
```python
distance = params.get('distance', 1.30)
basis = params.get('basis', 'sto-3g')
n_occ = int(params.get('num_occupied', 3))
n_vac = int(params.get('num_vacant', 3))

general.file_stub = f"Be-H_{distance}_{basis}_occ{n_occ}_vac{n_vac}"

hamiltonian.add_atom("Be", 0.00, 0, 0)
hamiltonian.add_atom("H", distance, 0, 0)
hamiltonian.basis = basis
hamiltonian.num_active_occupied = n_occ
hamiltonian.num_active_vacant = n_vac
```

**Run**:
```bash
python3.11 driver.py config.py -p distance=1.5 -p basis=sto-6g -p num_occupied=4 -p num_vacant=4
```

### Example 4: Scripting parametric studies

```bash
#!/bin/bash
# Run a sweep over different distances

for dist in 1.0 1.2 1.4 1.6 1.8 2.0; do
    echo "Running with distance = $dist"
    python3.11 driver.py my_config.py -p distance=$dist
done
```

## Advanced Usage

### Type conversion and validation

```python
# Ensure integer types
num_steps = int(params.get('num_steps', 100))

# Convert string to float
energy = float(params.get('energy', '0.0'))

# Parse complex structures
if 'atoms' in params:
    for atom_type, position in params['atoms']:
        hamiltonian.add_atom(atom_type, position, 0, 0)
```

### Conditional logic based on parameters

```python
mode = params.get('mode', 'production')

if mode == 'debug':
    general.print_debug()
elif mode == 'verbose':
    general.print_verbose()
else:
    general.print_default()

# Adjust other settings based on mode
algorithm.num_phase_qubits = 8 if mode == 'debug' else 12
```

## Error Handling

### Invalid format
```bash
# Missing '=' will raise ValueError
python3.11 driver.py config.py -p bad_param
# Error: Parameter must be in KEY=VALUE format, got: bad_param
```

### Using undefined parameters without defaults
```python
# In config file:
x = params['missing_key']  # KeyError if not provided
# Better:
x = params.get('missing_key', default_value)  # Safe with default
```

## See Also

- `examples/cmdline_params_demo_config.py` - Full demonstration of parameter usage
- `tests/test_cmdline_params.py` - Test suite showing various use cases
