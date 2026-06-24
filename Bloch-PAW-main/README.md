# bloch-paw

Quantum resource estimation using Projector Augmented Wave (PAW) methods. Extracts Hamiltonian ingredients from GPAW DFT calculations and estimates quantum computing costs (qubits and Toffoli gates) for quantum phase estimation.

## Pipeline

```
1. GPAW: run DFT calculation
2. PawExtractor: extract PAW terms, export to HDF5
3. PawReader + OneNormCalculator: load HDF5, compute λ (one-norm)
4. ResourceEstimator: estimate qubit and gate counts
```

## Install

Requires Python 3.12+. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

For development/testing:

```bash
uv sync --extra test
```

## Usage

```python
from ase.build import bulk
from gpaw import GPAW

from bloch_paw.extractor import PawExtractor
from bloch_paw.reader import PawReader
from bloch_paw.one_norm import OneNormCalculator
from bloch_paw.resources import ResourceEstimator

# === Phase 1 (materials): run DFT, extract Hamiltonian, save to disk ===

# 1. Run a GPAW DFT calculation
lattice_constant = 3.67  # lattice constant in Å
k_mesh_size = (3, 3, 3)  # Monkhorst-Pack k-point grid dimensions
h_cut = 0.30             # real-space grid spacing (Å); sets planewave cutoff via ecut ∝ 1/h²
nbands = 5               # number of Kohn-Sham bands to compute
xc = 'PBE'               # exchange-correlation functional

atoms = bulk('H', 'fcc', a=lattice_constant)
calc = GPAW(mode="lcao", basis='dzp', h=h_cut,
            kpts={'size': k_mesh_size, 'gamma': True},
            xc=xc, nbands=nbands,
            symmetry={'point_group': False, 'time_reversal': False})
atoms.calc = calc
calc.get_potential_energy(atoms)

# 2. Extract PAW ingredients to HDF5
thr_rho = 1e-3    # threshold for zeroing small pseudo-density elements
thr_D = 1e-2      # threshold for zeroing small projector density matrix elements
thr_C = 1e-2      # threshold for zeroing small on-site Coulomb tensor elements
thr_h = 1e-5      # threshold for zeroing small one-body matrix elements
thr_kappa = 1e-5  # threshold for zeroing small two-body integral elements

extractor = PawExtractor(calc, nbands=nbands,
                         thr_rho=thr_rho, thr_D=thr_D, thr_C=thr_C,
                         thr_h=thr_h, thr_kappa=thr_kappa)
extractor.export_hdf5("paw_ingredients.h5")

# === Hamiltonian saved to paw_ingredients.h5 — hand off to QRE team ===

# === Phase 2 (QRE): load Hamiltonian, compute resource estimates ===

# 3. Load and compute one-norm
thr_rank = 3e-5    # eigenvalue threshold for rank truncation in SVD
sv_floor = 1e-12   # floor for zeroing small singular values
scale_floor = 1e-12  # floor for Frobenius-norm rescaling

reader = PawReader("paw_ingredients.h5")
with reader:
    inputs = reader.to_calculator_inputs(lazy=True)
    one_norm_calc = OneNormCalculator(**inputs, thr_rank=thr_rank,
                                     sv_floor=sv_floor, scale_floor=scale_floor)
    lam = one_norm_calc.lambda_one_norm()
    avg_rank_l, avg_rank_zero = one_norm_calc.compute_average_rank()

# 4. Estimate resources
eps_chem = 1.6e-3    # chemical accuracy in Hartree (1 kcal/mol)
eps_qpe = eps_chem / 5  # QPE precision: 1/5 of error budget
# Error budget split: QPE, LCU truncation, coefficient loading, finite basis, finite k-mesh

est = ResourceEstimator.from_hdf5("paw_ingredients.h5")
toffoli = est.toffoli_count_per_be(Rl=avg_rank_l, R0=avg_rank_zero)
qubits = est.total_qubits(Rl=avg_rank_l, R0=avg_rank_zero,
                           lam=lam, eps_qpe=eps_qpe)
```

See `examples/` for Jupyter notebooks walking through each pipeline stage.

## Tests

```bash
uv run pytest tests/ -v
uv run pytest tests/ -v -m "not integration"  # skip tests requiring GPAW
```

## Project structure

```
bloch_paw/          Core library
  extractor.py        Extract PAW ingredients from GPAW
  reader.py           Read HDF5 ingredient files
  one_norm.py         Compute λ one-norm for quantum phase estimation
  resources.py        Estimate qubit and gate counts

examples/           Jupyter notebooks and scripts demonstrating the pipeline
tests/              Unit and integration tests
```
