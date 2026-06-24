"""Integration tests requiring a live GPAW installation.

Run with::

    uv run pytest -m integration

Skipped automatically in normal ``pytest tests/`` runs.

Gold-standard numbers (LCBO mode, metallic hydrogen FCC, 3×3×3 k-mesh,
5 bands, GPAW 25.7.0):

    One-norm lambda  :  319.51
    Toffoli / query  :  206244
    Total qubits     :   36805

All three are asserted within 1 % to tolerate minor GPAW version drift.
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.integration
def test_lcbo_end_to_end(tmp_path):
    """Full LCBO pipeline: DFT → ingredients → one-norm → resource estimates.

    Metallic hydrogen FCC, 3×3×3 k-mesh, 5 bands, LCAO/dzp, PBE.
    """
    from ase.build import bulk
    from gpaw import GPAW

    from bloch_paw.extractor import PawExtractor
    from bloch_paw.one_norm import OneNormCalculator
    from bloch_paw.reader import PawReader
    from bloch_paw.resources import ResourceEstimator

    # GPAW enables np.seterr(divide='raise') when it detects pytest (its
    # intentional debug mode).  The full pipeline contains known harmless
    # divide-by-zeros in GPAW's PAW-setup and in our ingredients code;
    # restore them to warnings for this test.
    hdf5_path = str(tmp_path / "paw_ingredients.h5")
    with np.errstate(divide="warn", over="warn", invalid="warn"):
        # ------------------------------------------------------------------ #
        # 1. DFT calculation
        # ------------------------------------------------------------------ #
        atoms = bulk("H", "fcc", a=3.67)
        calc = GPAW(
            mode="lcao",
            basis="dzp",
            h=0.30,
            kpts={"size": (3, 3, 3), "gamma": True},
            xc="PBE",
            nbands=5,
            symmetry={"point_group": False, "time_reversal": False},
        )
        atoms.calc = calc
        calc.get_potential_energy(atoms)

        # ------------------------------------------------------------------ #
        # 2. Extract PAW ingredients and export to HDF5
        # ------------------------------------------------------------------ #
        extractor = PawExtractor(
            calc,
            thr_rho=1.0e-3,
            thr_D=1.0e-2,
            thr_C=1.0e-2,
            thr_h=1e-5,
            thr_kappa=1e-5,
            nbands=5,
        )
        extractor.export_hdf5(filepath=hdf5_path)

        # ------------------------------------------------------------------ #
        # 3. One-norm calculation
        # ------------------------------------------------------------------ #
        reader = PawReader(hdf5_path)
        one_norm_calc = OneNormCalculator.from_reader(
            reader, thr_rank=3e-5, sv_floor=1e-12, scale_floor=1e-12
        )
        lam = one_norm_calc.lambda_one_norm()
        avg_rank_l, avg_rank_zero = one_norm_calc.compute_average_rank()

        # ------------------------------------------------------------------ #
        # 4. Resource estimation
        # ------------------------------------------------------------------ #
        est = ResourceEstimator.from_hdf5(hdf5_path)
        params = dict(br=20, Rl=avg_rank_l, R0=avg_rank_zero,
                      N1=10, N2=10, B=20)
        t_per_query = est.toffoli_count_per_be(**params)
        eps_chem = 1.6e-3    # chemical accuracy (Hartree)
        eps_qpe = eps_chem / 5  # QPE share of error budget
        q_total = est.total_qubits(**params, lam=lam, eps_qpe=eps_qpe)

    # ------------------------------------------------------------------ #
    # 5. Gold-standard assertions (1 % tolerance)
    # ------------------------------------------------------------------ #
    assert lam == pytest.approx(319.51, rel=0.01), (
        f"One-norm lambda {lam:.4f} deviates >1% from gold standard 319.51"
    )
    assert t_per_query == pytest.approx(206244, rel=0.01), (
        f"Toffoli/query {t_per_query} deviates >1% from gold standard 206244"
    )
    assert q_total == pytest.approx(36805, rel=0.01), (
        f"Total qubits {q_total} deviates >1% from gold standard 36805"
    )
