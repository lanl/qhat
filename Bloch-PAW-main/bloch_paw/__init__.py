"""Bloch--UPAW quantum resource estimation pipeline.

Extracts PAW integrals from GPAW, computes the Hamiltonian one-norm λ
(Eq. 37 of the paper), and estimates Toffoli gate and qubit counts
(Eqs. 98--99) for fault-tolerant quantum phase estimation.

Typical workflow::

    # 1. Extract (requires GPAW)
    extractor = PawExtractor(calc)
    extractor.export_hdf5("integrals.h5")

    # 2. Read + compute one-norm
    reader = PawReader("integrals.h5")
    calc = OneNormCalculator.from_reader(reader)
    lam = calc.lambda_one_norm()

    # 3. Resource estimate
    est = ResourceEstimator.from_hdf5("integrals.h5")
    toffolis = est.toffoli_count_per_be(...)
    qubits = est.total_qubits(...)
"""

from .reader import PawReader
from .one_norm import OneNormCalculator
from .resources import ResourceEstimator

__all__ = [
    "PawReader",
    "OneNormCalculator",
    "ResourceEstimator",
]
