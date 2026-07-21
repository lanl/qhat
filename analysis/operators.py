"""
Unified representation of quantum operators with lazy conversion between forms.

This module provides the OperatorRepresentation class for managing quantum operators
that can be represented as either Hamiltonians or time-evolution operators, with or
without energy shifts, in various computational representations (dense matrices, eigendecompositions).

Phase 2 scope: Dense matrices and eigendecompositions only.
Phase 3 scope: Matrix-free operators (LinearOperator).
"""

import numpy as np
import scipy.linalg
import logging
from typing import Dict, Optional, Union, Literal

logger = logging.getLogger(__name__)


class OperatorRepresentation:
    """
    Unified representation of quantum operators with lazy conversion between forms.

    Represents a quantum operator that can be viewed as:
    - Hamiltonian H or time-evolution operator U = exp(-i*H*t)
    - Energy-shifted or unshifted
    - Dense matrix or eigendecomposition

    Internal storage: Stores original data and caches converted forms.
    Converts lazily on demand and caches results for efficiency.

    Phase 2 scope: Dense matrices and eigendecompositions only.
    Phase 3 scope: Matrix-free operators (LinearOperator).

    Examples
    --------
    >>> # Create from exact Hamiltonian
    >>> H_exact = np.array([[1.0, 0.5], [0.5, 2.0]])
    >>> exact_op = OperatorRepresentation(
    ...     data=H_exact,
    ...     operator_type='hamiltonian',
    ...     energy_shifted=False,
    ...     representation='dense_matrix',
    ...     timestep=1.0
    ... )
    >>>
    >>> # Get as unshifted time-evolution operator
    >>> U_exact = exact_op.get(
    ...     operator_type='time_evolution',
    ...     energy_shifted=False,
    ...     representation='dense_matrix'
    ... )
    >>>
    >>> # Get eigendecomposition
    >>> eigendata = exact_op.get(
    ...     operator_type='hamiltonian',
    ...     representation='eigendecomposition'
    ... )
    >>> eigenvalues = eigendata['eigenvalues']
    >>> eigenvectors = eigendata['eigenvectors']
    """

    def __init__(
        self,
        data: Union[np.ndarray, Dict],
        operator_type: Literal['hamiltonian', 'time_evolution'],
        energy_shifted: bool = False,
        representation: Literal['dense_matrix', 'eigendecomposition'] = 'dense_matrix',
        timestep: Optional[float] = None,
        energy_shift: float = 0.0,
        hbar: float = 1.0
    ):
        """
        Initialize from a known representation.

        Parameters
        ----------
        data : np.ndarray or dict
            The operator data:
            - If representation='dense_matrix': 2D numpy array
            - If representation='eigendecomposition': dict with keys 'eigenvalues', 'eigenvectors'
        operator_type : {'hamiltonian', 'time_evolution'}
            Whether this data represents a Hamiltonian H or time-evolution operator U
        energy_shifted : bool, optional
            Whether this representation includes an energy shift (default: False)
        representation : {'dense_matrix', 'eigendecomposition'}, optional
            The form of the input data (default: 'dense_matrix')
        timestep : float, optional
            Time evolution parameter t, required for converting between H and U
        energy_shift : float, optional
            Energy shift value E (default: 0.0, meaning no shift)
        hbar : float, optional
            Value of ℏ in natural units (default: 1.0)

        Raises
        ------
        ValueError
            If data format doesn't match the specified representation
        """
        # Store original specification
        self._original_data = data
        self._original_type = operator_type
        self._original_shifted = energy_shifted
        self._original_repr = representation

        # Store parameters
        self.timestep = timestep
        self.energy_shift = energy_shift
        self.hbar = hbar

        # Validate data format
        if representation == 'dense_matrix':
            if not isinstance(data, np.ndarray):
                raise ValueError(
                    f"representation='dense_matrix' requires data to be np.ndarray, "
                    f"got {type(data)}"
                )
            if data.ndim != 2 or data.shape[0] != data.shape[1]:
                raise ValueError(
                    f"Dense matrix must be 2D square array, got shape {data.shape}"
                )
        elif representation == 'eigendecomposition':
            if not isinstance(data, dict):
                raise ValueError(
                    f"representation='eigendecomposition' requires data to be dict, "
                    f"got {type(data)}"
                )
            if 'eigenvalues' not in data or 'eigenvectors' not in data:
                raise ValueError(
                    "Eigendecomposition dict must contain 'eigenvalues' and 'eigenvectors' keys"
                )

        # Initialize cache: maps (operator_type, energy_shifted, representation) -> data
        self._cache: Dict[tuple, Union[np.ndarray, Dict]] = {}

        # Store original in cache
        cache_key = (operator_type, energy_shifted, representation)
        self._cache[cache_key] = data

        logger.debug(
            f"Created OperatorRepresentation: type={operator_type}, "
            f"shifted={energy_shifted}, repr={representation}, "
            f"timestep={timestep}, E_shift={energy_shift}"
        )

    def get(
        self,
        operator_type: Optional[Literal['hamiltonian', 'time_evolution']] = None,
        energy_shifted: Optional[bool] = None,
        representation: Optional[Literal['dense_matrix', 'eigendecomposition']] = None
    ) -> Union[np.ndarray, Dict]:
        """
        Get operator in requested form, converting lazily as needed.

        Parameters
        ----------
        operator_type : {'hamiltonian', 'time_evolution'}, optional
            Desired operator type. If None, keeps current type.
        energy_shifted : bool, optional
            Whether to include energy shift. If None, keeps current shift state.
        representation : {'dense_matrix', 'eigendecomposition'}, optional
            Desired representation. If None, keeps current representation.

        Returns
        -------
        np.ndarray or dict
            The operator in the requested form:
            - If representation='dense_matrix': 2D numpy array
            - If representation='eigendecomposition': dict with 'eigenvalues', 'eigenvectors'

        Raises
        ------
        ValueError
            If conversion requires timestep but it was not provided

        Notes
        -----
        Results are cached, so repeated calls with the same parameters are fast.
        """
        # Use original values if not specified
        if operator_type is None:
            operator_type = self._original_type
        if energy_shifted is None:
            energy_shifted = self._original_shifted
        if representation is None:
            representation = self._original_repr

        # Check cache
        cache_key = (operator_type, energy_shifted, representation)
        if cache_key in self._cache:
            logger.debug(f"Cache hit: {cache_key}")
            return self._cache[cache_key]

        logger.debug(f"Cache miss: {cache_key}, computing conversion")

        # Need to convert - strategy: convert through eigendecomposition (universal intermediate)
        # Path: current form -> eigendecomp -> target form

        # Step 1: Get eigendecomposition in current operator_type and energy_shifted state
        # Start from original, convert operator_type and energy_shift at eigenvalue level

        # Get eigendecomposition (may need to compute it)
        eigendata = self._get_eigendecomposition(
            self._original_type,
            self._original_shifted,
            self._original_repr,
            self._original_data
        )

        eigenvalues = eigendata['eigenvalues'].copy()  # Copy so we don't modify cache
        eigenvectors = eigendata['eigenvectors']

        # Step 2: Convert eigenvalues to desired operator_type and energy_shifted state
        eigenvalues = self._convert_eigenvalues(
            eigenvalues,
            from_type=self._original_type,
            from_shifted=self._original_shifted,
            to_type=operator_type,
            to_shifted=energy_shifted
        )

        # Step 3: Convert to desired representation
        if representation == 'eigendecomposition':
            result = {
                'eigenvalues': eigenvalues,
                'eigenvectors': eigenvectors
            }
        else:  # dense_matrix
            # Reconstruct matrix: M = V @ diag(λ) @ V†
            result = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T

        # Cache and return
        self._cache[cache_key] = result
        return result

    def _get_eigendecomposition(
        self,
        operator_type: str,
        energy_shifted: bool,
        representation: str,
        data: Union[np.ndarray, Dict]
    ) -> Dict:
        """
        Get eigendecomposition for given operator state.

        Returns dict with 'eigenvalues' and 'eigenvectors' keys.
        """
        cache_key = (operator_type, energy_shifted, 'eigendecomposition')
        if cache_key in self._cache:
            return self._cache[cache_key]

        if representation == 'eigendecomposition':
            # Already have it
            return data

        # Need to compute eigendecomposition from dense matrix
        if operator_type == 'hamiltonian' and not energy_shifted:
            # Unshifted Hamiltonian - should be Hermitian, use eigh for efficiency
            eigenvalues, eigenvectors = scipy.linalg.eigh(data)
        else:
            # Everything else - use general eig
            eigenvalues, eigenvectors = scipy.linalg.eig(data)

        result = {
            'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors
        }

        self._cache[cache_key] = result
        return result

    def _convert_eigenvalues(
        self,
        eigenvalues: np.ndarray,
        from_type: str,
        from_shifted: bool,
        to_type: str,
        to_shifted: bool
    ) -> np.ndarray:
        """
        Convert eigenvalues between different operator representations.

        Conversion path:
        1. Remove energy shift if present: λ_shifted → λ_unshifted
        2. Convert operator type: H ↔ U
        3. Apply energy shift if desired: λ_unshifted → λ_shifted

        Parameters
        ----------
        eigenvalues : np.ndarray
            Input eigenvalues
        from_type : str
            Source operator type ('hamiltonian' or 'time_evolution')
        from_shifted : bool
            Whether source eigenvalues include energy shift
        to_type : str
            Target operator type
        to_shifted : bool
            Whether target eigenvalues should include energy shift

        Returns
        -------
        np.ndarray
            Converted eigenvalues
        """
        result = eigenvalues.copy()

        # Step 1: Remove energy shift if present (go to unshifted state)
        if from_shifted:
            result = self._remove_energy_shift(result, from_type)

        # Step 2: Convert operator type (both unshifted at this point)
        if from_type != to_type:
            if from_type == 'hamiltonian' and to_type == 'time_evolution':
                result = self._hamiltonian_to_time_evolution(result)
            elif from_type == 'time_evolution' and to_type == 'hamiltonian':
                result = self._time_evolution_to_hamiltonian(result)

        # Step 3: Apply energy shift if desired
        if to_shifted:
            result = self._apply_energy_shift(result, to_type)

        return result

    def _hamiltonian_to_time_evolution(self, eigenvalues: np.ndarray) -> np.ndarray:
        """
        Convert Hamiltonian eigenvalues to time-evolution eigenvalues.

        U = exp(-i*H*t/ℏ), so λ_U = exp(-i*λ_H*t/ℏ)
        """
        if self.timestep is None:
            raise ValueError(
                "Conversion from Hamiltonian to time-evolution operator requires timestep"
            )
        return np.exp(-1j * eigenvalues * self.timestep / self.hbar)

    def _time_evolution_to_hamiltonian(self, eigenvalues: np.ndarray) -> np.ndarray:
        """
        Convert time-evolution eigenvalues to Hamiltonian eigenvalues.

        H = i*ℏ*log(U)/t, so λ_H = i*ℏ*log(λ_U)/t

        Note: This involves logarithm branch cuts. We choose the principal branch
        such that phases are in (-π, π].
        """
        if self.timestep is None:
            raise ValueError(
                "Conversion from time-evolution to Hamiltonian operator requires timestep"
            )

        # log(exp(-i*E*t/ℏ)) should give -i*E*t/ℏ, so E = i*ℏ*log(...)/t
        # np.log uses principal branch with phase in (-π, π]
        return 1j * self.hbar * np.log(eigenvalues) / self.timestep

    def _apply_energy_shift(self, eigenvalues: np.ndarray, operator_type: str) -> np.ndarray:
        """
        Apply energy shift to eigenvalues.

        CONVENTION CHANGE (2026-07-21): Now matches hamiltonian.py driver convention.
        For Hamiltonians: H_shifted = H + E*I, so λ_shifted = λ + E
        For time-evolution: U_shifted = exp(-i*E*t/ℏ)*U, so λ_U_shifted = exp(-i*E*t/ℏ)*λ_U

        This matches the driver.py convention where energy_shift(+E) adds E to eigenvalues.
        """
        if operator_type == 'hamiltonian':
            return eigenvalues + self.energy_shift  # FIXED: was - (sign error)
        else:  # time_evolution
            phase_factor = np.exp(-1j * self.energy_shift * self.timestep / self.hbar)  # FIXED: was +1j
            return phase_factor * eigenvalues

    def _remove_energy_shift(self, eigenvalues: np.ndarray, operator_type: str) -> np.ndarray:
        """
        Remove energy shift from eigenvalues (inverse of _apply_energy_shift).

        CONVENTION CHANGE (2026-07-21): Now matches hamiltonian.py driver convention.
        For Hamiltonians: H = H_shifted - E*I, so λ = λ_shifted - E
        For time-evolution: U = exp(+i*E*t/ℏ)*U_shifted, so λ_U = exp(+i*E*t/ℏ)*λ_U_shifted

        This matches the driver.py convention where energy_shift(+E) adds E to eigenvalues.
        """
        if operator_type == 'hamiltonian':
            return eigenvalues - self.energy_shift  # FIXED: was + (sign error)
        else:  # time_evolution
            phase_factor = np.exp(1j * self.energy_shift * self.timestep / self.hbar)  # FIXED: was -1j
            return phase_factor * eigenvalues
