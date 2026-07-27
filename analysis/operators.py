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


# =================================================================================================
# Conversion Utilities
# =================================================================================================

def convert_unitary_eigenvalues_to_eigenenergies(unitary_eigenvalues, timestep, hbar=1.0):
    """
    Convert unitary eigenvalues e^(-iφ) to Hamiltonian eigenenergies E.

    The time evolution operator is U = exp(-iHt/ℏ), so if H has eigenenergy E,
    then U has eigenvalue exp(-iEt/ℏ) = exp(-iφ) where φ = Et/ℏ is the eigenphase.

    Assumption: The Hamiltonian has been shifted and scaled (by existing code) such that
    all eigenenergies produce eigenphases in the range [0, 2π), preventing aliasing.

    Parameters:
        unitary_eigenvalues: Complex eigenvalues of unitary U = exp(-iHt/ℏ) (on unit circle)
        timestep: Time evolution parameter t (units: ℏ/energy, e.g., ℏ/Hartree)
        hbar: Value of ℏ (default: 1.0 for atomic units)

    Returns:
        tuple: (eigenenergies, eigenphases)
            eigenenergies: Real eigenvalues of Hamiltonian (shifted/scaled, units: energy)
            eigenphases: Phases φ ∈ [0, 2π) where φ = Et/ℏ
    """
    # Extract phases using np.angle, which returns [-π, π]
    # Since U = exp(-iφ), we have arg(U) = -φ (mod 2π)
    phases_neg_pi_to_pi = np.angle(unitary_eigenvalues)

    # Eigenphase φ = -arg(U), then map to [0, 2π) convention
    eigenphases_raw = -phases_neg_pi_to_pi
    eigenphases = np.where(eigenphases_raw < 0,
                           eigenphases_raw + 2*np.pi,
                           eigenphases_raw)

    # Convert to eigenenergies: E = φ * ℏ / t
    # These are shifted/scaled eigenenergies that correspond to the phases in [0, 2π)
    eigenenergies = eigenphases * hbar / timestep

    return eigenenergies, eigenphases


# =================================================================================================
# OperatorRepresentation Class
# =================================================================================================

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

        # Get matrix dimension for logging
        if representation == 'dense_matrix':
            dim = data.shape[0]
        else:
            dim = len(data['eigenvalues'])

        logger.info(
            f"OperatorRepresentation created: type={operator_type}, "
            f"shifted={energy_shifted}, repr={representation}, dim={dim}, "
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
            logger.debug(
                f"OperatorRepresentation.get() cache hit: "
                f"type={operator_type}, shifted={energy_shifted}, repr={representation}"
            )
            return self._cache[cache_key]

        logger.info(
            f"OperatorRepresentation.get() cache miss - converting: "
            f"from (type={self._original_type}, shifted={self._original_shifted}, repr={self._original_repr}) "
            f"to (type={operator_type}, shifted={energy_shifted}, repr={representation})"
        )

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
        needs_conversion = (
            self._original_type != operator_type or
            self._original_shifted != energy_shifted
        )
        if needs_conversion:
            logger.debug(
                f"Converting eigenvalues: "
                f"from (type={self._original_type}, shifted={self._original_shifted}) "
                f"to (type={operator_type}, shifted={energy_shifted})"
            )
        eigenvalues = self._convert_eigenvalues(
            eigenvalues,
            from_type=self._original_type.lower(),
            from_shifted=self._original_shifted,
            to_type=operator_type.lower(),
            to_shifted=energy_shifted
        )

        # Step 3: Convert to desired representation
        if representation == 'eigendecomposition':
            result = {
                'eigenvalues': eigenvalues,
                'eigenvectors': eigenvectors
            }
        else:  # dense_matrix
            # Reconstruct matrix: M = V @ diag(λ) @ V^(-1)
            # Note: For Hermitian matrices (from eigh), V is unitary so V† = V^(-1)
            # For general matrices (from eig), V is not necessarily unitary, so we must use V^(-1)
            logger.debug(
                f"Reconstructing dense matrix from eigendecomposition: "
                f"type={operator_type}, shifted={energy_shifted}, dim={len(eigenvalues)}"
            )
            if operator_type == 'hamiltonian':
                # Eigenvectors from eigh are orthonormal, use V† (faster)
                result = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T
            else:
                # Eigenvectors from eig are not necessarily orthonormal, use V^(-1)
                result = eigenvectors @ np.diag(eigenvalues) @ np.linalg.inv(eigenvectors)

        # Cache and return
        self._cache[cache_key] = result
        logger.info(
            f"Conversion complete and cached: "
            f"type={operator_type}, shifted={energy_shifted}, repr={representation}"
        )
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
            logger.debug(
                f"Eigendecomposition cache hit: type={operator_type}, shifted={energy_shifted}"
            )
            return self._cache[cache_key]

        if representation == 'eigendecomposition':
            # Already have it
            logger.debug(
                f"Eigendecomposition already available: type={operator_type}, shifted={energy_shifted}"
            )
            return data

        # Need to compute eigendecomposition from dense matrix
        dim = data.shape[0]
        if operator_type == 'hamiltonian':
            # Hamiltonians (shifted or unshifted) are Hermitian, use eigh
            # IMPORTANT: Must use eigh (not eig) for Hermitian matrices to get
            # properly orthonormal eigenvectors. Using eig() causes reconstruction errors.
            logger.info(
                f"Computing eigendecomposition (eigh): type={operator_type}, "
                f"shifted={energy_shifted}, dim={dim} [O(N³) operation]"
            )
            eigenvalues, eigenvectors = scipy.linalg.eigh(data)
        else:
            # Time-evolution operators - use general eig
            logger.info(
                f"Computing eigendecomposition (eig): type={operator_type}, "
                f"shifted={energy_shifted}, dim={dim} [O(N³) operation]"
            )
            eigenvalues, eigenvectors = scipy.linalg.eig(data)

        result = {
            'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors
        }

        self._cache[cache_key] = result
        logger.info(
            f"Eigendecomposition computed and cached: type={operator_type}, shifted={energy_shifted}"
        )
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
        1. Apply energy shift if not present: λ_unshifted → λ_shifted
        2. Convert operator type: H ↔ U or U ↔ H
        3. Remove energy shift if not desired: λ_shifted → λ_unshifted

        The reason for this is that the energy shift is chosen so that the shifted energies are
        non-negative.  When converting from U to H with the shift applied, we know that the phases
        should be in the range [0, 2π).  That allows us to correctly handle aliasing due to the
        fact that θ + 2π k is indistinguishable from θ when taking the logarithm.

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

        # Step 1: Apply energy shift if not present (go to shifted state)
        if not from_shifted:
            result = self._apply_energy_shift(result, from_type)

        # Step 2: Convert operator type (both shifted at this point)
        if from_type != to_type:
            if from_type == 'hamiltonian' and to_type == 'time_evolution':
                result = self._hamiltonian_to_time_evolution(result)
            elif from_type == 'time_evolution' and to_type == 'hamiltonian':
                result = self._time_evolution_to_hamiltonian(result)

        # Step 3: Remove energy shift if not desired
        if not to_shifted:
            result = self._remove_energy_shift(result, to_type)

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

        # get phase in [0, 2π) instead of (-π, π]
        phase = -1 * np.angle(eigenvalues)
        phase = np.where(phase < 0, phase + 2*np.pi, phase)

        # convert phase angle to energy eigenvalue
        H_eigen = phase * self.hbar / self.timestep

        # ensure Hamiltonian eigenvalues are real
        assert(all(abs(H_eigen.imag) < 1.0e-9 * abs(H_eigen)))
        return H_eigen.real

    def _apply_energy_shift(self, eigenvalues: np.ndarray, operator_type: str) -> np.ndarray:
        """
        Apply energy shift to eigenvalues.
        """
        if operator_type == 'hamiltonian':
            return eigenvalues - self.energy_shift
        else:  # time_evolution
            phase_factor = np.exp(1j * self.energy_shift * self.timestep / self.hbar)
            print(f"phase factor = {phase_factor}")
            print(f"energy shift = {self.energy_shift}")
            print(f"time step    = {self.timestep}")
            print(f"hbar         = {self.hbar}")
            return phase_factor * eigenvalues

    def _remove_energy_shift(self, eigenvalues: np.ndarray, operator_type: str) -> np.ndarray:
        """
        Remove energy shift from eigenvalues (inverse of _apply_energy_shift).
        """
        if operator_type == 'hamiltonian':
            return eigenvalues + self.energy_shift
        else:  # time_evolution
            phase_factor = np.exp(-1j * self.energy_shift * self.timestep / self.hbar)
            print(f"phase factor = {phase_factor}")
            print(f"energy shift = {self.energy_shift}")
            print(f"time step    = {self.timestep}")
            print(f"hbar         = {self.hbar}")
            return phase_factor * eigenvalues
