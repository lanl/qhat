#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import json
from typing import Optional, Tuple, Dict

import numpy as np
import h5py


class PawReader:
    """Read PAW integrals from an HDF5 file produced by PawExtractor.export_hdf5().

    Deserialises the k-mesh, lattice vectors, smooth pseudo pair-density ρ̃,
    C^a and D^a tensors, and optional one-body / two-body integrals. The
    loaded data can be passed directly to OneNormCalculator via
    ``to_calculator_inputs()``.

    Supports both eager loading (all data into memory) and lazy loading
    (rho returned as an h5py dataset handle for block-by-block access).

    Args:
        filepath: Path to the HDF5 file.
        kmesh_region: Which mesh to load --- 'bz' (full Brillouin zone)
            or 'ibz' (irreducible wedge).
        kmesh_space: Coordinate system --- 'cart' (Å⁻¹) or 'reduced'.

    Attributes (populated after load):
        N_atoms: Number of atoms (int or None if absent from HDF5).
        h_pq: One-body matrix elements (Nk, Nb, Nk, Nb) or None.
        one_body_unit: Unit string ('eV' or 'Ha') or None.
        kappa_pqrs: Full rank-8 two-body integrals or None.
        kappa_unit: Unit string or None.
    """

    def __init__(self,
                 filepath: str,
                 *,
                 kmesh_region: str = "bz",
                 kmesh_space: str = "cart"):
        self.filepath = filepath
        if kmesh_region not in {"bz", "ibz"}:
            raise ValueError("kmesh_region must be 'bz' or 'ibz'")
        if kmesh_space not in {"cart", "reduced"}:
            raise ValueError("kmesh_space must be 'cart' or 'reduced'")
        self.kmesh_region = kmesh_region
        self.kmesh_space = kmesh_space

        # HDF5 file handle for lazy reads (managed by context manager)
        self._h5: Optional[h5py.File] = None

        # Metadata (populated by .load())
        self.N_atoms: Optional[int] = None

        # One-body (populated by .load() / .load_one_body())
        self.h_pq: Optional[np.ndarray] = None
        self.one_body_unit: Optional[str] = None
        self.one_body_k_diagonal: Optional[bool] = None

        # Two-body (populated by .load() / .load_two_body())
        self.kappa_pqrs: Optional[np.ndarray] = None
        self.kappa_unit: Optional[str] = None
        self.kappa_g0: Optional[str] = None
        self.kappa_split_written: Optional[bool] = None
        self.kappa_soft_rank8: Optional[np.ndarray] = None
        self.kappa_paw_rank8: Optional[np.ndarray] = None

        # ---- reduced two-body block-diagonal tensors (populated by .load() / .load_two_body()) ----
        self.kappa_A: Optional[np.ndarray] = None
        self.kappa_B: Optional[np.ndarray] = None
        self.kappa_A_soft: Optional[np.ndarray] = None
        self.kappa_A_paw: Optional[np.ndarray] = None
        self.kappa_B_soft: Optional[np.ndarray] = None
        self.kappa_B_paw: Optional[np.ndarray] = None

    # ----------------------------- context manager -----------------------------

    def _ensure_open(self):
        """Open the HDF5 file if not already open."""
        if self._h5 is None:
            self._h5 = h5py.File(self.filepath, 'r')

    def close(self):
        """Close the HDF5 file if open."""
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None

    def __enter__(self):
        self._ensure_open()
        return self

    def __exit__(self, *args):
        self.close()

    # ------------------------------- public API -------------------------------

    def load(self):
        """
        Returns
        -------
        k_mesh : (Nk, 3) float64
        L_size : (3,) int64
        a_vectors : (3, 3) float64
        N_pw : int
        rho : complex128 ndarray, shape (Nk, Nb, Nk, Nb, Nx, Ny, Nz)
        Ca_dict : dict[int, float64 ndarray]
        Da_dict : dict[int, complex128 ndarray]
        meta : dict
        """
        with h5py.File(self.filepath, "r") as h5:
            # --- meta (optional but handy) ---
            meta = {}
            if "meta_json" in h5.attrs:
                try:
                    meta = json.loads(h5.attrs["meta_json"])
                except Exception:
                    meta = {}

            # --- core datasets ---
            N_pw = int(np.asarray(h5["Npw"][()]))

            L_size = np.asarray(h5["supercell_size"][()], dtype=np.int64)

            a_vectors = np.asarray(h5["lattice"]["A_direct_ang"][()], dtype=np.float64)

            k_mesh = self._read_kmesh(h5)

            rho = np.asarray(h5["rho_tilde"]["data"][()], dtype=np.complex128)

            Ca_dict = self._read_atom_dict(h5, group_name="C_tensor", dtype=np.float64)

            Da_dict = self._read_atom_dict(h5, group_name="D_tensor", dtype=np.complex128)

            # --- optional metadata ---
            if "N_atoms" in h5:
                self.N_atoms = int(np.asarray(h5["N_atoms"][()]))
            else:
                # Backward compat: infer from Ca_dict keys
                self.N_atoms = len(Ca_dict) if Ca_dict else None

            # --- optional blocks ---
            self._maybe_read_one_body(h5, meta)
            self._maybe_read_two_body(h5, meta)

        return k_mesh, L_size, a_vectors, N_pw, rho, Ca_dict, Da_dict, meta

    def _load_lazy(self):
        """Load all datasets except rho, which is returned as an h5py dataset handle.

        The HDF5 file must remain open for the lifetime of the returned handle.
        Use the context manager (``with reader:``) to manage file lifetime.
        """
        self._ensure_open()
        h5 = self._h5

        meta = {}
        if "meta_json" in h5.attrs:
            try:
                meta = json.loads(h5.attrs["meta_json"])
            except Exception:
                meta = {}

        N_pw = int(np.asarray(h5["Npw"][()]))
        L_size = np.asarray(h5["supercell_size"][()], dtype=np.int64)
        a_vectors = np.asarray(h5["lattice"]["A_direct_ang"][()], dtype=np.float64)
        k_mesh = self._read_kmesh(h5)
        rho = h5["rho_tilde"]["data"]  # h5py dataset handle (not loaded)
        Ca_dict = self._read_atom_dict(h5, group_name="C_tensor", dtype=np.float64)
        Da_dict = self._read_atom_dict(h5, group_name="D_tensor", dtype=np.complex128)

        if "N_atoms" in h5:
            self.N_atoms = int(np.asarray(h5["N_atoms"][()]))
        else:
            self.N_atoms = len(Ca_dict) if Ca_dict else None

        self._maybe_read_one_body(h5, meta)
        self._maybe_read_two_body(h5, meta)

        return k_mesh, L_size, a_vectors, N_pw, rho, Ca_dict, Da_dict, meta

    def to_calculator_inputs(self,
                             include_one_body: bool = True,
                             include_two_body: bool = False,
                             lazy: bool = False) -> Dict[str, object]:
        """Return a dict of keyword arguments for ``OneNormCalculator(**inputs)``.

        Args:
            include_one_body: Include h_pq in the output if present in the file.
            include_two_body: Include kappa_pqrs in the output if present.
            lazy: If True, rho is returned as an h5py dataset handle (the HDF5
                file must remain open via the context manager).

        Returns:
            Dictionary with keys: k_mesh, L_size, a_vectors, N_pw, rho,
            Ca_dict, Da_dict, and optionally h_pq / kappa_pqrs.
        """
        if lazy:
            k_mesh, L_size, a_vectors, N_pw, rho, Ca, Da, _ = self._load_lazy()
        else:
            k_mesh, L_size, a_vectors, N_pw, rho, Ca, Da, _ = self.load()
        out = {
            "k_mesh": k_mesh,
            "L_size": tuple(int(x) for x in L_size),
            "a_vectors": a_vectors,
            "N_pw": int(N_pw),
            "rho": rho,
            "Ca_dict": Ca,
            "Da_dict": Da,
        }
        if include_one_body and (self.h_pq is not None):
            out["h_pq"] = self.h_pq

        # NOTE: The two-body payload changed. We keep backward-compat (rank-8)
        # but prefer reduced tensors if available.
        if include_two_body:
            if (self.kappa_A is not None) and (self.kappa_B is not None):
                out["kappa_A"] = self.kappa_A
                out["kappa_B"] = self.kappa_B
            elif self.kappa_pqrs is not None:
                out["kappa_pqrs"] = self.kappa_pqrs
        return out

    def load_one_body(self) -> Tuple[Optional[np.ndarray], Optional[str], Optional[bool]]:
        """
        Read only the one-body integral block if present.
        Returns (h_pq, unit, k_diagonal).
        """
        with h5py.File(self.filepath, "r") as h5:
            meta = {}
            if "meta_json" in h5.attrs:
                try:
                    meta = json.loads(h5.attrs["meta_json"])
                except Exception:
                    pass
            self._maybe_read_one_body(h5, meta)
        return self.h_pq, self.one_body_unit, self.one_body_k_diagonal

    def load_two_body(self) -> Tuple[
        Optional[np.ndarray], Optional[str], Optional[str], Optional[bool],
        Optional[np.ndarray], Optional[np.ndarray]
    ]:
        """
        Read only the two-body integral block if present.

        Returns
        -------
        kappa_pqrs : (Nk,Nb,Nk,Nb,Nk,Nb,Nk,Nb) complex128 or None
        unit        : str or None
        g0          : str or None
        split_written : bool or None
        kappa_soft_rank8 : ndarray or None
        kappa_paw_rank8  : ndarray or None
        """
        with h5py.File(self.filepath, "r") as h5:
            meta = {}
            if "meta_json" in h5.attrs:
                try:
                    meta = json.loads(h5.attrs["meta_json"])
                except Exception:
                    pass
            self._maybe_read_two_body(h5, meta)
        return (self.kappa_pqrs, self.kappa_unit, self.kappa_g0,
                self.kappa_split_written, self.kappa_soft_rank8, self.kappa_paw_rank8)

    # ------------------------------ helpers ------------------------------

    def _read_kmesh(self, h5: h5py.File) -> np.ndarray:
        g = h5["kmesh"]
        if self.kmesh_space == "cart":
            ds = g[self.kmesh_region]["cart_1_perA"]
        else:
            ds = g[self.kmesh_region]["reduced"]
        k = np.asarray(ds[()], dtype=np.float64)
        if k.ndim != 2 or k.shape[1] != 3:
            raise ValueError(f"k-mesh must be (Nk,3); got {k.shape}")
        return k

    @staticmethod
    def _read_atom_dict(h5: h5py.File, *, group_name: str, dtype) -> Dict[int, np.ndarray]:
        if group_name not in h5:
            return {}
        out: Dict[int, np.ndarray] = {}
        for name, ds in h5[group_name].items():
            if not isinstance(ds, h5py.Dataset):
                continue
            # Expect names like 'atom_0000'
            try:
                a = int(name.split("_")[-1])
            except Exception as e:
                raise ValueError(f"Unexpected dataset name in '{group_name}': {name}") from e
            out[a] = np.asarray(ds[()], dtype=dtype)
        return out

    # --- read optional one-body block and mirror to attributes/meta ---
    def _maybe_read_one_body(self, h5: h5py.File, meta: dict) -> None:
        if "one_body" not in h5 or "H_kikj" not in h5["one_body"]:
            self.h_pq = None
            self.one_body_unit = None
            self.one_body_k_diagonal = None
            return
        ds = h5["one_body"]["H_kikj"]
        self.h_pq = np.asarray(ds[()], dtype=np.complex128)
        unit = ds.attrs.get("unit", None)
        if unit is None and isinstance(meta.get("one_body"), dict):
            unit = meta["one_body"].get("unit", None)
        self.one_body_unit = str(unit) if unit is not None else None
        kdiag = ds.attrs.get("k_diagonal", None)
        if kdiag is None and isinstance(meta.get("one_body"), dict):
            kdiag = meta["one_body"].get("k_diagonal", None)
        self.one_body_k_diagonal = bool(kdiag) if kdiag is not None else None

    # --- read optional two-body block (rank-8) and mirror to attributes/meta ---
    def _maybe_read_two_body(self, h5: h5py.File, meta: dict) -> None:
        # defaults if absent
        self.kappa_pqrs = None
        self.kappa_soft_rank8 = None
        self.kappa_paw_rank8 = None
        self.kappa_unit = None
        self.kappa_g0 = None
        self.kappa_split_written = None

        # defaults for reduced blocks if absent
        self.kappa_A = None
        self.kappa_B = None
        self.kappa_A_soft = None
        self.kappa_A_paw = None
        self.kappa_B_soft = None
        self.kappa_B_paw = None

        if "two_body" not in h5:
            return
        g2 = h5["two_body"]

        # ---- NEW FORMAT (preferred): reduced block-diagonal tensors ----
        # These datasets are written by the updated PawExtractor:
        #   kappa_A_kikplkjkpl, kappa_B_kikjkplkpl
        if ("kappa_A_kikplkjkpl" in g2) and ("kappa_B_kikjkplkpl" in g2):
            dsA = g2["kappa_A_kikplkjkpl"]
            dsB = g2["kappa_B_kikjkplkpl"]

            self.kappa_A = np.asarray(dsA[()], dtype=np.complex128)
            self.kappa_B = np.asarray(dsB[()], dtype=np.complex128)

            # attrs / meta (prefer dataset attrs, fallback to meta)
            unit = dsA.attrs.get("unit", None) or dsB.attrs.get("unit", None)
            g0 = dsA.attrs.get("g0", None) or dsB.attrs.get("g0", None)
            split = None
            if isinstance(meta.get("two_body"), dict):
                unit = unit or meta["two_body"].get("unit", None)
                g0 = g0 or meta["two_body"].get("g0", None)
                split = meta["two_body"].get("split_written", None)

            self.kappa_unit = str(unit) if unit is not None else None
            self.kappa_g0 = str(g0) if g0 is not None else None
            self.kappa_split_written = bool(split) if split is not None else None

            # optional parts (reduced splits)
            if "kappa_A_soft" in g2:
                self.kappa_A_soft = np.asarray(g2["kappa_A_soft"][()], dtype=np.complex128)
            if "kappa_A_paw" in g2:
                self.kappa_A_paw = np.asarray(g2["kappa_A_paw"][()], dtype=np.complex128)
            if "kappa_B_soft" in g2:
                self.kappa_B_soft = np.asarray(g2["kappa_B_soft"][()], dtype=np.complex128)
            if "kappa_B_paw" in g2:
                self.kappa_B_paw = np.asarray(g2["kappa_B_paw"][()], dtype=np.complex128)

            # Backward compat attributes remain None in the reduced format:
            #   self.kappa_pqrs, self.kappa_soft_rank8, self.kappa_paw_rank8
            return

        # ---- OLD FORMAT (backward compat): rank-8 tensor ----
        if "kappa" not in g2:
            return

        dsK = g2["kappa"]
        self.kappa_pqrs = np.asarray(dsK[()], dtype=np.complex128)

        # attrs / meta
        unit = dsK.attrs.get("unit", None)
        g0 = dsK.attrs.get("g0", None)
        split = None
        if isinstance(meta.get("two_body"), dict):
            unit = unit or meta["two_body"].get("unit", None)
            g0 = g0 or meta["two_body"].get("g0", None)
            split = meta["two_body"].get("split_written", None)

        self.kappa_unit = str(unit) if unit is not None else None
        self.kappa_g0 = str(g0) if g0 is not None else None
        self.kappa_split_written = bool(split) if split is not None else None

        # optional parts
        if "kappa_soft_rank8" in g2:
            self.kappa_soft_rank8 = np.asarray(g2["kappa_soft_rank8"][()], dtype=np.complex128)
        if "kappa_paw_rank8" in g2:
            self.kappa_paw_rank8 = np.asarray(g2["kappa_paw_rank8"][()], dtype=np.complex128)
