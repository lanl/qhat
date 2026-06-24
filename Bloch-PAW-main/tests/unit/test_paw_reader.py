from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import numpy as np
import h5py
import pytest

from bloch_paw.reader import PawReader


# ---------------------------------------------------------------------------
# Helpers to build a mock HDF5 file
# ---------------------------------------------------------------------------

class _MockDataset:
    """Lightweight stand-in for h5py.Dataset that passes isinstance checks."""

    def __init__(self, data, attrs=None):
        self._data = data
        self.attrs = attrs or {}

    def __getitem__(self, key):
        return self._data


def _make_mock_dataset(data, attrs=None):
    """Return a _MockDataset that behaves like an h5py.Dataset."""
    return _MockDataset(data, attrs=attrs)


def _make_mock_group(datasets: dict, subgroups: dict | None = None):
    """Return a MagicMock that behaves like an h5py.Group."""
    g = MagicMock()
    children = dict(datasets)
    if subgroups:
        children.update(subgroups)

    def getitem(key):
        return children[key]

    def contains(key):
        return key in children

    g.__getitem__ = MagicMock(side_effect=getitem)
    g.__contains__ = MagicMock(side_effect=contains)

    all_ds = {k: v for k, v in datasets.items()}
    g.items = MagicMock(return_value=list(all_ds.items()))
    return g


class _FakeH5File:
    """Minimal dict-like stand-in for h5py.File used as a context manager."""

    def __init__(self, tree: dict, attrs: dict | None = None):
        self._tree = tree
        self.attrs = attrs or {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __getitem__(self, key):
        return self._tree[key]

    def __contains__(self, key):
        return key in self._tree


def _build_fake_h5(
    Nk: int = 2,
    Nb: int = 2,
    grid: tuple[int, int, int] = (3, 3, 3),
    na: int = 2,
    n_atoms: int = 1,
    include_one_body: bool = False,
    include_two_body: bool = False,
    meta: dict | None = None,
) -> _FakeH5File:
    """Construct a _FakeH5File with the structure that PawReader expects."""
    Nx, Ny, Nz = grid
    rng = np.random.default_rng(0)

    Npw_ds = _make_mock_dataset(np.int64(42))
    supercell_ds = _make_mock_dataset(np.array([1, 1, 1], dtype=np.int64))

    a_direct = np.eye(3, dtype=np.float64) * 3.0
    lattice_grp = _make_mock_group({"A_direct_ang": _make_mock_dataset(a_direct)})

    k_cart = rng.random((Nk, 3)).astype(np.float64)
    k_reduced = rng.random((Nk, 3)).astype(np.float64)
    bz_grp = _make_mock_group({
        "cart_1_perA": _make_mock_dataset(k_cart),
        "reduced": _make_mock_dataset(k_reduced),
    })
    ibz_grp = _make_mock_group({
        "cart_1_perA": _make_mock_dataset(k_cart[:1]),
        "reduced": _make_mock_dataset(k_reduced[:1]),
    })
    kmesh_grp = _make_mock_group({}, subgroups={"bz": bz_grp, "ibz": ibz_grp})

    rho_data = rng.random((Nk, Nb, Nk, Nb, Nx, Ny, Nz)) + 0j
    rho_grp = _make_mock_group({"data": _make_mock_dataset(rho_data)})

    c_datasets = {}
    d_datasets = {}
    for i in range(n_atoms):
        name = f"atom_{i:04d}"
        c_datasets[name] = _make_mock_dataset(rng.random((na, na, na, na)))
        d_datasets[name] = _make_mock_dataset(
            rng.random((Nk, Nb, Nk, Nb, na, na)) + 0j
        )
    c_grp = _make_mock_group(c_datasets)
    d_grp = _make_mock_group(d_datasets)

    tree: dict = {
        "Npw": Npw_ds,
        "supercell_size": supercell_ds,
        "lattice": lattice_grp,
        "kmesh": kmesh_grp,
        "rho_tilde": rho_grp,
        "C_tensor": c_grp,
        "D_tensor": d_grp,
    }

    if include_one_body:
        h_kikj_ds = _make_mock_dataset(
            rng.random((Nk, Nb, Nk, Nb)) + 0j,
            attrs={"unit": "Ha", "k_diagonal": True},
        )
        ob_grp = _make_mock_group({"H_kikj": h_kikj_ds})
        tree["one_body"] = ob_grp

    if include_two_body:
        kappa_ds = _make_mock_dataset(
            rng.random((Nk, Nb, Nk, Nb, Nk, Nb, Nk, Nb)) + 0j,
            attrs={"unit": "eV", "g0": "wigner-seitz"},
        )
        tb_grp = _make_mock_group({"kappa": kappa_ds})
        tree["two_body"] = tb_grp

    attrs = {}
    if meta is not None:
        attrs["meta_json"] = json.dumps(meta)
    return _FakeH5File(tree, attrs=attrs)


# ---------------------------------------------------------------------------
# Construction tests
# ---------------------------------------------------------------------------

class TestPawReaderConstruction:
    """Test PawReader constructor validation."""

    def test_valid_construction(self):
        """Test that valid parameters are stored, including non-default region/space.

        PawReader supports bz/ibz regions and cart/reduced coordinate spaces.
        Defaults must be bz + cart.
        """
        reader = PawReader("/fake/path.h5")
        assert reader.filepath == "/fake/path.h5"
        assert reader.kmesh_region == "bz"
        assert reader.kmesh_space == "cart"

        reader2 = PawReader("/f.h5", kmesh_region="ibz", kmesh_space="reduced")
        assert reader2.kmesh_region == "ibz"
        assert reader2.kmesh_space == "reduced"

    def test_rejects_invalid_params(self):
        """Test that invalid kmesh_region or kmesh_space raises ValueError.

        Catches typos like 'full' or 'spherical' before they cause cryptic
        KeyError inside the HDF5 file.
        """
        with pytest.raises(ValueError, match="kmesh_region"):
            PawReader("/f.h5", kmesh_region="full")
        with pytest.raises(ValueError, match="kmesh_space"):
            PawReader("/f.h5", kmesh_space="spherical")


# ---------------------------------------------------------------------------
# load() tests
# ---------------------------------------------------------------------------

class TestPawReaderLoad:
    """Test the main load() method with mock HDF5 data."""

    @patch("bloch_paw.reader.h5py")
    def test_load_core_fields(self, mock_h5py):
        """Test that load() returns an 8-tuple with correct shapes and types.

        Validates k_mesh, L_size, a_vectors, N_pw, and rho — the five
        mandatory fields for OneNormCalculator construction.
        """
        Nk, Nb, grid = 2, 2, (3, 3, 3)
        fake = _build_fake_h5(Nk=Nk, Nb=Nb, grid=grid)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        result = reader.load()
        assert len(result) == 8

        k_mesh, L_size, a_vectors, N_pw, rho, Ca_dict, Da_dict, meta = result
        assert k_mesh.shape == (Nk, 3)
        assert k_mesh.dtype == np.float64
        np.testing.assert_array_equal(L_size, [1, 1, 1])
        assert a_vectors.shape == (3, 3)
        np.testing.assert_allclose(a_vectors, np.eye(3) * 3.0)
        assert N_pw == 42
        assert rho.shape == (Nk, Nb, Nk, Nb, *grid)
        assert rho.dtype == np.complex128

    @patch("bloch_paw.reader.h5py")
    def test_load_kmesh_shape(self, mock_h5py):
        """Test that k_mesh has shape (Nk, 3) for a non-trivial Nk.

        Uses Nk=4 to verify the reader doesn't hardcode the mesh size.
        """
        Nk = 4
        fake = _build_fake_h5(Nk=Nk)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        k_mesh, *_ = reader.load()
        assert k_mesh.shape == (Nk, 3)


# ---------------------------------------------------------------------------
# Atom dict parsing
# ---------------------------------------------------------------------------

class TestAtomDictParsing:
    """Test parsing of per-atom C^a and D^a tensors."""

    @patch("bloch_paw.reader.h5py")
    def test_Ca_Da_dicts(self, mock_h5py):
        """Test that atom dicts have integer keys and correct dtypes.

        Ca must be float64 (symmetric real tensor), Da must be complex128
        (Bloch orbital projections). Keys must be 0-indexed integers, not
        the raw 'atom_0000' strings.
        """
        n_atoms = 3
        fake = _build_fake_h5(n_atoms=n_atoms)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        *_, Ca_dict, Da_dict, _ = reader.load()
        assert set(Ca_dict.keys()) == {0, 1, 2}
        assert set(Da_dict.keys()) == {0, 1, 2}
        for arr in Ca_dict.values():
            assert arr.dtype == np.float64
        for arr in Da_dict.values():
            assert arr.dtype == np.complex128


# ---------------------------------------------------------------------------
# Meta / optional blocks
# ---------------------------------------------------------------------------

class TestOptionalBlocks:
    """Test parsing of optional HDF5 groups (metadata, one-body, two-body)."""

    @patch("bloch_paw.reader.h5py")
    def test_metadata_handling(self, mock_h5py):
        """Test that JSON metadata is parsed from file attrs, or {} if absent.

        Metadata carries provenance information (Nk, xc functional, etc.)
        but is not required for the calculation.
        """
        meta_in = {"description": "test run", "Nk": 2}
        fake = _build_fake_h5(meta=meta_in)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset
        reader = PawReader("/f.h5")
        *_, meta_out = reader.load()
        assert meta_out == meta_in

        # No meta -> empty dict
        fake2 = _build_fake_h5(meta=None)
        mock_h5py.File.return_value = fake2
        reader2 = PawReader("/f.h5")
        *_, meta_out2 = reader2.load()
        assert meta_out2 == {}

    @patch("bloch_paw.reader.h5py")
    def test_one_body_present_and_absent(self, mock_h5py):
        """Test that one-body data is loaded when present and None when absent.

        When present, h_pq must be complex128 and metadata (unit, k_diagonal)
        must be extracted from dataset attributes.
        """
        fake = _build_fake_h5(include_one_body=True)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset
        reader = PawReader("/f.h5")
        reader.load()
        assert reader.h_pq is not None
        assert reader.h_pq.dtype == np.complex128
        assert reader.one_body_unit == "Ha"
        assert reader.one_body_k_diagonal is True

        fake2 = _build_fake_h5(include_one_body=False)
        mock_h5py.File.return_value = fake2
        reader2 = PawReader("/f.h5")
        reader2.load()
        assert reader2.h_pq is None

    @patch("bloch_paw.reader.h5py")
    def test_two_body_present_and_absent(self, mock_h5py):
        """Test that two-body kappa is loaded when present and None when absent.

        When present, kappa_pqrs must be complex128 with unit and g0 metadata.
        """
        fake = _build_fake_h5(include_two_body=True)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset
        reader = PawReader("/f.h5")
        reader.load()
        assert reader.kappa_pqrs is not None
        assert reader.kappa_pqrs.dtype == np.complex128
        assert reader.kappa_unit == "eV"
        assert reader.kappa_g0 == "wigner-seitz"

        fake2 = _build_fake_h5(include_two_body=False)
        mock_h5py.File.return_value = fake2
        reader2 = PawReader("/f.h5")
        reader2.load()
        assert reader2.kappa_pqrs is None


# ---------------------------------------------------------------------------
# to_calculator_inputs()
# ---------------------------------------------------------------------------

class TestToCalculatorInputs:
    """Test the to_calculator_inputs() convenience method."""

    @patch("bloch_paw.reader.h5py")
    def test_base_keys_always_present(self, mock_h5py):
        """Test that the 7 mandatory keys for OneNormCalculator are always returned.

        to_calculator_inputs() is the bridge between reader and calculator —
        missing keys cause TypeError on construction.
        """
        fake = _build_fake_h5()
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        out = reader.to_calculator_inputs()
        expected_keys = {"k_mesh", "L_size", "a_vectors", "N_pw", "rho", "Ca_dict", "Da_dict"}
        assert expected_keys.issubset(set(out.keys()))

    @patch("bloch_paw.reader.h5py")
    def test_optional_fields_inclusion(self, mock_h5py):
        """Test that h_pq and kappa_pqrs are included/excluded based on flags.

        Users may want to exclude the expensive two-body term for fast
        approximate calculations.
        """
        fake = _build_fake_h5(include_one_body=True, include_two_body=True)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        out_with = reader.to_calculator_inputs(include_one_body=True, include_two_body=True)
        assert "h_pq" in out_with
        assert "kappa_pqrs" in out_with

        mock_h5py.File.return_value = _build_fake_h5(include_one_body=True)
        reader2 = PawReader("/f.h5")
        out_without = reader2.to_calculator_inputs(include_one_body=False)
        assert "h_pq" not in out_without

    @patch("bloch_paw.reader.h5py")
    def test_L_size_is_tuple_of_ints(self, mock_h5py):
        """Test that L_size is a tuple of Python ints, not numpy int64.

        OneNormCalculator unpacks L_size as (Lx, Ly, Lz) — numpy int64
        can cause issues with downstream integer arithmetic in some contexts.
        """
        fake = _build_fake_h5()
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        out = reader.to_calculator_inputs()
        assert isinstance(out["L_size"], tuple)
        assert all(isinstance(x, int) for x in out["L_size"])


# ---------------------------------------------------------------------------
# load_one_body / load_two_body standalone
# ---------------------------------------------------------------------------

class TestStandaloneLoadMethods:
    """Test standalone load methods that skip the full load() pipeline."""

    @patch("bloch_paw.reader.h5py")
    def test_load_one_body_present(self, mock_h5py):
        """Test that load_one_body returns (h_pq, unit, k_diagonal) when present."""
        fake = _build_fake_h5(include_one_body=True)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        h_pq, unit, k_diag = reader.load_one_body()
        assert h_pq is not None
        assert unit == "Ha"
        assert k_diag is True

    @patch("bloch_paw.reader.h5py")
    def test_load_one_body_absent(self, mock_h5py):
        """Test that load_one_body returns (None, None, None) when group is missing."""
        fake = _build_fake_h5(include_one_body=False)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        h_pq, unit, k_diag = reader.load_one_body()
        assert h_pq is None
        assert unit is None
        assert k_diag is None

    @patch("bloch_paw.reader.h5py")
    def test_load_two_body_returns_6_tuple(self, mock_h5py):
        """Test that load_two_body returns a 6-element tuple."""
        fake = _build_fake_h5(include_two_body=True)
        mock_h5py.File.return_value = fake
        mock_h5py.Dataset = _MockDataset

        reader = PawReader("/f.h5")
        result = reader.load_two_body()
        assert len(result) == 6


# ---------------------------------------------------------------------------
# Lazy read + context manager (uses real HDF5 files via tmp_path)
# ---------------------------------------------------------------------------

def _write_test_hdf5(filepath, *, Nk=2, Nb=2, grid=(3, 3, 3), na=2, seed=42):
    """Write a minimal HDF5 file with the structure PawReader expects."""
    Nx, Ny, Nz = grid
    rng = np.random.default_rng(seed)
    rho_data = rng.standard_normal((Nk, Nb, Nk, Nb, Nx, Ny, Nz)) + \
               1j * rng.standard_normal((Nk, Nb, Nk, Nb, Nx, Ny, Nz))

    a_lat = 3.0
    a_vectors = np.eye(3) * a_lat
    L_size = np.array([Nk, 1, 1], dtype=np.int64)
    A = np.diag([Nk * a_lat, a_lat, a_lat])
    vol = np.linalg.det(A)
    B = np.vstack([
        2 * np.pi * np.cross(A[1], A[2]) / vol,
        2 * np.pi * np.cross(A[2], A[0]) / vol,
        2 * np.pi * np.cross(A[0], A[1]) / vol,
    ])
    c_mesh = np.zeros((Nk, 3))
    for i in range(Nk):
        c_mesh[i, 0] = i / Nk
    k_cart = c_mesh @ B

    with h5py.File(filepath, "w") as h5:
        h5.attrs["meta_json"] = json.dumps({"Nk": Nk})
        h5.create_dataset("Npw", data=np.int64(Nx * Ny * Nz))
        h5.create_dataset("supercell_size", data=L_size)

        g_lat = h5.create_group("lattice")
        g_lat.create_dataset("A_direct_ang", data=a_vectors)

        g_k = h5.create_group("kmesh")
        g_k.create_dataset("bz/cart_1_perA", data=k_cart)
        g_k.create_dataset("bz/reduced", data=k_cart)

        g_rho = h5.create_group("rho_tilde")
        g_rho.create_dataset("data", data=rho_data,
                             chunks=(1, Nb, 1, Nb, Nx, Ny, Nz))

        g_C = h5.create_group("C_tensor")
        Ca = rng.standard_normal((na, na, na, na))
        Ca_sym = Ca.reshape(na * na, na * na)
        Ca_sym = 0.5 * (Ca_sym + Ca_sym.T)
        g_C.create_dataset("atom_0000", data=Ca_sym.reshape(na, na, na, na))

        g_D = h5.create_group("D_tensor")
        g_D.create_dataset("atom_0000", data=(
            rng.standard_normal((Nk, Nb, Nk, Nb, na, na)) +
            1j * rng.standard_normal((Nk, Nb, Nk, Nb, na, na))
        ))

    return rho_data


class TestLazyRead:
    """Test lazy-read mode where rho is returned as an h5py dataset handle."""

    def test_lazy_returns_h5py_dataset(self, tmp_path):
        """Test that lazy mode returns an h5py.Dataset for rho instead of ndarray.

        Lazy reading avoids loading the full (Nk,Nb,Nk,Nb,Nx,Ny,Nz) tensor
        into memory — critical for large systems where rho can exceed 10 GB.
        """
        filepath = str(tmp_path / "test.h5")
        rho_data = _write_test_hdf5(filepath)

        with PawReader(filepath) as reader:
            inputs = reader.to_calculator_inputs(lazy=True)
            assert isinstance(inputs["rho"], h5py.Dataset)
            block = np.asarray(inputs["rho"][0, :, 0, :, :, :, :])
            np.testing.assert_allclose(block, rho_data[0, :, 0, :, :, :, :])

    def test_lazy_vs_eager_same_data(self, tmp_path):
        """Test that lazy and eager modes return identical rho data.

        This is the round-trip correctness test — if lazy reads produce
        different values, downstream lambda computations will be wrong.
        """
        filepath = str(tmp_path / "test.h5")
        _write_test_hdf5(filepath)

        reader1 = PawReader(filepath)
        inputs1 = reader1.to_calculator_inputs()
        rho_eager = inputs1["rho"]

        with PawReader(filepath) as reader2:
            inputs2 = reader2.to_calculator_inputs(lazy=True)
            rho_lazy = np.asarray(inputs2["rho"][()])

        np.testing.assert_array_equal(rho_eager, rho_lazy)

    def test_context_manager_closes_file(self, tmp_path):
        """Test that the file handle is released after exiting the context manager.

        Unclosed HDF5 files can cause data corruption and resource leaks.
        """
        filepath = str(tmp_path / "test.h5")
        _write_test_hdf5(filepath)

        reader = PawReader(filepath)
        with reader:
            reader.to_calculator_inputs(lazy=True)
            assert reader._h5 is not None
        assert reader._h5 is None

    def test_lazy_rho_works_with_calculator(self, tmp_path):
        """Test that a lazy rho handle produces the same lambda as eager loading.

        End-to-end validation that the OneNormCalculator can consume h5py
        dataset handles transparently.
        """
        from bloch_paw.one_norm import OneNormCalculator

        filepath = str(tmp_path / "test.h5")
        _write_test_hdf5(filepath, Nk=2, Nb=2, grid=(3, 3, 3))

        reader1 = PawReader(filepath)
        inputs1 = reader1.to_calculator_inputs()
        inputs1["h_pq"] = None
        inputs1["kappa_pqrs"] = None
        calc1 = OneNormCalculator(**inputs1)
        lam1 = calc1.lambda_one_norm()

        with PawReader(filepath) as reader2:
            inputs2 = reader2.to_calculator_inputs(lazy=True)
            inputs2["h_pq"] = None
            inputs2["kappa_pqrs"] = None
            calc2 = OneNormCalculator(**inputs2)
            lam2 = calc2.lambda_one_norm()

        np.testing.assert_allclose(lam1, lam2)
