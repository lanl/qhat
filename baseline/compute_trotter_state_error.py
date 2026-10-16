#!/usr/bin/env python3
"""
compute_trotter_state_error.py
------------------------------
Trotter state vector error for the library Hamiltonians:

    err = || psi_trotter(T) - psi_exact(T) ||

Both states start from the HF determinant and are evolved by the SAME operator
(normalized, pi-scaled, identity excluded), matching the Julia state-error path
in trotter_expts.jl Part 2. No 0.5*I centering and no global phase: those belong
to the energy/operator-error path.

n_steps is PER U, where U = exp(-i H_terms) is one unit of evolution, so
delta_t = 1/n_steps and total_steps = T/delta_t = T*n_steps. delta_t is
therefore identical across active spaces.

T is swept, identical for every Hamiltonian -- an independent axis of the study,
not a per-molecule lookup. Total Trotter steps per row:

    T \ delta_t      1      1/2      1/4       (n_steps/U = 1, 2, 4)
    2^6  =   64      64      128      256
    2^8  =  256     256      512     1024
    2^9  =  512     512     1024     2048

27 rows per Hamiltonian (3 T x 3 delta_t x 3 orders). Results are written per
molecule to --out-dir and flushed row by row; rows already in the CSV are
skipped, so an interrupted run loses nothing.

Usage:
    python compute_trotter_state_error.py library/ --out-dir state_errors
    python compute_trotter_state_error.py library/ --molecule Li-Li \
        --T 64,256,512 --n-steps 1,2,4 --order 1,2,4 --n-L 10
"""

import argparse
import csv
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, kron as spkron
from scipy.sparse.linalg import expm_multiply

# ---------------------------------------------------------------------------
# Pauli matrices
# ---------------------------------------------------------------------------
_I = coo_matrix(np.array([[1,0],[0,1]], dtype=complex))
_X = coo_matrix(np.array([[0,1],[1,0]], dtype=complex))
_Y = coo_matrix(np.array([[0,-1j],[1j,0]], dtype=complex))
_Z = coo_matrix(np.array([[1,0],[0,-1]], dtype=complex))
_PAULI = {"I":_I, "X":_X, "Y":_Y, "Z":_Z}


def pauli_string_matrix(ps):
    mat = _PAULI[ps[0]]
    for c in ps[1:]:
        mat = spkron(mat, _PAULI[c], format="csr")
    return mat


# ---------------------------------------------------------------------------
# Parse .dat
# ---------------------------------------------------------------------------

def parse_dat(path):
    metadata, terms = {}, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                if "=" in body:
                    k, _, v = body.partition("=")
                    metadata[k.strip()] = v.strip()
            else:
                parts = line.split()
                if len(parts) == 2:
                    # coefficients may be real ("0.15") or complex
                    # ("9.8e-03+0.0e+00j", as in some generated files); complex()
                    # parses both, and the Hamiltonian is Hermitian so imag~0.
                    z = complex(parts[0])
                    if abs(z.imag) > 1e-9:
                        raise ValueError(
                            f"{path}: nonzero imaginary coefficient {parts[0]!r}")
                    terms.append((z.real, parts[1].upper()))
    return metadata, terms


# ---------------------------------------------------------------------------
# Normalization (Julia-matching: hamiltonian_utils.jl normalize_hamiltonian)
# ---------------------------------------------------------------------------

def normalize_hamiltonian(terms, metadata):
    """
    shift         = identity coefficient
    norm_bound    = one_norm - |shift|         (one_norm from metadata)
    normalization = 2 * max(1, norm_bound)
    Julia uses metadata key 'one-norm of sum of Pauli strings'.
    """
    n      = len(terms[0][1])
    id_str = "I" * n
    key = next((k for k in ("one-norm of sum of Pauli strings",
                            "one-norm of sum of Pauli strings (Hartrees)")
                if k in metadata), None)
    if key is None:
        raise ValueError(
            "metadata lacks 'one-norm of sum of Pauli strings' - refusing to "
            "guess the normalization")
    one_norm = float(metadata[key])
    shift         = next((c for c, ps in terms if ps == id_str), 0.0)
    norm_bound    = one_norm - abs(shift)
    normalization = 2.0 * max(1.0, norm_bound)
    return shift, normalization, one_norm


# ---------------------------------------------------------------------------
# HF state (matches hamgen.get_initial_state)
# ---------------------------------------------------------------------------

def _encoder_bk(n_modes):
    """
    Binary-tree encoder for the Bravyi-Kitaev transform. Vendored verbatim from
    openfermion.transforms.opconversions.binary_codes, which is the same
    function hamgen.py imports, so this script stays numpy/scipy-only.
    """
    reps = int(np.ceil(np.log2(n_modes)))
    mtx  = np.array([[1, 0], [1, 1]])
    for repetition in np.arange(1, reps + 1):
        mtx = np.kron(np.eye(2, dtype=int), mtx)
        for column in np.arange(0, 2 ** repetition):
            mtx[2 ** (repetition + 1) - 1, column] = 1
    return mtx[0:n_modes, 0:n_modes]


def hf_state(metadata, terms, mapping="jw"):
    """
    HF determinant as a computational basis state, following
    hamgen.get_initial_state: the lowest n_occ spin orbitals are occupied, JW
    stores that occupation directly, BK stores (encoder_bk @ occ) mod 2. The
    occupation vector is packed MSB-first, matching the kron order in
    pauli_string_matrix, so qubit 0 is Pauli-string position 0.

    n_occ comes from metadata, same key as Julia trotter_expts.jl line 44.
    """
    n   = len(terms[0][1])
    key = "number of active, occupied, single-occupancy orbitals"
    if key not in metadata:
        raise KeyError(
            f"HF state requires metadata key '{key}' (matches Julia). "
            f"Refusing to guess n_occ.")
    n_occ = int(metadata[key])

    occ = np.zeros(n, dtype=int)
    occ[0:n_occ] = 1

    m = mapping.lower()
    if m in ("jw", "jordan-wigner"):
        occ_qubit = occ
    elif m in ("bk", "bravyi-kitaev"):
        occ_qubit = (_encoder_bk(n) @ occ) % 2
    else:
        raise ValueError(f"unknown fermion-to-qubit mapping {mapping!r}")

    idx      = sum(int(occ_qubit[i]) << (n - 1 - i) for i in range(n))
    psi      = np.zeros(2 ** n, dtype=complex)
    psi[idx] = 1.0
    return psi


# ---------------------------------------------------------------------------
# Term ordering
#
# Trotter error depends strongly on term order, so it is isolated in one
# swappable function. An ordering takes and returns a list of
# (coeff, pauli_string) pairs, already normalized, pi-scaled and identity-free:
#
#     mats = make_qpe_mats(terms, normalization, ordering=my_ordering)
# ---------------------------------------------------------------------------

def order_ascending_lex(terms):
    """
    Ascending signed coefficient, ties broken by ascending lexicographic Pauli
    string. The direction matches Julia's sort!(H_terms, by = x -> x[1]); the
    tie-break is an addition, since Julia's bare sort leaves equal coefficients
    (common among symmetry-related terms) in an unspecified order.
    """
    return sorted(terms, key=lambda t: (t[0], t[1]))


# ---------------------------------------------------------------------------
# QPE Trotter terms (Julia-matching build_hamiltonian_terms:
#   coeff = (c/normalization)*pi, identity dropped, then ordered)
# ---------------------------------------------------------------------------

def make_qpe_mats(terms, normalization, ordering=order_ascending_lex):
    """Drop the identity, rescale to pi*(c/normalization), order, build matrices."""
    n      = len(terms[0][1])
    id_str = "I" * n
    scaled = [((c / normalization) * np.pi, ps)
              for c, ps in terms if ps != id_str]
    scaled = ordering(scaled)
    return [(c, pauli_string_matrix(ps)) for c, ps in scaled]


def mats_to_sparse(mats, dim):
    """
    Sum (coeff * P) into one sparse Hermitian matrix -- the same operator that
    generates the Trotter evolution, so both states share a generator and no
    global phase appears. Forming it densely enough is fine at <= 12 qubits.
    """
    H = coo_matrix((dim, dim), dtype=complex).tocsr()
    for c, P in mats:
        H = H + c * P
    return H


def apply_rotation(P, coeff_dt, psi):
    # exp(-i*coeff_dt*P)|psi> = cos(coeff_dt)|psi> - i sin(coeff_dt) P|psi>
    # Julia: apply_pauli_rotation!(psi, P, coeff, dt) with theta = coeff*dt.
    return np.cos(coeff_dt) * psi - 1j * np.sin(coeff_dt) * (P @ psi)


# ---------------------------------------------------------------------------
# Trotter state vector evolution (single vector, not full unitary).
# These mirror Julia statevector_simulators.jl term for term.
# ---------------------------------------------------------------------------

def trotter1_state(mats, psi0, dt, n_steps):
    """Julia first_order_trotter_statevec."""
    psi = psi0.copy()
    for _ in range(n_steps):
        for c, P in mats:
            psi = apply_rotation(P, c * dt, psi)
    return psi


def _s2_step(mats, psi, dt):
    """One Strang step: 1..m-1 at dt/2, m at dt, m-1..1 at dt/2."""
    for c, P in mats[:-1]:
        psi = apply_rotation(P, c * dt / 2, psi)
    c_l, P_l = mats[-1]
    psi = apply_rotation(P_l, c_l * dt, psi)
    for c, P in reversed(mats[:-1]):
        psi = apply_rotation(P, c * dt / 2, psi)
    return psi


def trotter2_state(mats, psi0, dt, n_steps):
    """Julia second_order_trotter_statevec."""
    psi = psi0.copy()
    for _ in range(n_steps):
        psi = _s2_step(mats, psi, dt)
    return psi


def trotter4_state(mats, psi0, dt, n_steps):
    """Julia fourth_order_trotter_statevec: S2(p dt)^2 S2((1-4p)dt) S2(p dt)^2."""
    p = 1.0 / (4.0 - 4.0**(1.0/3.0))
    psi = psi0.copy()
    for _ in range(n_steps):
        psi = _s2_step(mats, psi, p*dt)
        psi = _s2_step(mats, psi, p*dt)
        psi = _s2_step(mats, psi, (1-4*p)*dt)
        psi = _s2_step(mats, psi, p*dt)
        psi = _s2_step(mats, psi, p*dt)
    return psi


# ---------------------------------------------------------------------------
# Find .dat files
# ---------------------------------------------------------------------------

_PAT     = re.compile(r"^(?P<stub>.+)_as-(?P<occ>\d{3})-(?P<vac>\d{3})_(?P<mapping>[a-zA-Z]+)\.dat$")
_L_PAT   = re.compile(r"_(?P<L>\d+\.\d+)_")
_MOL_PAT = re.compile(r"^(?P<mol>[A-Za-z]+(?:-[A-Za-z]+)?)_")

LIBRARY_MOLECULES = {
    "H-H","He-H","He-He","Li-H","Li-Li",
    "Be-H","Be-Be","B-H","B-B",
    # monatomic systems, stub '<X>_atom_<basis>' -- no bond length
    "H","He","Li","Be","B","C","N","O","F","Ne",
}


def find_dat_files(root, molecule, basis, mappings, max_qubits, L_min, L_max, active=None):
    results = []
    for p in sorted(Path(root).rglob("*.dat")):
        m = _PAT.fullmatch(p.name)
        if not m:
            continue
        n_q = int(m.group("occ")) + int(m.group("vac"))
        if n_q > max_qubits:
            continue
        if basis   and basis.lower()   not in p.name.lower():
            continue
        active_str = f"{m.group('occ')}-{m.group('vac')}"
        if active and active_str != active:
            continue
        if mappings and m.group("mapping").lower() not in mappings:
            continue
        lm  = _L_PAT.search(m.group("stub"))
        mm  = _MOL_PAT.match(m.group("stub"))
        mol = mm.group("mol") if mm else "unknown"
        if molecule and mol.lower() != molecule.lower():
            continue
        if not molecule and mol not in LIBRARY_MOLECULES:
            continue
        # Monatomic systems have no bond length. nan sorts, formats with
        # .4f, survives the L window (nan comparisons are False) and
        # round-trips through the CSV as 'nan', so the resume key still
        # matches -- none of which None would do.
        L_val = float(lm.group("L")) if lm else float("nan")
        if L_min is not None and L_val < L_min:
            continue
        if L_max is not None and L_val > L_max:
            continue
        results.append({
            "path"    : p,
            "molecule": mol,
            "L"       : L_val,
            "basis"   : m.group("stub").split("_")[-1],
            "active"  : f"{m.group('occ')}-{m.group('vac')}",
            "mapping" : m.group("mapping").lower(),
            "n_qubits": n_q,
        })
    results.sort(key=lambda x: (x["molecule"], x["basis"], x["active"], x["L"]))
    return results


def subsample_L(files, n_L):
    """
    Keep at most n_L bond lengths per (molecule, basis, mapping, active) family,
    evenly
    spaced over the sorted L range with both endpoints always kept, so the
    near-equilibrium region and the dissociation tail are both represented.
    n_L=None -> no subsampling.
    """
    if n_L is None:
        return files
    from collections import defaultdict
    groups = defaultdict(list)
    for f in files:
        groups[(f["molecule"], f["basis"], f["mapping"], f["active"])].append(f)

    kept = []
    for key, g in groups.items():
        g = sorted(g, key=lambda x: x["L"])
        m = len(g)
        if m <= n_L:
            kept.extend(g)
            continue
        # evenly spaced indices over [0, m-1], inclusive of both endpoints
        idx = sorted(set(int(round(i)) for i in np.linspace(0, m - 1, n_L)))
        kept.extend(g[i] for i in idx)
    kept.sort(key=lambda x: (x["molecule"], x["basis"], x["active"], x["L"]))
    return kept


# ---------------------------------------------------------------------------
# CSV fields
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "molecule", "L", "basis", "active", "mapping", "n_qubits",
    "T", "n_steps", "total_steps", "delta_t", "dt_norm", "order",
    "err",           # state error ||psi_trotter - psi_exact||  (expm reference)
    "seconds",       # wall time for this row
    "one_norm", "normalization",
]


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def provenance_lines(args, basis_tag, mapping_tag):
    """
    Conditions that produced the rows, written as leading '#' comment lines so
    the CSV is self-documenting. A challenger comparing a new Trotter technique
    must match all of these; only then is a lower `err` a real improvement.
    Load the data with pandas.read_csv(path, comment='#').
    """
    return [
        "BASELINE: Trotter state error for QPE Hamiltonian simulation.",
        "",
        "metric        : err = || psi_trotter - psi_exact ||_2",
        "initial state : Hartree-Fock |11..100..0> (top n_occ bits, JW)",
        "reference     : scipy expm_multiply, exp(-i H T) psi0",
        "generator     : H_terms = pi * (H - shift*I) / normalization,",
        "                identity term excluded",
        "normalization : 2 * max(1, one_norm - |shift|)  [Julia "
        "hamiltonian_utils.jl]",
        "term ordering : ascending signed coefficient, ties broken by ascending",
        "                lexicographic Pauli string (deterministic; the",
        "                ascending direction follows Julia, the tie-break is an",
        "                addition)",
        "delta_t       : 1 / n_steps, where n_steps is PER U and",
        "                U = exp(-i H_terms) is one unit of evolution",
        "total_steps   : T / delta_t = T * n_steps",
        "dt_norm       : delta_t * sum|c_j|; equals pi/2, pi/4, pi/8 wherever",
        "                the shifted one-norm is at least 1, which fixes the",
        "                dimensionless step size across systems",
        "T             : swept, identical for every Hamiltonian. NOT derived",
        "                from any per-molecule circuit-depth table; it is an",
        "                independent axis of the study.",
        f"mapping       : {mapping_tag}   (HF determinant encoded per mapping:",
        "                JW stores the occupation directly, BK stores",
        "                (encoder_bk @ occ) mod 2, as in hamgen.get_initial_state)",
        f"basis         : {basis_tag}",
        f"active space  : {args.active if args.active else 'all'}",
        f"T swept       : {args.T_values}",
        f"n_steps swept : {args.n_steps}   (delta_t = "
        f"{', '.join(str(round(1.0/int(s), 6)) for s in args.n_steps.split(','))})",
        f"orders        : {args.order}",
        f"bond lengths  : {args.n_L if args.n_L else 'all'} per active space, "
        "evenly spaced over the sorted L range (deterministic)",
        "",
        "Generated by compute_trotter_state_error.py.",
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(args):
    n_steps_list = [int(s) for s in args.n_steps.split(",")]
    T_list       = [int(t) for t in args.T_values.split(",")]
    mappings     = [m.strip().lower() for m in args.mapping.split(",")]
    bad_map = sorted(set(mappings) - {"jw", "bk"})
    if bad_map:
        print(f"ERROR: unsupported mappings {bad_map}; use jw, bk, or jw,bk.")
        sys.exit(1)
    orders   = [int(o) for o in args.order.split(",")]
    invalid = sorted(set(orders) - {1, 2, 4})
    if invalid:
        print(f"ERROR: unsupported Trotter orders {invalid}; use 1, 2, or 4.")
        sys.exit(1)
    out_dir  = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if any(s <= 0 for s in n_steps_list):
        print("ERROR: all n_steps must be positive.")
        sys.exit(1)
    if any(t <= 0 for t in T_list):
        print("ERROR: all T must be positive.")
        sys.exit(1)

    molecules = [args.molecule] if args.molecule else sorted(LIBRARY_MOLECULES)

    print(f"T swept   : {T_list}   (same for every Hamiltonian)")
    print(f"n_steps/U : {n_steps_list}   -> delta_t = "
          f"{[round(1.0/n, 6) for n in n_steps_list]}  (same for every "
          f"active space)")
    print(f"Orders    : {orders}")
    print(f"Mappings  : {mappings}")
    print()
    print("  total Trotter steps = T / delta_t:")
    hdr = "    " + "T".rjust(8) + "".join(
        f"{'dt=' + str(round(1.0/n, 6)):>16}" for n in n_steps_list)
    print(hdr)
    for T in T_list:
        print("    " + f"{T:>8}" + "".join(f"{T*n:>16,}" for n in n_steps_list))
    print()
    print(f"rows/Ham  : {len(T_list)*len(n_steps_list)*len(orders)}")
    print(f"Max qubits: {args.max_qubits}")
    print(f"Molecules : {molecules}")
    print(f"Output dir: {out_dir}\n")

    for molecule in molecules:
        print(f"\n{'='*60}")
        print(f"Molecule: {molecule}")
        print(f"{'='*60}")
        run_molecule(args, molecule, T_list, n_steps_list, orders, mappings, out_dir)

    print(f"\nAll done. CSVs saved under {out_dir}/<basis>/<mapping>/")


def run_molecule(args, molecule, T_list, n_steps_list, orders, mappings, out_root):
    """
    Find every .dat for this molecule, then split by (basis, mapping) and write
    one CSV per combination under out_root/<basis>/<mapping>/. Each CSV holds a
    single basis and a single mapping, so its provenance header describes
    exactly its own contents.
    """
    files = find_dat_files(
        args.root, molecule, args.basis, mappings,
        args.max_qubits, args.L_min, args.L_max, args.active
    )
    if not files:
        print(f"  No .dat files found for {molecule}.")
        return

    n_found = len(files)
    files = subsample_L(files, args.n_L)
    if args.n_L is not None:
        print(f"Found {n_found} files; kept {len(files)} after L subsample "
              f"(<= {args.n_L} per basis/mapping/active space)\n")
    else:
        print(f"Found {len(files)} files\n")

    from collections import defaultdict
    groups = defaultdict(list)
    for f in files:
        groups[(f["basis"], f["mapping"])].append(f)

    for basis, mapping in sorted(groups):
        out_csv = out_root / basis / mapping / f"state_errors_{molecule}.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        print(f"--- {molecule}  basis={basis}  mapping={mapping}  "
              f"({len(groups[(basis, mapping)])} files) -> {out_csv}")
        run_group(args, molecule, basis, mapping, groups[(basis, mapping)],
                  T_list, n_steps_list, orders, out_csv)


def run_group(args, molecule, basis_tag, mapping_tag, files,
              T_list, n_steps_list, orders, out_csv):
    existing = set()
    if out_csv.exists():
        with open(out_csv) as f:
            # skip the leading '#' provenance lines before the column header
            reader = csv.DictReader(line for line in f if not line.startswith("#"))
            existing_fields = reader.fieldnames or []
            if set(existing_fields) != set(CSV_FIELDS):
                print("WARNING: existing CSV has different columns - starting fresh.")
                out_csv.unlink()
            else:
                for row in reader:
                    key = (row["molecule"], row["L"], row["basis"],
                           row["active"], row["mapping"],
                           row["T"], row["n_steps"], row["order"])
                    existing.add(key)
                print(f"Resuming - {len(existing)} existing rows found.\n")

    write_header = not out_csv.exists()
    out_f  = open(out_csv, "a", newline="")
    writer = csv.DictWriter(out_f, fieldnames=CSV_FIELDS)
    if write_header:
        for line in provenance_lines(args, basis_tag, mapping_tag):
            out_f.write(f"# {line}\n")
        writer.writeheader()

    processed = skipped = 0
    fn_map = {1: trotter1_state, 2: trotter2_state, 4: trotter4_state}

    for info in files:
        mol      = info["molecule"]
        L        = info["L"]
        basis    = info["basis"]
        active   = info["active"]
        mapping  = info["mapping"]
        n_qubits = info["n_qubits"]
        dim      = 2 ** n_qubits

        all_done = all(
            (mol, str(L), basis, active, mapping,
             str(T), str(s), str(o)) in existing
            for T in T_list for s in n_steps_list for o in orders
        )
        if all_done:
            skipped += 1
            continue

        print(f"\n[{processed+skipped+1}/{len(files)}] "
              f"{mol}  L={L:.4f}  {active}  {n_qubits}q")

        try:
            meta, terms = parse_dat(info["path"])
            shift, normalization, one_norm = normalize_hamiltonian(terms, meta)

            if one_norm <= 0:
                print(f"  WARNING: one_norm={one_norm} - skipping")
                skipped += 1
                continue

            psi0 = hf_state(meta, terms, mapping)

            mats = make_qpe_mats(terms, normalization,
                                 ordering=order_ascending_lex)
            print(f"  one_norm={one_norm:.4f}  norm={normalization:.4f}"
                  f"  n_terms={len(mats)}")

            # Same operator as the Trotter steps: assembled once per
            # Hamiltonian, exponentiated once per T below.
            H_mats = mats_to_sparse(mats, dim)
            l1 = sum(abs(c) for c, _ in mats)
            print(f"  dt_norm at delta_t=1 is {l1:.6f}"
                  f"{'  (= pi/2)' if abs(l1 - np.pi/2) < 1e-9 else ''}")

        except Exception as e:
            print(f"  ERROR: {e}")
            skipped += 1
            continue

        for T in T_list:
            # Skip the reference entirely if every row at this T is already in
            # the CSV -- expm_multiply is the expensive part of a resumed run.
            if all((mol, str(L), basis, active, mapping,
                    str(T), str(s), str(o)) in existing
                   for s in n_steps_list for o in orders):
                continue
            try:
                psi_exact = expm_multiply(-1j * float(T) * H_mats, psi0)
            except Exception as e:
                print(f"  ERROR building reference at T={T}: {e}")
                continue

            for order in orders:
                fn = fn_map[order]
                for n_steps in n_steps_list:
                    key = (mol, str(L), basis, active, mapping,
                           str(T), str(n_steps), str(order))
                    if key in existing:
                        continue
                    try:
                        delta_t = 1.0 / n_steps
                        total_steps = T * n_steps   # both powers of 2: exact
                        dt_norm = delta_t * l1
                        t0 = time.time()
                        # NO global phase: both states share the generator.
                        psi_tr = fn(mats, psi0, delta_t, total_steps)
                        err = float(np.linalg.norm(psi_tr - psi_exact))
                        elapsed = time.time() - t0

                        writer.writerow({
                            "molecule"      : mol,
                            "L"             : L,
                            "basis"         : basis,
                            "active"        : active,
                            "mapping"       : mapping,
                            "n_qubits"      : n_qubits,
                            "T"             : T,
                            "n_steps"       : n_steps,
                            "total_steps"   : total_steps,
                            "delta_t"       : delta_t,
                            "dt_norm"       : dt_norm,
                            "order"         : order,
                            "err"           : err,
                            "seconds"       : round(elapsed, 3),
                            "one_norm"      : one_norm,
                            "normalization" : normalization,
                        })
                        out_f.flush()

                        print(f"  order={order}  T={T:<6d} "
                              f"n_steps/U={n_steps:4d}  dt={delta_t:.4g}  "
                              f"total={total_steps:8d}  err={err:.3e}  "
                              f"({elapsed:.1f}s)")

                    except Exception as e:
                        print(f"  ERROR order={order} n_steps={n_steps}: {e}")

        processed += 1

    out_f.close()
    print(f"\nDone. {processed} processed, {skipped} skipped -> {out_csv.name}")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Compute Trotter state vector errors over a grid of "
                    "evolution times T and step sizes delta_t."
    )
    ap.add_argument("root",            help="Library root e.g. library/")
    ap.add_argument("--out-dir",       default="state_errors", dest="out_dir")
    ap.add_argument("--molecule",      default=None, help="Single molecule (default: all 9 library molecules)")
    ap.add_argument("--basis",         default=None,
                    help="hgbs-5 or sto-6g (default: all bases)")
    ap.add_argument("--active",        default=None,
                    help="filter to one active space, e.g. 004-008 "
                         "(default: all active spaces)")
    ap.add_argument("--mapping",       default="jw,bk",
                    help="Fermion-to-qubit mappings to include: jw, bk, or "
                         "jw,bk (default: both). Both land in one CSV, "
                         "distinguished by the mapping column.")
    ap.add_argument("--T",             default="64,256,512", dest="T_values",
                    help="Total evolution times to sweep, in applications of "
                         "U = exp(-i H_terms) (default: 64,256,512, i.e. "
                         "2^6, 2^8, 2^9). The same values are used for every "
                         "Hamiltonian; T is an independent axis of the study, "
                         "not a per-molecule property.")
    ap.add_argument("--n-steps",       default="1,2,4", dest="n_steps",
                    help="Trotter steps PER U, where U = exp(-i H_terms) is "
                         "one unit of evolution (default: 1,2,4, i.e. "
                         "delta_t = 1, 1/2, 1/4). delta_t = 1/n_steps, so it "
                         "is identical across active spaces. Total steps for "
                         "a row is T*n_steps = T/delta_t, which sets the cost.")
    ap.add_argument("--order",         default="1,2,4")
    ap.add_argument("--max-qubits",    type=int, default=12, dest="max_qubits")
    ap.add_argument("--L-min",         type=float, default=None, dest="L_min")
    ap.add_argument("--L-max",         type=float, default=None, dest="L_max")
    ap.add_argument("--n-L",           type=int,   default=None, dest="n_L",
                    help="Subsample to at most this many bond lengths per "
                         "active space, evenly spaced over the L range "
                         "(endpoints always kept). Default: use all L.")
    return ap.parse_args()


if __name__ == "__main__":
    run(parse_args())