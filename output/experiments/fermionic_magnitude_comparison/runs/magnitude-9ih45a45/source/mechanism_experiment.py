#!/usr/bin/env python3
"""Fixed-input ordering mechanism audit. See PROTOCOL.md before interpreting."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

# Set before NumPy/QHAT imports; all runs are single-threaded.
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
sys.dont_write_bytecode = True

REPO = Path('/Users/albertlee0125/Repos/qhat')
PRIOR = REPO / 'codex_first_experiment/repro_runs/first-repro-y23sjf9v'
BASE = {
    'fermionic_signed': 'fermionic_signed_coefficient_lexicographic',
    'jw_signed': 'signed_coefficient_lexicographic',
    'jw_magnitude': 'jw_magnitude_descending_lexicographic',
}
SEEDS = [1101, 1102, 1103, 1104, 1105]
FLOOR = 1e-12


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def obj_sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def save(path, obj):
    with Path(path).open('x') as stream:
        json.dump(obj, stream, indent=2, allow_nan=False)
        stream.write('\n')


def git(*args):
    return subprocess.check_output(['git', '-C', str(REPO), *args], text=True).strip()


def source_hashes():
    return {p: sha(REPO / p) for p in git('ls-files', '--', '*.py', 'pyproject.toml').splitlines()}


def init_numerics(cache):
    global np, sparse, expm_multiply, baseline, matrixfree, ablation
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    os.environ['MPLCONFIGDIR'] = str(cache / 'mpl')
    os.environ['NUMBA_CACHE_DIR'] = str(cache / 'numba')
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import expm_multiply
    sys.path[:0] = [str(REPO.parent), str(REPO)]
    from qhat.analysis import benchmark_b2_signed_coefficient_baseline as baseline
    from qhat.analysis import benchmark_b2_active_spaces_matrix_free as matrixfree
    from qhat.analysis import benchmark_fermionic_structure_ablation as ablation
    for module in (baseline, matrixfree, ablation):
        if not Path(module.__file__).resolve().is_relative_to(REPO):
            raise RuntimeError('Imported source outside the selected QHAT repository')


def affine_space(terms, hf_index):
    """Exact closure under every flip; coordinate bit k is basis[k]."""
    pivots = {}
    for _, flip, _, _ in terms:
        x = int(flip)
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    basis = [pivots[p] for p in sorted(pivots)]
    states = [int(hf_index)]
    for b in basis:
        states += [x ^ b for x in states]
    indices = np.asarray(states, dtype=np.int64)
    coordinates = {int(x) ^ int(hf_index): i for i, x in enumerate(indices)}
    reduced = []
    for c, x, z, y in terms:
        zz = sum(((b & int(z)).bit_count() % 2) << k for k, b in enumerate(basis))
        cc = c * (-1 if (int(hf_index) & int(z)).bit_count() % 2 else 1)
        reduced.append((cc, np.uint64(coordinates[int(x)]), np.uint64(zz), y))
    return indices, reduced, basis


def apply_p(state, term, parity):
    _, x, z, y = term
    src = np.arange(state.size, dtype=np.int64) ^ int(x)
    signs = 1 - 2 * parity[src & int(z)]
    return (1j ** int(y)) * signs * state[src]


def assemble_hamiltonian(terms, dimension):
    """Assemble exactly the retained Pauli sum; discard only exact zero entries."""
    coords = np.arange(dimension, dtype=np.int64)
    parity = np.array([int(x).bit_count() % 2 for x in coords], dtype=np.int8)
    grouped = {}
    for c, x, z, y in terms:
        x = int(x)
        if x not in grouped:
            grouped[x] = np.zeros(dimension, dtype=np.complex128)
        grouped[x] += c * (1j ** int(y)) * (1 - 2 * parity[coords & int(z)])
    rows, cols, data = [], [], []
    for x, values in grouped.items():
        keep = np.flatnonzero(values != 0)
        rows.append(keep ^ x)
        cols.append(keep)
        data.append(values[keep])
    H = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(dimension, dimension)).tocsr()
    if sparse.linalg.norm(H - H.getH()) > 1e-10:
        raise RuntimeError('Assembled Hamiltonian is not Hermitian')
    return H, parity


def evolve(state, terms, duration, steps=1):
    return matrixfree.evolve_trotter_state(state, terms, 1, steps, duration, 2**30)[0]


def leading_defect(terms, initial, H, parity):
    """Coefficient of dt^2 in S(dt)|HF>-exp(-iHdt)|HF>."""
    first = np.zeros_like(initial)
    second = np.zeros_like(initial)
    for term in terms:
        c = term[0]
        second -= 1j * c * apply_p(first, term, parity)
        second -= 0.5 * c * c * initial
        first -= 1j * c * apply_p(initial, term, parity)
    return second + 0.5 * (H @ (H @ initial))


def norm2(v):
    return float(np.vdot(v, v).real)


def physical_metrics(exact, approximate, number_mask, spin_mask):
    ne, na = norm2(exact), norm2(approximate)
    if abs(ne-1) > 1e-9 or abs(na-1) > 1e-9:
        raise RuntimeError(f'Excessive state norm drift: exact={ne}, approximate={na}')
    phi = exact / math.sqrt(ne)
    overlap = np.vdot(phi, approximate)
    transverse = approximate - phi * overlap
    I = norm2(transverse) / na
    E = I / (1 + math.sqrt(max(0.0, 1-I)))
    in_weight = norm2(approximate[spin_mask])
    exact_in = exact[spin_mask]
    exact_in = exact_in / np.linalg.norm(exact_in)
    approx_in = approximate[spin_mask]
    conditional_residual = approx_in - exact_in * np.vdot(exact_in, approx_in)
    out = {
        'one_minus_overlap': E,
        'infidelity': I,
        'naive_normalized_one_minus_overlap': 1 - float(abs(overlap)) / math.sqrt(na),
        'inside_sector_error': norm2(transverse[spin_mask]) / na,
        'outside_sector_error': norm2(transverse[~spin_mask]) / na,
        'particle_number_leakage': norm2(approximate[~number_mask]) / na,
        'spin_sector_leakage': norm2(approximate[~spin_mask]) / na,
        'conditional_sector_infidelity': norm2(conditional_residual) / in_weight,
        'exact_spin_leakage': norm2(exact[~spin_mask]) / ne,
        'exact_norm_drift': abs(ne-1),
        'trotter_norm_drift': abs(na-1),
        'below_analysis_floor': E <= FLOOR,
    }
    if out['exact_spin_leakage'] > 1e-18:
        raise RuntimeError('Exact Hamiltonian leaves the intended symmetry sector')
    return out


def decompose(trajectory, terms, duration, approximate):
    """O(r) exact projected-defect decomposition and vector reconstruction."""
    steps = len(trajectory)-1
    dt = duration / steps
    phi = trajectory[-1] / np.linalg.norm(trajectory[-1])
    bra = phi.copy()
    reverse_terms = tuple(reversed(terms))
    defects = [None] * steps
    trace = [None] * steps
    overlaps = []
    raw_power = projected_power = 0.0
    for j in range(steps-1, -1, -1):
        delta = evolve(trajectory[j], terms, dt) - trajectory[j+1]
        defects[j] = delta
        a = np.vdot(bra, delta)
        power = norm2(delta)
        projected = power - float(abs(a)**2)
        if projected < -1e-16:
            raise RuntimeError('Invalid negative projected local power')
        projected = max(0.0, projected)
        raw_power += power
        projected_power += projected
        overlaps.append(a)
        trace[j] = {'step': j, 'time': j*dt, 'local_defect_norm2': power,
                    'transported_projected_norm2': projected,
                    'overlap_contribution_real': float(a.real),
                    'overlap_contribution_imag': float(a.imag)}
        if j:
            bra = evolve(bra, reverse_terms, -dt)
    reconstruction = np.zeros_like(phi)
    for delta in defects:
        reconstruction = evolve(reconstruction, terms, dt) + delta
    direct_delta = approximate - trajectory[-1]
    residual = float(np.linalg.norm(reconstruction-direct_delta))
    actual_overlap = np.vdot(phi, approximate)
    reconstructed_overlap = np.linalg.norm(trajectory[-1]) + sum(overlaps)
    overlap_residual = float(abs(actual_overlap-reconstructed_overlap))
    if residual > 1e-10 or overlap_residual > 1e-10:
        raise RuntimeError(f'Defect reconstruction failed: {residual}, {overlap_residual}')
    final_norm2 = norm2(approximate)
    qerror = approximate - phi * actual_overlap
    I = norm2(qerror) / final_norm2
    D = projected_power / final_norm2
    X = I-D
    return {
        'projected_local_power': D,
        'projected_interference': X,
        'coherence_factor': I/D if D else None,
        'raw_local_power': raw_power / final_norm2,
        'raw_final_error_norm2': norm2(direct_delta) / final_norm2,
        'state_reconstruction_residual': residual,
        'overlap_reconstruction_residual': overlap_residual,
        'reconstructed_projected_infidelity': norm2(reconstruction-phi*np.vdot(phi,reconstruction))/final_norm2,
    }, trace


def canonical_terms(coefficients):
    return [[k, float(complex(v).real).hex(), float(complex(v).imag).hex()]
            for k, v in sorted(coefficients.items())]


def prepare_case(case, tensor, old, controls):
    from openfermion import get_fermion_operator, jordan_wigner
    interaction, n = baseline.load_interaction_operator(tensor)
    metadata = baseline.parse_case_metadata(tensor, n)
    fermion = baseline.clean_fermion_operator(get_fermion_operator(interaction), 1e-12)
    jw = jordan_wigner(fermion)
    jw.compress(abs_tol=1e-12)
    coeff = {k:v for k,v in jw.terms.items() if k and abs(v)>1e-12}
    keys = list(coeff)
    orders = baseline.build_deterministic_orderings(fermion, keys, coeff, n, 1e-12)
    key_index = {k:i for i,k in enumerate(keys)}
    indices = {short: [key_index[k] for k in orders[name]] for short,name in BASE.items()}
    parents = baseline.build_hermitian_fermion_terms(fermion, 1e-12)
    mapping = ablation.precompute_fermion_to_pauli_indices(parents, coeff, key_index, 1e-12)
    parent_order = baseline.fermionic_term_order_indices(parents, 'signed_ascending', 1e-12)
    buckets, owners, fallback = ablation.build_reference_parent_buckets(parent_order, mapping, len(keys))
    if ablation.flatten_reference_buckets(buckets, fallback) != indices['fermionic_signed']:
        raise RuntimeError('Frozen bucket ownership does not reproduce baseline ordering')
    if controls:
        indices['round_robin'] = ablation.round_robin_descendants(buckets, fallback)[0]
        indices['within_blocks_20260917'] = ablation.randomized_within_blocks(
            buckets, fallback, np.random.default_rng(20260917))[0]
        for seed in SEEDS:
            indices[f'blocks_random_{seed}'] = ablation.randomized_block_order(
                buckets, fallback, np.random.default_rng(seed))[0]
    for name, order in indices.items():
        if sorted(order) != list(range(len(keys))):
            raise RuntimeError(f'Not a permutation: {name}')
    fingerprints = {'identity_free_jw_sha256': obj_sha(canonical_terms(coeff)),
        'hermitian_parents_sha256': obj_sha(sorted([canonical_terms(p.operator.terms) for p in parents], key=repr)),
        'physical_order_sha256': {name: obj_sha([[keys[i], complex(coeff[keys[i]]).real.hex(),
                                                complex(coeff[keys[i]]).imag.hex()] for i in order])
                                  for name,order in indices.items()}}
    if old:
        for key in ('identity_free_jw_sha256','hermitian_parents_sha256'):
            if fingerprints[key] != old['fingerprints'][key]:
                raise RuntimeError(f'Changed frozen operator {case}: {key}')
        for short,name in BASE.items():
            if fingerprints['physical_order_sha256'][short] != old['fingerprints']['physical_order_sha256'][name]:
                raise RuntimeError('Physical baseline order changed')
    original_terms = matrixfree.compile_ordered_terms(keys, coeff, n, 1e-12)
    initial_full = matrixfree.build_hartree_fock_state(n, metadata.active_occupied)
    hf = int(np.argmax(abs(initial_full)))
    states, reduced, basis = affine_space(original_terms, hf)
    H, parity = assemble_hamiltonian(reduced, len(states))
    initial = np.zeros(len(states), dtype=np.complex128)
    initial[0] = 1
    number_mask = np.array([int(x).bit_count()==metadata.active_occupied for x in states])
    spin_set = set(map(int, matrixfree.spin_sector_basis_indices(n, metadata.active_occupied)))
    spin_mask = np.array([int(x) in spin_set for x in states])
    terms = {name: [reduced[i] for i in order] for name,order in indices.items()}
    # Compare three full-space baseline slices with exact compressed execution.
    validation = {}
    for name in BASE:
        full = evolve(initial_full, [original_terms[i] for i in indices[name]], 0.01)
        small = evolve(initial, terms[name], 0.01)
        embedded = np.zeros_like(full)
        embedded[states] = small
        err = float(np.linalg.norm(full-embedded))
        if err > 1e-10:
            raise RuntimeError(f'Full-space validation failed: {case}/{name}/{err}')
        validation[name] = err
    # Independent Hamiltonian action test on a random vector in the closed space.
    rng = np.random.default_rng(17092026)
    probe = rng.normal(size=len(states)) + 1j*rng.normal(size=len(states))
    probe /= np.linalg.norm(probe)
    expected = sum((term[0]*apply_p(probe,term,parity) for term in reduced), np.zeros_like(probe))
    h_residual = float(np.linalg.norm(H@probe-expected))
    if h_residual > 1e-10:
        raise RuntimeError('Hamiltonian action validation failed')
    leak_norm = float(sparse.linalg.norm(H[~spin_mask][:,spin_mask]))
    if leak_norm > 1e-9:
        raise RuntimeError('Hamiltonian spin-sector invariance check failed')
    leading = {}
    for name, compiled in terms.items():
        v = leading_defect(compiled, initial, H, parity)
        qv = v-initial*np.vdot(initial,v)
        leading[name] = {'bch_hf_norm': 2*float(np.linalg.norm(v)),
                         'projected_bch_hf_norm': 2*float(np.linalg.norm(qv)),
                         'leading_defect_inside_norm2': norm2(qv[spin_mask]),
                         'leading_defect_outside_norm2': norm2(qv[~spin_mask])}
        if old and name in BASE:
            ref = next(r['bch2_hf_state_norm'] for r in old['rows'] if r['ordering']==BASE[name])
            if not math.isclose(leading[name]['bch_hf_norm'], ref, rel_tol=1e-7, abs_tol=1e-9):
                raise RuntimeError(f'Initial BCH validation failed {case}/{name}')
    info = {'case':case, 'n_qubits':n, 'n_electrons':metadata.active_occupied,
            'pauli_terms':len(keys), 'parents':len(parents), 'parent_blocks':len(buckets),
            'fallback_count':len(fallback), 'closed_dimension':len(states), 'flip_basis':basis,
            'hamiltonian_nnz':H.nnz, 'fingerprints':fingerprints,
            'full_space_slice_residuals':validation, 'hamiltonian_action_residual':h_residual,
            'hamiltonian_sector_offblock_frobenius':leak_norm,
            'jw_identity_coefficient':float(complex(jw.terms.get((),0)).real),
            'old_exact_minus_trotter_scalar':old.get('exact_minus_trotter_scalar') if old else None,
            'leading_defects':leading}
    snapshot = {'coefficients':canonical_terms(coeff), 'orders':indices,
                'parent_buckets':buckets, 'fallback':fallback, 'closed_basis_indices':states.tolist()}
    return H, initial, terms, number_mask, spin_mask, info, snapshot


def run_campaign(args, campaign, manifest):
    init_numerics(campaign/'cache')
    from codex_first_experiment.first_reproducibility_experiment import CASES
    rows, traces = [], {}
    started = time.monotonic()
    old_manifest = json.loads((PRIOR/'manifest.json').read_text())
    for case in args.cases:
        print(f'PREPARE {case}', flush=True)
        if case=='smoke':
            source=REPO/CASES[case]
            expected_sha=sha(source)
            old=None
        else:
            item=old_manifest['cases'][case]
            source=PRIOR/'inputs'/item['relative_path']
            expected_sha=item['sha256']
            old=json.loads((PRIOR/f'{case}.repeat1.json').read_text())
        tensor=campaign/'inputs'/source.name
        tensor.parent.mkdir(exist_ok=True)
        with source.open('rb') as src, tensor.open('xb') as dst:
            shutil.copyfileobj(src,dst)
        if sha(tensor)!=expected_sha or sha(source)!=expected_sha:
            raise RuntimeError('Frozen tensor mismatch')
        H,initial,orders,nmask,smask,info,snapshot=prepare_case(case,tensor,old,args.controls)
        info['tensor_sha256']=expected_sha
        info['tensor_file']=str(tensor.relative_to(campaign))
        save(campaign/f'{case}.operators.json',snapshot)
        save(campaign/f'{case}.setup.json',info)
        row_number=0
        for duration in args.times:
            maximum=max(args.steps)
            if any(maximum % r for r in args.steps):
                raise ValueError('Every step count must divide maximum step count')
            print(f'EXACT {case} T={duration}, dimension={len(initial)}',flush=True)
            exact=expm_multiply(-1j*H,initial,start=0,stop=duration,num=maximum+1,
                                endpoint=True,traceA=-1j*H.diagonal().sum())
            for steps in args.steps:
                trajectory=exact[::maximum//steps]
                selected=list(BASE)
                if args.controls and duration==1 and steps==100:
                    selected=list(orders)
                for name in selected:
                    if time.monotonic()-started > args.max_minutes*60:
                        raise TimeoutError('Campaign exceeded declared runtime budget')
                    begin=time.monotonic()
                    tr=evolve(initial,orders[name],duration,steps)
                    metrics=physical_metrics(trajectory[-1],tr,nmask,smask)
                    row={'case':case,'ordering':name,'T':duration,'steps':steps,'dt':duration/steps,
                         'nominal_exponential_count':steps*info['pauli_terms'],**metrics,
                         **info['leading_defects'][name]}
                    should_decompose=(case in ('f2','smoke') or (duration==1 and steps==100))
                    if should_decompose:
                        details,trace=decompose(trajectory,orders[name],duration,tr)
                        row.update(details)
                        trace_name=f'{case}.{name}.T{duration:g}.r{steps}.trace.json'
                        save(campaign/trace_name,trace)
                        row['trace_file']=trace_name
                    if old and name in BASE and duration==1 and steps==100:
                        reference=next(r for r in old['rows'] if r['ordering']==BASE[name])
                        row['old_one_minus_overlap']=reference['one_minus_overlap']
                        row['old_anchor_match']=math.isclose(metrics['one_minus_overlap'],reference['one_minus_overlap'],
                                                               rel_tol=.005,abs_tol=2e-12)
                        if not row['old_anchor_match']:
                            raise RuntimeError('New implementation failed old physical-error anchor')
                    row['wall_seconds']=time.monotonic()-begin
                    rows.append(row)
                    row_number+=1
                    save(campaign/f'{case}.row{row_number:03d}.json',row)
                    print(f"DONE {case} {name} T={duration:g} r={steps} E={metrics['one_minus_overlap']:.7g} "
                          f"leak={metrics['spin_sector_leakage']:.5g} C={row.get('coherence_factor','-')} "
                          f"{row['wall_seconds']:.1f}s",flush=True)
                    if source_hashes()!=manifest['source_sha256'] or sha(Path(__file__))!=manifest['driver_sha256']:
                        raise RuntimeError('Sources changed during campaign')
                    if sha(tensor)!=expected_sha:
                        raise RuntimeError('Frozen input changed during campaign')
    save(campaign/'results.json',rows)
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+',choices=['f2','nh3','li2','smoke'],default=['f2','nh3','li2'])
    parser.add_argument('--times',type=float,nargs='+',default=[.25,.5,1.])
    parser.add_argument('--steps',type=int,nargs='+',default=[25,50,100,200])
    parser.add_argument('--controls',action='store_true')
    parser.add_argument('--max-minutes',type=float,default=45)
    parser.add_argument('--output-root',type=Path,default=Path(__file__).resolve().parent/'runs')
    args=parser.parse_args()
    if any(v<=0 or not math.isfinite(v) for v in args.times+[args.max_minutes]) or any(r<=0 for r in args.steps):
        parser.error('Positive finite times and positive integer step counts required')
    if len(set(args.cases))!=len(args.cases) or len(set(args.times))!=len(args.times) or len(set(args.steps))!=len(args.steps):
        parser.error('No duplicate cases, times or step counts')
    args.output_root.mkdir(parents=True,exist_ok=True)
    campaign=Path(tempfile.mkdtemp(prefix='mechanism-',dir=args.output_root))
    print(f'CAMPAIGN {campaign}',flush=True)
    manifest={'repo':str(REPO),'commit':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),
              'tracked_worktree_status':git('status','--porcelain','--untracked-files=no'),
              'source_sha256':source_hashes(),'driver_sha256':sha(Path(__file__)),
              'protocol_sha256':sha(Path(__file__).with_name('PROTOCOL.md')),
              'prior_campaign':str(PRIOR),'prior_manifest_sha256':sha(PRIOR/'manifest.json'),
              'settings':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
              'python':sys.version,'versions':{x:importlib.metadata.version(x) for x in ['numpy','scipy','openfermion','numba']}}
    save(campaign/'manifest.json',manifest)
    # Freeze this new driver and protocol; original QHAT is pinned by commit and hashes.
    for name in ('mechanism_experiment.py','PROTOCOL.md','test_mechanism_experiment.py'):
        source=Path(__file__).with_name(name)
        with source.open('rb') as src,(campaign/name).open('xb') as dst:
            shutil.copyfileobj(src,dst)
    try:
        rows=run_campaign(args,campaign,manifest)
        save(campaign/'completion.json',{'status':'complete','rows':len(rows),
              'decomposed_rows':sum('coherence_factor' in r for r in rows)})
    except Exception as exc:
        save(campaign/'failure.json',{'type':type(exc).__name__,'error':str(exc)})
        raise
    print(f'COMPLETE {len(rows)} rows: {campaign}',flush=True)


if __name__=='__main__':
    main()
