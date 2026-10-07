#!/usr/bin/env python3
"""Small fixed-input comparison of fermionic signed and magnitude orders.

Each invocation writes a fresh run. Existing physics code and results are read
only. Uses the same stable metric and leakage-retaining exact closure as the
September mechanism study.
"""
from pathlib import Path
import importlib.metadata
import json
import math
import os
import shutil
import signal
import sys
import tempfile
import time

REPO = Path('/Users/albertlee0125/Repos/qhat')
HELPERS = REPO / 'codex_mechanism_experiment'
PRIOR = HELPERS / 'runs/mechanism-t06qcgac'
sys.path.insert(0, str(HELPERS))
import mechanism_experiment as m
import validate_structure_and_prediction as v

CASES = ['f2', 'nh3', 'li2']
TIMES = [0.25, 1.0]
STEPS = 100
TOL = 1e-12
METHODS = {'fermionic_signed': 'signed_ascending',
           'fermionic_magnitude': 'magnitude_descending'}


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as a, destination.open('xb') as b:
        shutil.copyfileobj(a, b)
    if m.sha(source) != m.sha(destination):
        raise RuntimeError('Copy hash mismatch')


def structure(keys, coefficients, buckets, fallback, n):
    from openfermion import QubitOperator
    counts = []
    for parity in (0, 1):
        op = QubitOperator()
        op.terms = {((q, 'Z'),): -0.5 for q in range(parity, n, 2)}
        counts.append(op)
    noncommuting = 0
    maximum = 0.0
    for _, bucket in buckets:
        for pos, i in enumerate(bucket):
            a = dict(keys[i])
            for j in bucket[pos+1:]:
                b = dict(keys[j])
                noncommuting += sum(a[q] != b[q] for q in a.keys() & b.keys()) % 2
        block = QubitOperator()
        block.terms = {keys[i]: coefficients[keys[i]] for i in bucket}
        maximum = max(maximum, *(v.commutator_bound(block, count) for count in counts))
    return {'nonempty_blocks': len(buckets), 'fallback_count': len(fallback),
            'within_block_noncommuting_pairs': noncommuting,
            'max_block_Nalpha_or_Nbeta_commutator_l1': maximum}


def main():
    runs = Path(__file__).resolve().parent / 'runs'
    runs.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='magnitude-', dir=runs))
    print('RUN', run, flush=True)
    started = time.monotonic()
    def timeout(*_):
        raise TimeoutError('10-minute numerical comparison budget exceeded')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(600)
    source_hashes = m.source_hashes()
    helper_hashes = {name: m.sha(HELPERS/name) for name in
                     ['mechanism_experiment.py', 'validate_structure_and_prediction.py']}
    manifest = {'commit': m.git('rev-parse', 'HEAD'), 'prior_campaign': str(PRIOR),
                'prior_results_sha256': m.sha(PRIOR/'results.json'),
                'source_sha256': source_hashes, 'helper_sha256': helper_hashes,
                'driver_sha256': m.sha(Path(__file__)), 'cases': CASES,
                'times': TIMES, 'steps': STEPS, 'coefficient_tolerance': TOL,
                'methods': METHODS, 'metric': 'stable normalized 1 - abs(overlap)',
                'order_rule': 'recompute first-occurrence ownership for each parent order',
                'python': sys.version,
                'versions': {x: importlib.metadata.version(x) for x in
                             ['numpy', 'scipy', 'openfermion', 'numba']}}
    m.save(run/'manifest.json', manifest)
    copy_file(Path(__file__), run/'source/compare.py')
    for name in helper_hashes:
        copy_file(HELPERS/name, run/'source'/name)
    rows, comparisons, validations = [], [], {}
    try:
        m.init_numerics(run/'cache')
        np = m.np
        from openfermion import get_fermion_operator
        prior_rows = json.loads((PRIOR/'results.json').read_text())
        for case in CASES:
            print('PREPARE', case, flush=True)
            info = json.loads((PRIOR/f'{case}.setup.json').read_text())
            tensor_source = PRIOR/info['tensor_file']
            if m.sha(tensor_source) != info['tensor_sha256']:
                raise RuntimeError('Frozen input hash mismatch')
            copy_file(tensor_source, run/info['tensor_file'])
            for suffix in ['setup.json', 'operators.json']:
                copy_file(PRIOR/f'{case}.{suffix}', run/f'{case}.{suffix}')
            snapshot, info = v.load_snapshot(run, case)
            n, electrons = info['n_qubits'], info['n_electrons']
            keys, coefficients, states, reduced = v.restore(snapshot, n, electrons)
            H, parity = m.assemble_hamiltonian(reduced, len(states))
            initial = np.zeros(len(states), complex)
            initial[0] = 1
            nmask = np.array([int(x).bit_count() == electrons for x in states])
            spin_states = set(map(int, m.matrixfree.spin_sector_basis_indices(n, electrons)))
            smask = np.array([int(x) in spin_states for x in states])
            interaction, _ = m.baseline.load_interaction_operator(run/info['tensor_file'])
            fermion = m.baseline.clean_fermion_operator(get_fermion_operator(interaction), TOL)
            parents = m.baseline.build_hermitian_fermion_terms(fermion, TOL)
            mapping = m.baseline.precompute_fermion_to_pauli_indices(
                parents, coefficients, {k:i for i,k in enumerate(keys)}, TOL)
            all_orders = m.baseline.build_deterministic_orderings(fermion, keys, coefficients, n, TOL)
            original = m.matrixfree.compile_ordered_terms(keys, coefficients, n, TOL)
            full_initial = m.matrixfree.build_hartree_fock_state(n, electrons)
            terms, details, owners = {}, {}, {}
            for name, method in METHODS.items():
                parent_order = m.baseline.fermionic_term_order_indices(parents, method, TOL)
                weights = [m.baseline.fermionic_term_signed_weight(parents[i], TOL) for i in parent_order]
                ordering_weights = weights if method == 'signed_ascending' else [-abs(x) for x in weights]
                if any(a > b for a,b in zip(ordering_weights, ordering_weights[1:])):
                    raise RuntimeError('Coefficient sorting validation failed')
                order = m.baseline.induced_pauli_order_indices(parent_order, mapping, len(keys))
                expected_name = ('fermionic_signed_coefficient_lexicographic' if name == 'fermionic_signed'
                                 else 'fermionic_magnitude_descending_lexicographic')
                if [keys[i] for i in order] != all_orders[expected_name]:
                    raise RuntimeError('Mismatch with existing ordering implementation')
                if sorted(order) != list(range(len(keys))):
                    raise RuntimeError('Pauli terms are not a complete unchanged permutation')
                if name == 'fermionic_signed' and order != snapshot['orders'][name]:
                    raise RuntimeError('Signed order changed from frozen campaign')
                buckets, owner, fallback = m.ablation.build_reference_parent_buckets(parent_order, mapping, len(keys))
                if m.ablation.flatten_reference_buckets(buckets, fallback) != order:
                    raise RuntimeError('Parent ownership disagrees with induced Pauli order')
                owners[name] = owner
                terms[name] = [reduced[i] for i in order]
                full = m.evolve(full_initial, [original[i] for i in order], .01)
                small = m.evolve(initial, terms[name], .01)
                embedded = np.zeros_like(full)
                embedded[states] = small
                residual = float(np.linalg.norm(full-embedded))
                if residual > 1e-10:
                    raise RuntimeError('Full-vs-closed-space one-step validation failed')
                lead = m.leading_defect(terms[name], initial, H, parity)
                details[name] = {
                    'first_12_parent_coefficients': weights[:12],
                    'parent_order': parent_order, 'pauli_order': order,
                    'physical_order_sha256': m.obj_sha([[keys[i], complex(coefficients[keys[i]]).real.hex(),
                                                        complex(coefficients[keys[i]]).imag.hex()] for i in order]),
                    'full_space_one_step_residual': residual,
                    'bch_hf_norm': 2*float(np.linalg.norm(lead)),
                    **structure(keys, coefficients, buckets, fallback, n)}
                print('ORDER', case, name, 'BCH', details[name]['bch_hf_norm'],
                      'commutator', details[name]['max_block_Nalpha_or_Nbeta_commutator_l1'], flush=True)
            details['pauli_terms_changing_first_parent_owner'] = int(np.count_nonzero(
                owners['fermionic_signed'] != owners['fermionic_magnitude']))
            details['same_pauli_sequence'] = details['fermionic_signed']['pauli_order'] == details['fermionic_magnitude']['pauli_order']
            m.save(run/f'{case}.comparison_setup.json', details)
            validations[case] = {'n_qubits':n, 'n_electrons':electrons,
                                'tensor_sha256':info['tensor_sha256'],
                                'identity_free_jw_sha256':info['fingerprints']['identity_free_jw_sha256'],
                                'orders':details}
            for T in TIMES:
                exact = m.expm_multiply((-1j*T)*H, initial, traceA=(-1j*T)*H.diagonal().sum())
                pair = {}
                for name in METHODS:
                    before = time.monotonic()
                    approx = m.evolve(initial, terms[name], T, STEPS)
                    metric = m.physical_metrics(exact, approx, nmask, smask)
                    row = {'case':case, 'ordering':name, 'T':T, 'steps':STEPS,
                           'bch_hf_norm':details[name]['bch_hf_norm'], **metric,
                           'seconds':time.monotonic()-before}
                    if name == 'fermionic_signed':
                        old = next(x for x in prior_rows if x['case']==case and x['ordering']==name
                                   and x['T']==T and x['steps']==STEPS)
                        row['prior_signed_E'] = old['one_minus_overlap']
                        row['signed_anchor_match'] = math.isclose(metric['one_minus_overlap'],
                            old['one_minus_overlap'], rel_tol=1e-6, abs_tol=1e-16)
                        if not row['signed_anchor_match']:
                            raise RuntimeError('Current signed error fails existing stable-metric anchor')
                    rows.append(row)
                    pair[name] = row
                    m.save(run/f'{case}.{name}.T{T:g}.json', row)
                    print(f"DONE {case} {name} T={T:g} E={metric['one_minus_overlap']:.9g} "
                          f"leak={metric['spin_sector_leakage']:.3g} seconds={row['seconds']:.1f}", flush=True)
                a,b = pair['fermionic_signed'], pair['fermionic_magnitude']
                comparisons.append({'case':case,'T':T,'steps':STEPS,
                    'signed_E':a['one_minus_overlap'],'magnitude_E':b['one_minus_overlap'],
                    'magnitude_over_signed_E':b['one_minus_overlap']/a['one_minus_overlap'],
                    'signed_spin_leakage':a['spin_sector_leakage'],
                    'magnitude_spin_leakage':b['spin_sector_leakage'],
                    'above_floor':not(a['below_analysis_floor'] or b['below_analysis_floor'])})
            if m.sha(run/info['tensor_file']) != info['tensor_sha256']:
                raise RuntimeError('Input changed during calculation')
        if m.source_hashes() != source_hashes or any(m.sha(HELPERS/x)!=h for x,h in helper_hashes.items()):
            raise RuntimeError('Source files changed during calculation')
        m.save(run/'results.json', rows)
        m.save(run/'comparisons.json', comparisons)
        m.save(run/'validation.json', validations)
        m.save(run/'completion.json', {'status':'complete','rows':len(rows),
            'paired_comparisons':len(comparisons),'seconds':time.monotonic()-started,
            'signed_anchors_passed':sum(r.get('signed_anchor_match',False) for r in rows)})
        lines = ['Fermionic parent order: signed ascending versus magnitude descending',
                 'Frozen HGBS-5 inputs, HF, first-order Trotter, r=100, coefficient cutoff=1e-12.',
                 'E = 1 - abs(normalized overlap), evaluated using a stable orthogonal residual.',
                 'Ratio below 1 favors magnitude descending. First-occurrence ownership is recomputed.',
                 '', 'case T signed_E magnitude_E magnitude/signed']
        for c in comparisons:
            lines.append(f"{c['case']} {c['T']:g} {c['signed_E']:.9e} {c['magnitude_E']:.9e} {c['magnitude_over_signed_E']:.6g}")
        lines += ['', 'All six signed anchors agree with the frozen mechanism run.',
                  'Both methods contain exactly the same final merged Pauli terms and coefficients.',
                  'Full-space one-step checks retain leakage and agree with the exact closed-space calculation.',
                  'Three fixed molecules and two times are a small case study, not a general ranking claim.']
        with (run/'SUMMARY.txt').open('x') as f:
            f.write('\n'.join(lines)+'\n')
        print('COMPLETE',run,flush=True)
    except Exception as exc:
        m.save(run/'failure.json',{'type':type(exc).__name__,'message':str(exc),
                                 'completed_rows':len(rows)})
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
