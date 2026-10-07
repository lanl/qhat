#!/usr/bin/env python3
"""Prospective transfer test; read PROTOCOL.md before interpreting results."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
HELPERS = REPO / 'codex_mechanism_experiment'
sys.path.insert(0, str(HELPERS))
import mechanism_experiment as m
import validate_structure_and_prediction as v

TIMES = (.25, .75, 1.)
STEPS = (50, 100, 200)
SCALES = (1., 1.25)


def specifications():
    a = 1.48 / math.sqrt(3)
    x = 1.34 * math.sin(math.radians(46))
    z = 1.34 * math.cos(math.radians(46))
    geometries = {
        'HCN': [('C', (0.,0.,0.)), ('N', (0.,0.,1.16)), ('H', (0.,0.,-1.06))],
        'CH2O': [('C', (0.,0.,0.)), ('O', (0.,0.,1.21)), ('H', (.94,0.,-.59)), ('H', (-.94,0.,-.59))],
        'H2S': [('S', (0.,0.,0.)), ('H', (x,0.,z)), ('H', (-x,0.,z))],
        'SiH4': [('Si', (0.,0.,0.)), ('H', (a,a,a)), ('H', (a,-a,-a)), ('H', (-a,a,-a)), ('H', (-a,-a,a))],
    }
    return [{'case': f'{name}_s{scale:.2f}', 'family': name, 'scale':scale,
             'geometry':[(atom,[scale*q for q in xyz]) for atom,xyz in geometry],
             'basis':'6-31g', 'charge':0, 'spin':0,
             'active_electrons':4, 'active_spatial_orbitals':5}
            for name,geometry in geometries.items() for scale in SCALES]


def tensor_relative(spec):
    return Path('inputs')/spec['family']/f"s-{spec['scale']:.2f}"/spec['basis']/f"{spec['case']}_6-31g_as-004-006.tensors.npz"


def generate(spec, destination):
    from pyscf import gto, scf, lib
    from openfermion import MolecularData, get_sparse_operator, jordan_wigner, get_fermion_operator
    from openfermionpyscf._run_pyscf import compute_integrals
    np=m.np
    lib.num_threads(1)
    mol=gto.M(atom=spec['geometry'], basis=spec['basis'], charge=0, spin=0,
              unit='Angstrom', symmetry=False, verbose=0)
    mf=scf.RHF(mol);mf.conv_tol=1e-12;mf.max_cycle=200
    mf.chkfile=None
    energy=float(mf.kernel())
    if not mf.converged:raise RuntimeError(f"SCF did not converge: {spec['case']}")
    # A deterministic sign convention; arbitrary degenerate rotations are not fixed.
    C=mf.mo_coeff.copy()
    for j in range(C.shape[1]):
        if C[np.argmax(abs(C[:,j])),j]<0:C[:,j]*=-1
    mf.mo_coeff=C
    one,two=compute_integrals(mol,mf)
    data=MolecularData(spec['geometry'],spec['basis'],1,0,
                       filename=str(destination/'unused_molecular_data'))
    data.n_orbitals=C.shape[1];data.n_qubits=2*C.shape[1]
    data.nuclear_repulsion=float(mol.energy_nuc())
    data.one_body_integrals=one;data.two_body_integrals=two
    nocc=mol.nelectron//2
    frozen=list(range(nocc-2));active=list(range(nocc-2,nocc+3))
    if active[-1]>=data.n_orbitals:raise RuntimeError('Insufficient virtual orbitals')
    interaction=data.get_molecular_hamiltonian(occupied_indices=frozen,active_indices=active)
    hf=m.matrixfree.build_hartree_fock_state(10,4)
    H=get_sparse_operator(jordan_wigner(get_fermion_operator(interaction)),n_qubits=10)
    residual=abs(float(m.np.vdot(hf,H@hf).real)-energy)
    if residual>1e-8:raise RuntimeError(f'HF tensor conversion mismatch {residual}')
    filename=destination/tensor_relative(spec);filename.parent.mkdir(parents=True,exist_ok=True)
    with filename.open('xb') as f:
        np.savez_compressed(f,constant=interaction.constant,one_body=interaction.one_body_tensor,
                           two_body=interaction.two_body_tensor)
    with filename.with_suffix('.orbitals.npz').open('xb') as f:
        np.savez_compressed(f,mo_coeff=C,mo_energy=mf.mo_energy,mo_occ=mf.mo_occ,
                           overlap=mf.get_ovlp(),active_spatial=active,frozen_spatial=frozen)
    metadata={**spec,'unit':'Angstrom','hf_energy_hartree':energy,'scf_converged':True,
              'total_electrons':mol.nelectron,'total_spatial_orbitals':data.n_orbitals,
              'active_spatial_indices':active,'frozen_occupied_spatial_indices':frozen,
              'hf_tensor_expectation_residual':residual,'tensor_sha256':m.sha(filename)}
    m.save(filename.with_suffix('.metadata.json'),metadata)
    return filename,metadata


def structure(snapshot):
    from openfermion import QubitOperator
    keys=[tuple(tuple(v) for v in k) for k in snapshot['raw_pauli_keys']]
    coeff={tuple(tuple(v) for v in k):complex(float.fromhex(re),float.fromhex(im))
           for k,re,im in snapshot['coefficients']}
    counts=[]
    for spin in (0,1):
        op=QubitOperator();op.terms={((q,'Z'),):-.5 for q in range(spin,10,2)};counts.append(op)
    total=QubitOperator();total.terms={**counts[0].terms,**counts[1].terms}
    bounds=[];noncommuting=0
    for _,bucket in snapshot['parent_buckets']:
        for pos,i in enumerate(bucket):
            a=dict(keys[i])
            for j in bucket[pos+1:]:
                b=dict(keys[j]);noncommuting+=sum(a[q]!=b[q] for q in a.keys()&b.keys())%2
        block=QubitOperator();block.terms={keys[i]:coeff[keys[i]] for i in bucket}
        bounds.append([v.commutator_bound(block,N) for N in (*counts,total)])
    return {'noncommuting_pairs':noncommuting,'parent_blocks':len(bounds),
            'fallback_count':len(snapshot['fallback']),
            'max_Nalpha_Nbeta_N_commutator_l1':[max(b[j] for b in bounds) for j in range(3)]}


def predict(case,H,initial,orders,smask,info):
    np=m.np
    sector=np.flatnonzero(smask);Hs=H[sector][:,sector]
    s0=initial[sector];dimension=len(initial)
    _,parity=m.assemble_hamiltonian(orders['fermionic_signed'],dimension)
    records=[];checks={}
    rng=np.random.default_rng(1709)
    probe_sector=rng.normal(size=len(sector))+1j*rng.normal(size=len(sector))
    probe_sector/=np.linalg.norm(probe_sector)
    probe=np.zeros(dimension,complex);probe[sector]=probe_sector
    for name in m.BASE:
        D=v.defect_columns(orders[name],sector,H,parity)
        column_check=float(np.linalg.norm(D@s0-m.leading_defect(orders[name],initial,H,parity)))
        if column_check>1e-10:raise RuntimeError('Leading-defect column check failed')
        differences=[]
        for dt in (.002,.001):
            exact=m.expm_multiply(-1j*dt*H,probe,traceA=-1j*dt*H.diagonal().sum())
            finite=(m.evolve(probe,orders[name],dt)-exact)/dt**2
            differences.append(float(np.linalg.norm(finite-D@probe_sector)))
        checks[name]={'column_check':column_check,'finite_difference_residuals':differences,
                      'finite_difference_ratio':differences[1]/differences[0]}
        if differences[1]>.65*differences[0]:raise RuntimeError('Local derivative convergence failed')
        A=m.sparse.bmat([[-1j*H,m.sparse.csr_matrix(D)],[None,-1j*Hs]],format='csr')
        vector=np.concatenate([np.zeros(dimension,complex),s0])
        trajectory=m.expm_multiply(A,vector,start=0,stop=1,num=5,traceA=A.diagonal().sum())
        for T in TIMES:
            eta=trajectory[round(T*4),:dimension]
            exact=np.zeros(dimension,complex);exact[sector]=trajectory[round(T*4),dimension:]
            phi=exact/np.linalg.norm(exact);qeta=eta-phi*np.vdot(phi,eta)
            for r in STEPS:
                records.append({'case':case,'ordering':name,'T':T,'steps':r,
                    'predicted_infidelity':(T/r)**2*m.norm2(qeta),
                    'predicted_inside_sector_error':(T/r)**2*m.norm2(qeta[smask]),
                    'predicted_outside_sector_error':(T/r)**2*m.norm2(qeta[~smask]),
                    'static_bch_proxy_infidelity':T**4/r**2*(info['leading_defects'][name]['projected_bch_hf_norm']/2)**2})
    return records,checks


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs-from',type=Path,help='Reuse frozen inputs from a completed campaign')
    args=parser.parse_args()
    (ROOT/'runs').mkdir(exist_ok=True)
    destination=Path(tempfile.mkdtemp(prefix='transfer-',dir=ROOT/'runs'))
    print(f'CAMPAIGN {destination}',flush=True)
    sources=[Path(__file__),ROOT/'PROTOCOL.md',ROOT/'test_external_validation.py',
             HELPERS/'mechanism_experiment.py',HELPERS/'validate_structure_and_prediction.py']
    (destination/'source').mkdir()
    for p in sources:shutil.copy2(p,destination/'source'/p.name)
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'specifications':specifications(),
              'repository_commit':m.git('rev-parse','HEAD'),'repository_status':m.git('status','--short'),
              'source_hashes':{str(p):m.sha(p) for p in sources},'tracked_source_hashes':m.source_hashes(),
              'environment':{p:importlib.metadata.version(p) for p in ('numpy','scipy','pyscf','openfermion','openfermionpyscf','numba')},
              'inputs_from':str(args.inputs_from) if args.inputs_from else None,
              'expected_baseline_rows':216,'expected_control_rows':56}
    m.save(destination/'manifest.json',manifest)
    m.init_numerics(destination/'cache')
    os.environ['TMPDIR']=str(destination/'cache')
    tempfile.tempdir=str(destination/'cache')
    prepared={};all_predictions=[];structures={};fd_checks={};started=time.time()
    for spec in specifications():
        case=spec['case']
        if args.inputs_from:
            source=args.inputs_from/tensor_relative(spec)
            tensor=destination/tensor_relative(spec);tensor.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,tensor)
            metadata=json.loads(source.with_suffix('.metadata.json').read_text())
            if m.sha(tensor)!=metadata['tensor_sha256']:raise RuntimeError('Frozen input hash mismatch')
            shutil.copy2(source.with_suffix('.metadata.json'),tensor.with_suffix('.metadata.json'))
            shutil.copy2(source.with_suffix('.orbitals.npz'),tensor.with_suffix('.orbitals.npz'))
        else:tensor,metadata=generate(spec,destination)
        print('GENERATED',case,metadata['hf_tensor_expectation_residual'],flush=True)
        H,initial,orders,nmask,smask,info,snapshot=m.prepare_case(case,tensor,None,True)
        # Preserve the raw list whose indices the ordering permutations use.
        from openfermion import get_fermion_operator,jordan_wigner
        interaction,_=m.baseline.load_interaction_operator(tensor)
        fermion=m.baseline.clean_fermion_operator(get_fermion_operator(interaction),1e-12)
        jw=jordan_wigner(fermion);jw.compress(abs_tol=1e-12)
        snapshot['raw_pauli_keys']=[k for k,c in jw.terms.items() if k and abs(c)>1e-12]
        info.update({'tensor_file':str(tensor.relative_to(destination)),'tensor_sha256':m.sha(tensor),
                     'family':spec['family'],'scale':spec['scale'],'spin_sector_dimension':int(smask.sum())})
        m.save(destination/f'{case}.setup.json',info);m.save(destination/f'{case}.operators.json',snapshot)
        structures[case]=structure(snapshot)
        prediction,fd=predict(case,H,initial,orders,smask,info)
        all_predictions.extend(prediction);fd_checks[case]=fd
        prepared[case]=(H,initial,orders,nmask,smask,info)
        print('PREDICTED',case,'dimension',len(initial),'sector',int(smask.sum()),'terms',info['pauli_terms'],flush=True)
    m.save(destination/'structure_checks.json',structures)
    m.save(destination/'finite_difference_checks.json',fd_checks)
    prediction_file=destination/'predictions_before_target_evolution.json'
    m.save(prediction_file,all_predictions)
    prediction_hash=m.sha(prediction_file)
    m.save(destination/'prediction_freeze.json',{'sha256':prediction_hash,'rows':len(all_predictions),
             'frozen_utc':datetime.now(timezone.utc).isoformat(),'target_evolutions_completed':0})
    print('ALL_PREDICTIONS_FROZEN',prediction_hash,flush=True)
    rows=[]
    for case,(H,initial,orders,nmask,smask,info) in prepared.items():
        for T in TIMES:
            trajectory=m.expm_multiply(-1j*H,initial,start=0,stop=T,num=101,traceA=-1j*H.diagonal().sum())
            exact=trajectory[-1]
            for r in STEPS:
                names=list(m.BASE) if (T,r)!=(1.,100) else list(orders)
                for name in names:
                    approx=m.evolve(initial,orders[name],T,r)
                    row={'case':case,'family':info['family'],'scale':info['scale'],
                         'ordering':name,'T':T,'steps':r,
                         **m.physical_metrics(exact,approx,nmask,smask)}
                    if r==100:
                        decomposition,trace=m.decompose(trajectory,orders[name],T,approx)
                        row.update(decomposition)
                        m.save(destination/f'{case}.{len(rows):03d}.trace.json',trace)
                    m.save(destination/f'{case}.{len(rows):03d}.row.json',row);rows.append(row)
            print('MEASURED',case,'T',T,'rows',len(rows),flush=True)
        if m.sha(destination/info['tensor_file'])!=info['tensor_sha256']:raise RuntimeError('Input changed')
    if len(rows)!=272:raise RuntimeError('Wrong target row count')
    if m.sha(prediction_file)!=prediction_hash:raise RuntimeError('Predictions changed')
    for p,expected in manifest['source_hashes'].items():
        if m.sha(p)!=expected:raise RuntimeError('Source changed during run')
    if m.source_hashes()!=manifest['tracked_source_hashes']:raise RuntimeError('Tracked source changed')
    m.save(destination/'results.json',rows)
    m.save(destination/'completion.json',{'status':'complete','rows':len(rows),'baseline_rows':216,
           'control_rows':56,'seconds':time.time()-started,'prediction_sha256':prediction_hash,
           'completed_utc':datetime.now(timezone.utc).isoformat()})
    print('COMPLETE',destination,'seconds',time.time()-started,flush=True)


if __name__=='__main__':main()
