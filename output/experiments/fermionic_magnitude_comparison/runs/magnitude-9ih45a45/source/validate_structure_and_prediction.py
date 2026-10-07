#!/usr/bin/env python3
"""Exploratory structural proof checks and unfitted finite-time BCH prediction."""
import argparse
import json
import math
import shutil
from pathlib import Path
import tempfile
import mechanism_experiment as m


def restore(snapshot, n, electrons):
    coeff={tuple(tuple(v) for v in k):complex(float.fromhex(re),float.fromhex(im))
           for k,re,im in snapshot['coefficients']}
    # Snapshot orders index the original (not sorted) Pauli list. Reconstruct
    # its order from the first experiment's preserved order key list, not from
    # canonical coefficient sorting. The runner needs to record this mapping.
    raw_keys=snapshot.get('raw_pauli_keys')
    if raw_keys is None:
        # The immutable Hamiltonian is reloaded from the frozen tensor below;
        # this fallback is deliberately not an assumption about sorted order.
        raise ValueError('raw_pauli_keys required')
    keys=[tuple(tuple(v) for v in key) for key in raw_keys]
    original=m.matrixfree.compile_ordered_terms(keys,coeff,n,1e-12)
    hf=int(m.np.argmax(abs(m.matrixfree.build_hartree_fock_state(n,electrons))))
    states,reduced,_=m.affine_space(original,hf)
    m.np.testing.assert_array_equal(states,snapshot['closed_basis_indices'])
    return keys,coeff,states,reduced


def load_snapshot(campaign,case):
    from openfermion import get_fermion_operator,jordan_wigner
    snap=json.loads((campaign/f'{case}.operators.json').read_text())
    info=json.loads((campaign/f'{case}.setup.json').read_text())
    interaction,n=m.baseline.load_interaction_operator(campaign/info['tensor_file'])
    fermion=m.baseline.clean_fermion_operator(get_fermion_operator(interaction),1e-12)
    jw=jordan_wigner(fermion);jw.compress(abs_tol=1e-12)
    coeff={k:v for k,v in jw.terms.items() if k and abs(v)>1e-12}
    if m.obj_sha(m.canonical_terms(coeff))!=info['fingerprints']['identity_free_jw_sha256']:
        raise RuntimeError('Restored Hamiltonian fingerprint mismatch')
    snap['raw_pauli_keys']=list(coeff)
    for name,order in snap['orders'].items():
        keys=list(coeff)
        h=m.obj_sha([[keys[i],complex(coeff[keys[i]]).real.hex(),complex(coeff[keys[i]]).imag.hex()] for i in order])
        if h!=info['fingerprints']['physical_order_sha256'][name]:
            raise RuntimeError('Restored raw index mapping fails physical order hash')
    return snap,info


def commutator_bound(a,b):
    # SymbolicOperator +/- implicitly deletes small coefficients. Form the
    # two unpruned products and subtract their dictionaries ourselves.
    ab=(a*b).terms;ba=(b*a).terms
    return sum(abs(ab.get(k,0)-ba.get(k,0)) for k in ab.keys() | ba.keys())


def structure(campaign):
    from openfermion import QubitOperator
    output={}
    for case in ('f2','nh3','li2'):
        snapshot,info=load_snapshot(campaign,case)
        keys,coeff,_,_=restore(snapshot,info['n_qubits'],info['n_electrons'])
        counts=[]
        for parity in (0,1):
            N=QubitOperator()
            # Constants commute, so omit them without changing commutators.
            N.terms={((q,'Z'),):-.5 for q in range(parity,info['n_qubits'],2)}
            counts.append(N)
        total=QubitOperator();total.terms={**counts[0].terms,**counts[1].terms}
        bounds=[];noncommuting=0
        for _,bucket in snapshot['parent_buckets']:
            for pos,i in enumerate(bucket):
                a=dict(keys[i])
                for j in bucket[pos+1:]:
                    b=dict(keys[j])
                    noncommuting+=sum(a[q]!=b[q] for q in a.keys() & b.keys())%2
            H=QubitOperator();H.terms={keys[i]:coeff[keys[i]] for i in bucket}
            bounds.append([commutator_bound(H,N) for N in (counts[0],counts[1],total)])
        output[case]={'parent_blocks':len(bounds),'within_bucket_noncommuting_pairs':noncommuting,
                      'max_bucket_commutator_coefficient_l1_Nalpha':max(b[0] for b in bounds),
                      'max_bucket_commutator_coefficient_l1_Nbeta':max(b[1] for b in bounds),
                      'max_bucket_commutator_coefficient_l1_N':max(b[2] for b in bounds),
                      'fallback_count':len(snapshot['fallback'])}
    return output


def defect_columns(terms,sector,H,parity):
    np=m.np
    dim=H.shape[0];width=len(sector)
    E=np.zeros((dim,width),complex);E[sector,np.arange(width)]=1
    first=np.zeros_like(E);second=np.zeros_like(E)
    coords=np.arange(dim,dtype=np.int64)
    for c,x,z,y in terms:
        src=coords ^ int(x)
        phase=(1j**int(y))*(1-2*parity[src & int(z)])
        second-=1j*c*phase[:,None]*first[src]
        second[sector,np.arange(width)]-=.5*c*c
        destination=sector ^ int(x)
        initial_phase=(1j**int(y))*(1-2*parity[sector & int(z)])
        first[destination,np.arange(width)]-=1j*c*initial_phase
    return second+.5*(H@(H@E))


def predictions(campaign,destination):
    np=m.np
    snapshot,info=load_snapshot(campaign,'f2')
    _,_,states,reduced=restore(snapshot,info['n_qubits'],info['n_electrons'])
    H,parity=m.assemble_hamiltonian(reduced,len(states))
    n=info['n_qubits'];e=info['n_electrons']
    sector=np.array([i for i,x in enumerate(states) if
                     sum((int(x)>>k)&1 for k in range(0,n,2))==e//2 and
                     sum((int(x)>>k)&1 for k in range(1,n,2))==e//2],dtype=np.int64)
    Hs=H[sector][:,sector]
    number_mask=np.array([int(x).bit_count()==e for x in states])
    spin_mask=np.zeros(len(states),bool);spin_mask[sector]=True
    initial=np.zeros(len(states),complex);initial[0]=1
    s0=np.zeros(len(sector),complex);s0[np.flatnonzero(sector==0)[0]]=1
    predicted=[];finite_checks={};new_states={};terms_by_name={}
    for name in m.BASE:
        terms=[reduced[i] for i in snapshot['orders'][name]];terms_by_name[name]=terms
        D=defect_columns(terms,sector,H,parity)
        direct=m.leading_defect(terms,initial,H,parity)
        if np.linalg.norm(D@s0-direct)>1e-10:raise RuntimeError('Defect column check failed')
        rng=np.random.default_rng(1709)
        probe_sector=rng.normal(size=len(sector))+1j*rng.normal(size=len(sector))
        probe_sector/=np.linalg.norm(probe_sector)
        probe=np.zeros(len(states),complex);probe[sector]=probe_sector
        differences=[]
        for dt in (.002,.001):
            exact=m.expm_multiply(-1j*dt*H,probe,traceA=-1j*dt*H.diagonal().sum())
            finite=(m.evolve(probe,terms,dt)-exact)/dt**2
            differences.append(float(np.linalg.norm(finite-D@probe_sector)))
        finite_checks[name]={'dt_values':[.002,.001],'residuals':differences,
                             'ratio':differences[1]/differences[0]}
        if differences[1]>.65*differences[0]:raise RuntimeError('Defect finite-difference convergence failed')
        augmented=m.sparse.bmat([[-1j*H,m.sparse.csr_matrix(D)],[None,-1j*Hs]],format='csr')
        vector=np.concatenate([np.zeros(len(states),complex),s0])
        trajectory=m.expm_multiply(augmented,vector,start=0,stop=1,num=5,
                                   traceA=augmented.diagonal().sum())
        for index,T in enumerate((0,.25,.5,.75,1.)):
            if not index:continue
            eta=trajectory[index,:len(states)]
            exact=np.zeros(len(states),complex);exact[sector]=trajectory[index,len(states):]
            phi=exact/np.linalg.norm(exact)
            qeta=eta-phi*np.vdot(phi,eta)
            for r in ([100] if T==.75 else [25,50,100,200]):
                predicted.append({'case':'f2','ordering':name,'T':T,'steps':r,
                                  'predicted_infidelity':(T/r)**2*m.norm2(qeta)})
            if T==.75:new_states[name]=exact
    # Immutable prediction file is written BEFORE new T=0.75 product formulas.
    m.save(destination/'predictions_before_new_time.json',predicted)
    m.save(destination/'defect_finite_difference_checks.json',finite_checks)
    new_rows=[]
    for name,exact in new_states.items():
        approximate=m.evolve(initial,terms_by_name[name],.75,100)
        row={'case':'f2','ordering':name,'T':.75,'steps':100,
             **m.physical_metrics(exact,approximate,number_mask,spin_mask)}
        new_rows.append(row)
    m.save(destination/'new_time_results.json',new_rows)
    previous=json.loads((campaign/'results.json').read_text())
    comparisons=[]
    for prediction in predicted:
        actual=next(x for x in previous+new_rows if all(x[k]==prediction[k] for k in ('case','ordering','T','steps')))
        comparisons.append({**prediction,'actual_infidelity':actual['infidelity'],
                            'relative_difference':prediction['predicted_infidelity']/actual['infidelity']-1})
    m.save(destination/'prediction_comparisons.json',comparisons)
    return comparisons


def main():
    parser=argparse.ArgumentParser();parser.add_argument('campaign',type=Path);args=parser.parse_args()
    destination=Path(tempfile.mkdtemp(prefix='exploratory-',dir=args.campaign))
    print(f'EXPLORATORY {destination}',flush=True)
    for name in ('validate_structure_and_prediction.py','test_prediction.py','EXPLORATORY_ADDENDUM.md'):
        with Path(__file__).with_name(name).open('rb') as src,(destination/name).open('xb') as dst:
            shutil.copyfileobj(src,dst)
    m.save(destination/'manifest.json',{'main_campaign':str(args.campaign),
           'main_manifest_sha256':m.sha(args.campaign/'manifest.json'),
           'repository_commit':m.git('rev-parse','HEAD'),
           'source_sha256':m.sha(Path(__file__))})
    m.init_numerics(destination/'cache')
    structures=structure(args.campaign);m.save(destination/'structure_checks.json',structures)
    print('STRUCTURE',structures,flush=True)
    comparison=predictions(args.campaign,destination)
    summary={'prediction_rows':len(comparison),'new_time_rows':3,
             'max_relative_prediction_difference':max(abs(x['relative_difference']) for x in comparison),
             'r100_max_relative_prediction_difference':max(abs(x['relative_difference']) for x in comparison if x['steps']==100),
             'r200_max_relative_prediction_difference':max(abs(x['relative_difference']) for x in comparison if x['steps']==200),
             'source_sha256':m.sha(Path(__file__)),'addendum_sha256':m.sha(Path(__file__).with_name('EXPLORATORY_ADDENDUM.md'))}
    m.save(destination/'summary.json',summary)
    print('COMPLETE',summary,flush=True)


if __name__=='__main__':main()
