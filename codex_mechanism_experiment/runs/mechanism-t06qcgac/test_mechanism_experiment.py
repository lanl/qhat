import os
from pathlib import Path
import tempfile
import unittest
import mechanism_experiment as m

_cache=tempfile.TemporaryDirectory(prefix='qhat-mechanism-tests-')
m.init_numerics(Path(_cache.name))
np=m.np
from scipy.linalg import expm


class MechanismTests(unittest.TestCase):
    def setUp(self):
        # Real molecular-style even-Y operators, noncommuting and leaking sector.
        self.keys=[((0,'Z'),),((1,'Z'),),((0,'X'),(1,'X')),((0,'Y'),(1,'Y'))]
        self.coeff=dict(zip(self.keys,[.7,-.23,.13,-.09]))
        self.terms=m.matrixfree.compile_ordered_terms(self.keys,self.coeff,2,1e-12)
        self.states,self.reduced,_=m.affine_space(self.terms,2)
        self.H,self.parity=m.assemble_hamiltonian(self.reduced,len(self.states))
        self.psi=np.zeros(len(self.states),complex);self.psi[0]=1

    def test_affine_closed_and_actions(self):
        full=np.zeros(4,complex);full[self.states]=[.3+.4j,np.sqrt(.75)]
        for original,reduced in zip(self.terms,self.reduced):
            a=m.evolve(full,[original],.37)
            b=m.evolve(full[self.states],[reduced],.37)
            embedded=np.zeros(4,complex);embedded[self.states]=b
            np.testing.assert_allclose(a,embedded,atol=1e-14)

    def test_all_pauli_y_phases(self):
        keys=[((0,p),(1,q)) for p in 'XYZ' for q in 'XYZ']
        coeff={k:.11 for k in keys}
        original=m.matrixfree.compile_ordered_terms(keys,coeff,2,1e-12)
        states,reduced,_=m.affine_space(original,1)
        full=np.array([.2+.3j,.4,.1j,-.2],complex);full/=np.linalg.norm(full)
        for a,b in zip(original,reduced):
            np.testing.assert_allclose(m.evolve(full,[a],.31)[states],
                                       m.evolve(full[states],[b],.31),atol=1e-14)

    def test_hamiltonian_against_openfermion(self):
        from openfermion import QubitOperator,get_sparse_operator
        q=sum((QubitOperator(k,v) for k,v in self.coeff.items()),QubitOperator())
        full=get_sparse_operator(q,n_qubits=2).toarray()
        np.testing.assert_allclose(self.H.toarray(),full[np.ix_(self.states,self.states)],atol=1e-14)

    def test_product_against_dense(self):
        dt=.14
        dense=np.eye(len(self.psi),dtype=complex)
        for term in self.reduced:
            p=np.column_stack([m.apply_p(v,term,self.parity) for v in np.eye(len(self.psi))])
            dense=expm(-1j*dt*term[0]*p)@dense
        np.testing.assert_allclose(m.evolve(self.psi,self.reduced,dt),dense@self.psi,atol=1e-14)

    def test_leading_defect(self):
        v=m.leading_defect(self.reduced,self.psi,self.H,self.parity)
        errors=[]
        for dt in [.002,.001]:
            finite=(m.evolve(self.psi,self.reduced,dt)-expm(-1j*dt*self.H.toarray())@self.psi)/dt**2
            errors.append(np.linalg.norm(finite-v))
        self.assertLess(errors[1],.55*errors[0])

    def test_projected_decomposition_dense_cross_terms(self):
        T=.6;r=7;dt=T/r
        U=expm(-1j*dt*self.H.toarray())
        S=np.column_stack([m.evolve(v,self.reduced,dt) for v in np.eye(len(self.psi),dtype=complex)])
        trajectory=np.array([np.linalg.matrix_power(U,j)@self.psi for j in range(r+1)])
        approximate=np.linalg.matrix_power(S,r)@self.psi
        out,trace=m.decompose(trajectory,self.reduced,T,approximate)
        phi=trajectory[-1]/np.linalg.norm(trajectory[-1])
        Q=np.eye(len(phi))-np.outer(phi,phi.conj())
        vectors=[Q@np.linalg.matrix_power(S,r-1-j)@(S-U)@trajectory[j] for j in range(r)]
        D=sum(m.norm2(v) for v in vectors)/m.norm2(approximate)
        X=sum(2*np.vdot(vectors[i],vectors[j]).real for i in range(r) for j in range(i+1,r))/m.norm2(approximate)
        self.assertAlmostEqual(out['projected_local_power'],D,places=13)
        self.assertAlmostEqual(out['projected_interference'],X,places=13)
        self.assertLess(out['state_reconstruction_residual'],1e-13)

    def test_phase_invariance_and_sector_split(self):
        exact=np.array([1.,0.,0.,0.],complex)
        approximate=np.array([np.sqrt(.98),.1,.1,0],complex)
        mask=np.array([True,True,False,False])
        a=m.physical_metrics(exact,approximate,mask,mask)
        b=m.physical_metrics(exact*np.exp(.8j),approximate*np.exp(-.6j),mask,mask)
        for key in ('one_minus_overlap','infidelity','inside_sector_error','outside_sector_error'):
            self.assertAlmostEqual(a[key],b[key],places=14)
        self.assertAlmostEqual(a['infidelity'],a['inside_sector_error']+a['outside_sector_error'],places=14)

    def test_small_error_not_subtractive_floor(self):
        exact=np.array([1.,0.],complex)
        approximate=np.array([np.sqrt(1-1e-20),1e-10],complex)
        mask=np.array([True,True])
        metrics=m.physical_metrics(exact,approximate,mask,mask)
        self.assertAlmostEqual(metrics['one_minus_overlap']/5e-21,1,places=12)
        self.assertTrue(metrics['below_analysis_floor'])

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'result.json';m.save(p,{'ok':True})
            with self.assertRaises(FileExistsError):m.save(p,{})


if __name__=='__main__':
    unittest.main()
