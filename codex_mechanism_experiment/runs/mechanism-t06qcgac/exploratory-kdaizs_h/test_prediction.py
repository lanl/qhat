from pathlib import Path
import tempfile
import unittest
import mechanism_experiment as m
from validate_structure_and_prediction import defect_columns,commutator_bound

_cache=tempfile.TemporaryDirectory(prefix='qhat-prediction-tests-')
m.init_numerics(Path(_cache.name))
np=m.np
from scipy.linalg import expm
from openfermion import QubitOperator


class PredictionTests(unittest.TestCase):
    def test_small_commutator_not_silently_pruned(self):
        x=QubitOperator();x.terms={((0,'X'),):1e-14}
        z=QubitOperator();z.terms={((0,'Z'),):1.0}
        self.assertAlmostEqual(commutator_bound(x,z)/2e-14,1.0,places=14)

    def test_number_conserving_block(self):
        h=QubitOperator();h.terms={((0,'X'),(1,'X')):.2,((0,'Y'),(1,'Y')):.2}
        number=QubitOperator();number.terms={((0,'Z'),):-.5,((1,'Z'),):-.5}
        self.assertEqual(commutator_bound(h,number),0)

    def test_common_energy_shift_preserves_projected_decomposition(self):
        keys=[((0,'Z'),),((0,'X'),)]
        terms=m.matrixfree.compile_ordered_terms(keys,dict(zip(keys,[.7,.2])),1,1e-12)
        _,reduced,_=m.affine_space(terms,0)
        H,_=m.assemble_hamiltonian(reduced,2)
        initial=np.array([1.,0.],complex);T=.4;r=5
        results=[]
        for shift in (0.,3.4):
            shifted=reduced+[(shift,np.uint64(0),np.uint64(0),0)]
            trajectory=np.array([expm(-1j*(H.toarray()+shift*np.eye(2))*T*j/r)@initial for j in range(r+1)])
            tr=m.evolve(initial,shifted,T,r)
            out,_=m.decompose(trajectory,shifted,T,tr);results.append(out)
        for key in ('projected_local_power','projected_interference','coherence_factor'):
            self.assertAlmostEqual(results[0][key]/results[1][key],1.,places=9)

    def test_augmented_prediction_converges_to_product_error(self):
        keys=[((0,'Z'),),((0,'X'),)]
        coefficients=dict(zip(keys,[.7,.2]))
        terms=m.matrixfree.compile_ordered_terms(keys,coefficients,1,1e-12)
        states,reduced,_=m.affine_space(terms,0)
        H,parity=m.assemble_hamiltonian(reduced,len(states))
        D=defect_columns(reduced,np.arange(2),H,parity)
        A=np.block([[-1j*H.toarray(),D],[np.zeros((2,2)),-1j*H.toarray()]])
        initial=np.array([1.,0.],complex);T=.8
        propagated=expm(T*A)@np.concatenate([np.zeros(2,complex),initial])
        eta,exact=propagated[:2],propagated[2:]
        Q=np.eye(2)-np.outer(exact,exact.conj())
        errors=[]
        for r in (40,80):
            tr=m.evolve(initial,reduced,T,r)
            actual=Q@(tr-exact)
            predicted=T/r*(Q@eta)
            errors.append(np.linalg.norm(actual-predicted)/np.linalg.norm(actual))
        self.assertLess(errors[0],.05)
        self.assertLess(errors[1],.6*errors[0])


if __name__=='__main__':unittest.main()
