import inspect
import numpy as np

from pyLIQTR.ProblemInstances.ProblemInstance import ProblemInstance

def commute(term1, term2):
    # TODO
    # Each term has a Pauli string (term[0]) and a coefficient (term[1])
    return True

# TODO: OpenFermion's QubitOperator may be able to represent a group, and that may allow reusing
#       John's code without modification.
class LCPSHamiltonian(ProblemInstance):
    # Input format (TODO: check if this is the "right" format)
    # -- container of groups
    #    -- each group is a container of terms
    #       -- each term has a Pauli string and a coefficient
    #          -- entry 0: Pauli strings are dense string format (e.g., "IXZZYI")
    #          -- entry 1: coefficients are floats
    #       -- each term in a group must commute
    def __init__(self, groups):
        # Verify structure
        slen = None
        assert len(groups) > 0
        for group in groups:
            assert len(group) > 0
            for term in group:
                assert len(term) == 2
                pauli_string = term[0]
                if slen is None:
                    slen = len(pauli_string)
                else:
                    assert len(pauli_string) == slen
                for c in pauli_string:
                    assert c in "IXYZ"
                coefficient = term[1]
                assert coefficient == np.float64(coefficient)
        # Verify commutation
        for group in groups:
            for i in range(len(group)):
                for j in range(i+1, len(group)):
                    assert commute(group[i], group[j])
        # Store groups
        self._groups = groups

    def __str__(self):
        string = f"LCPSHamiltonian({self.n_qubits()}|{len(self._groups)}:"
        for group in self._groups:
            string = string + f"{len(group)},"
        return string[:-1] + ")"

    def n_qubits(self):
        return len(self._groups[0][0][0])

    def yield_PauliLCU_Info(self, return_as='arrays', do_pad=0, pad_value=1.0):
        # TODO: Currently don't actually support the arguments correctly
        assert return_as == 'strings'
        assert do_pad == 0
        assert pad_value == 1.0
        # Break things down into a single stream of Pauli terms (undoes grouping)
        for group in self._groups:
            for term in group:
                yield term

    def groups(self):
        return self._groups

# =================================================================================================

def test1():
    group1 = (
            ("IXYZ", 0.1),
            ("IXZY", 0.2),
            ("ZXZY", 0.3),
    )
    group2 = (
            ("XXXX", 0.4),
            ("YYYY", 0.5),
            ("ZZZZ", 0.6),
    )
    group3 = (
            ("IIII", 0.7),
            ("ZYXI", 0.8),
    )
    groups = (group1, group2, group3)
    problem_instance = LCPSHamiltonian(groups)
    print(problem_instance)

    assert problem_instance.n_qubits() == 4

    all_terms = list()
    for group in groups:
        all_terms.extend(group)
    n = 0
    for term in problem_instance.yield_PauliLCU_Info(return_as='strings'):
        assert term == all_terms[n]
        n += 1

    n = 0
    for group in problem_instance.groups():
        assert group == groups[n]
        n += 1

    print(f"test \"{inspect.currentframe().f_code.co_name}\" passed")

# =================================================================================================

def test2():
    group0 = (
            ("IIIIIIIII", 5.0),
    )
    group1 = (
            ("XIIIIIIII", 0.1),
            ("IXIIIIIII", 0.2),
            ("IIXIIIIII", 0.3),
            ("IIIXIIIII", 0.4),
            ("IIIIXIIII", 0.5),
            ("IIIIIXIII", 0.6),
            ("IIIIIIXII", 0.7),
            ("IIIIIIIXI", 0.8),
            ("IIIIIIIIX", 0.9),
    )
    group2 = (
            ("XXXXIXXXX", 0.25),
            ("YYYYIYYYY", 0.50),
            ("ZZZZIZZZZ", 0.75),
    )
    groups = (group0, group1, group2)
    problem_instance = LCPSHamiltonian(groups)
    print(problem_instance)

    assert problem_instance.n_qubits() == 9

    all_terms = list()
    for group in groups:
        all_terms.extend(group)
    n = 0
    for term in problem_instance.yield_PauliLCU_Info(return_as='strings'):
        assert term == all_terms[n]
        n += 1

    n = 0
    for group in problem_instance.groups():
        assert group == groups[n]
        n += 1

    print(f"test \"{inspect.currentframe().f_code.co_name}\" passed")

# =================================================================================================

def test3():
    import datetime
    from common.trotter import build_ramped_trotterized_unitary
    group1 = (
            ("IXYZ", 0.1),
            ("IXZY", 0.2),
            ("ZXZY", 0.3),
    )
    group2 = (
            ("XXXX", 0.4),
            ("YYYY", 0.5),
            ("ZZZZ", 0.6),
    )
    group3 = (
            ("IIII", 0.7),
            ("ZYXI", 0.8),
    )
    groups = (group1, group2, group3)
    problem_instance = LCPSHamiltonian(groups)
    print(problem_instance)

    print(f"{datetime.datetime.now()} | Calling build_ramped_trotterized_unitary")
    U0 = build_ramped_trotterized_unitary(problem_instance, "second order", 1.0, 1)
    print(f"{datetime.datetime.now()} | Calling U0.tensor_contract()")
    M0 = U0.tensor_contract()
    print(f"{datetime.datetime.now()} | Verifying M0.shape")
    assert M0.shape == (16,16)

    print(f"test \"{inspect.currentframe().f_code.co_name}\" passed")

# =================================================================================================

def main():
    test1()
    test2()
    test3()

# =================================================================================================

if __name__ == "__main__":
    main()
