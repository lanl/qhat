import cirq
import numpy
import scipy


# =================================================================================================

def build_Pauli(string):
    def char_to_cirq(c):
        if c == 'I':
            return cirq.I
        elif c == 'X':
            return cirq.X
        elif c == 'Y':
            return cirq.Y
        elif c == 'Z':
            return cirq.Z
        else:
            raise ValueError(f"Pauli strings may only contain 'I', 'X', 'Y', 'Z'.  Found '{c}'.")
    N = len(string)
    q = cirq.LineQubit.range(N)
    ps = cirq.PauliString(char_to_cirq(string[0])(q[0]))
    for n in range(1,N):
        ps = ps * char_to_cirq(string[n])(q[n])
    return ps

# =================================================================================================

def test1():
    PS1 = build_Pauli("IIXZZXII")
    PS2 = build_Pauli("XZZZXIII")
    M = 5 * PS1 + 2 * PS2
    print(type(M.matrix()))
    w, vr = scipy.linalg.eig(M.matrix())
    for i in range(len(w)):
        print(w[i])
        print(vr[:,i])

# =================================================================================================

def build_term(components):
    assert len(components) > 0
    c, s = components[0]
    term = c * build_Pauli(s)
    for n in range(1, len(components)):
        c, s = components[n]
        term = term + c * build_Pauli(s)
    return term

# =================================================================================================

def build_ramp(terms, coefficient, time, reverse):
    # Be aware that PauliSumExponential only works if all Pauli strings in the sum commute (which
    # the function checks for us)
    assert len(terms) > 0
    ramp = list()
    if reverse:
        for term in reversed(terms):
            ramp.append(cirq.PauliSumExponential(term, coefficient * time))
    else:
        for term in terms:
            ramp.append(cirq.PauliSumExponential(term, coefficient * time))
    ramp = cirq.Circuit(ramp)
    return ramp

# =================================================================================================

def Trotterize(terms, coefficients, time):
    reverse = True
    circuit = cirq.Circuit()
    for coef in coefficients:
        reverse = not reverse
        circuit.append(build_ramp(terms, coef, time, reverse))
    return circuit

# =================================================================================================

def test2():
    terms = list()
    terms.append([(10.0, "IIII"), ])
    terms.append([( 1.0, "IXYZ"), ( 2.0, "IXXY"), ( 3.0, "IXZX")])
    terms.append([(-1.0, "ZYXI"), (-2.0, "YXXI"), (-3.0, "XZXI"), (-4.0, "ZYXI")])
    print(terms)
    for nt, term in enumerate(terms):
        print(f"Term {nt}:")
        for ps in term:
            c = ps[0]
            s = ps[1]
            print(f"  {c:4.1f} {s:4s}")
        print(f"  {build_term(term)}")
        for i in range(len(term)):
            for j in range(i+1, len(term)):
                if cirq.commutes(build_term(term[i:i+1]), build_term(term[j:j+1])):
                    print(f"  {term[i]} and {term[j]} commute")
                else:
                    print(f"  {term[i]} and {term[j]} DO NOT commute")
    op_terms = list(build_term(term) for term in terms)
    print(Trotterize(op_terms, (0.5,0.5), 10.0).unitary())

# =================================================================================================

def main():
    #test1()
    test2()

# =================================================================================================

if __name__ == "__main__":
    main()
