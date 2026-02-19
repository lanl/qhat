from datetime import datetime
import glob
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import coo_array, kron

# =================================================================================================

def generate_I():
    dataI = [np.complex128(1), np.complex128(1)]
    coordsI = [ [0,1], [0,1] ]
    PauliI = coo_array((dataI, coordsI))
    return PauliI
def generate_X():
    dataX = [np.complex128(1), np.complex128(1)]
    coordsX = [ [0,1], [1,0] ]
    PauliX = coo_array((dataX, coordsX))
    return PauliX
def generate_Y():
    dataY = [np.complex128(-1j), np.complex128(1j)]
    coordsY = [ [0,1], [1,0] ]
    PauliY = coo_array((dataY, coordsY))
    return PauliY
def generate_Z():
    dataZ = [np.complex128(1), np.complex128(-1)]
    coordsZ = [ [0,1], [0,1] ]
    PauliZ = coo_array((dataZ, coordsZ))
    return PauliZ

class PauliMatrices:
    _I = generate_I()
    _X = generate_X()
    _Y = generate_Y()
    _Z = generate_Z()
    @property
    def I(self):
        return self._I
    @property
    def X(self):
        return self._X
    @property
    def Y(self):
        return self._Y
    @property
    def Z(self):
        return self._Z

Pauli = PauliMatrices()

# =================================================================================================

class PauliTerm:
    _coefficient = None
    _string = None
    def __init__(self, c, s):
        self._coefficient = np.complex128(c)
        self._string = s
    @property
    def coefficient(self):
        return self._coefficient
    @property
    def matrix(self):
        return string_to_matrix(self._string)

# =================================================================================================

def char_to_matrix(c):
    if c == 'I':
        return Pauli.I
    elif c == 'X':
        return Pauli.X
    elif c == 'Y':
        return Pauli.Y
    elif c == 'Z':
        return Pauli.Z
    else:
        raise ValueError(f"Unknown Pauli matrix \"{c}\".")

def string_to_matrix(pauli_string):
    generator = (char_to_matrix(c) for c in pauli_string)
    matrix = next(generator)
    for next_pauli in generator:
        matrix = kron(matrix, next_pauli, format='coo')
    return matrix

def sum_of_pauli_strings_to_matrix(sum_of_pauli_strings):
    matrix = 0
    for item in sum_of_pauli_strings:
        matrix += item.coefficient * item.matrix
    return matrix

# =================================================================================================

def test(Nmax):
    sum_of_pauli_strings = list()
    for n in range(Nmax+1):
        coef = 0.5**n
        pauli_string = "I" * (Nmax - n) + "X" * n
        #print(f"  {coef:0.15f}  {pauli_string}")
        sum_of_pauli_strings.append(PauliTerm(coef, pauli_string))
    Hamiltonian = sum_of_pauli_strings_to_matrix(sum_of_pauli_strings)
    #print(Hamiltonian)
    #num_eig = min(Hamiltonian.shape[0]/2, 3)
    num_eig = 1
    largest = eigsh(Hamiltonian, k=num_eig, which="LA")
    #print(largest)
    #vl = np.matrix(largest[1])
    #assert vl.H @ vl == 1
    smallest = eigsh(Hamiltonian, k=num_eig, which="SA")
    #print(smallest)
    #vs = np.matrix(smallest[1])
    #assert vs.H @ vs == 1
    print(largest[0] / smallest[0])
    print(datetime.now())

# =================================================================================================

def coefficient_sums(sum_of_pauli_strings):
    accum_c = 0
    accum_abs_c = 0
    accum_cinv = 0
    accum_abs_cinv = 0
    for term in sum_of_pauli_strings:
        c = term.coefficient
        assert c.imag == 0
        c = c.real
        cinv = 1.0 / c
        accum_c += c
        accum_abs_c += abs(c)
        accum_cinv += cinv
        accum_abs_cinv += abs(cinv)
    return accum_c, accum_abs_c, 1.0 / accum_cinv, 1.0 / accum_abs_cinv

# =================================================================================================

def main(filename):
    print(f"[{datetime.now()}] -- Beginning processing", flush=True)
    sum_of_pauli_strings = list()
    with open(filename, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                tokens = line.split(' ')
                sum_of_pauli_strings.append(PauliTerm(tokens[0], tokens[1]))
    print(f"[{datetime.now()}] -- Loaded the sum of Pauli strings", flush=True)
    #print(sum_of_pauli_strings)
    one_norm = sum(abs(item.coefficient) for item in sum_of_pauli_strings)
    print(f"[{datetime.now()}] -- Computed the one-norm: {one_norm:.9f}", flush=True)
    Hamiltonian = sum_of_pauli_strings_to_matrix(sum_of_pauli_strings)
    print(f"[{datetime.now()}] -- Constructed the Hamiltonian matrix", flush=True)
    # TODO: Do we want the full spectrum, the largest and smallest values on the spectrum, the full
    #       decomposition, or the largest and smallest components of the decomposition?
    largest = eigsh(Hamiltonian, k=1, which="LA")[0][0]
    print(f"[{datetime.now()}] -- Computed the largest eigenvalue", flush=True)
    smallest = eigsh(Hamiltonian, k=1, which="SA")[0][0]
    print(f"[{datetime.now()}] -- Computed the smallest eigenvalue", flush=True)
    print(f"largest  eigenvalue = {largest:12.9f} <= {one_norm:12.9f}")
    print(f"smallest eigenvalue = {smallest:12.9f} >= {-1*one_norm:12.9f}")
    print(f"ratio of ranges = {(largest - smallest) / (2 * one_norm):.6f}")

# =================================================================================================

def as_size(filename):
    # ./library/001-001_H-H/1.00/sto-6g/H-H_1.00_sto-6g_as-002-002_jw.dat
    return int(filename[-14:-11]) + int(filename[-10:-7])

# =================================================================================================

if __name__ == "__main__":
    files = glob.glob("./library/*/*/*/*.dat")
    files.sort(key=as_size)
    count = len(files)
    for index, file in enumerate(files):
        print(f"[{datetime.now()}] Processing file {index} of {count}: \"{file}\"", flush=True)
        main(file)
    print(f"[{datetime.now()}] Processing complete.", flush=True)
