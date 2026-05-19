using LinearAlgebra
using SparseArrays

⊗(a,b) = kron(a,b)
### Gates
M = "IXYZ"
sparse_la = true
#Pauli matrices
if sparse_la == true
     global  Id = [1. 0; 0 1] |> sparse
     global  σx = [0. 1; 1 0] |> sparse
     global  σy = [0. -1im; 1im 0] |> sparse
     global  σz = [1. 0; 0 -1] |> sparse
    else
    global Id = [1. 0; 0 1]
    global σx = [0. 1; 1 0]
    global σy = [0. -1im; 1im 0]
    global σz = [1. 0; 0 -1]
end

#List of Pauli
P = [Id, σx, σy, σz]
PIZ = [Id, σz]

#Normalized Pauli matrices
Id_n = [1. 0; 0 1]/sqrt(2)
σx_n = [0. 1; 1 0]/sqrt(2)
σy_n = [0. -1im; 1im 0]/sqrt(2)
σz_n = [1. 0; 0 -1]/sqrt(2)


#List of normalized Pauli
P_n = [Id_n, σx_n, σy_n, σz_n]
PIZ_n = [Id_n, σz_n]

#Pauli list as dictionary
PDict = Dict()
PDict['X'] = σx
PDict['Y'] = σy
PDict['Z'] = σz
PDict['I'] = Id

PDict_n = Dict()
PDict_n['X'] = σx/sqrt(2)
PDict_n['Y'] = σy/sqrt(2)
PDict_n['Z'] = σz/sqrt(2)
PDict_n['I'] = Id/sqrt(2)

#Rotations list as dictionary
RotDict = Dict()
[RotDict[a] = exp(-0.25*im*π*Matrix(PDict[a])) for a in "XYZ"]

#Elimination rotations list as dictionary
RPZDict = Dict()
RPZDict['X'] = RotDict['Y']'
RPZDict['Y'] = RotDict['X']
RPZDict['Z'] =  Id
RPZDict['I'] =  Id

RZPDict = Dict()
RZPDict['X'] = RotDict['Y']
RZPDict['Y'] = RotDict['X']'
RZPDict['Z'] =  Id
RZPDict['I'] =  Id

PVecDict = Dict()
PVecDict['X'] = [0 1 0 0]
PVecDict['Y'] = [0 0 1 0]
PVecDict['Z'] = [0 0 0 1]
PVecDict['I'] = [1 0 0 0]

#Projectors on |0> and |1>
M0 =[1 0. ; 0 0]
M1 =[0 0. ; 0 1]



#Get operator from string of Pauli (normalized basis and unnormalized basis)
OP_from_string_n(Pst::String) = reduce(kron, [PDict_n[x] for x in Pst])
OP_from_string(Pst::String) = reduce(kron, [PDict[x] for x in Pst])
function String_from_number(P,n, OS; b = length(OS))
    r = digits(P, base=b, pad=n) .+ 1 |> reverse
    return  String([OS[x] for x in r])
end




function CNOT_OP(i,j,n)
    U1 = fill(Id, n)
    U2 = fill(Id, n)
    U1[i] = M0
    U2[i] = M1
    U2[j] = σx
    return reduce(kron, U1) + reduce(kron, U2)
end

function meas_mat(nq)
    H = (σx + σz)/sqrt(2)
    H_l = fill(H,nq)
    return  reduce(kron, H_l)
end
function PB_list(n)
    OS = "IXYZ"
    out = []
    for R in 0:(4^n-1)
        r = digits(R, base=4, pad=n) .+ 1 |> reverse
        Rst = String([OS[x] for x in r])
        push!(out, OP_from_string_n(Rst))
    end
    return  out
end



function PTM(G,n)
    pau_bas_P = PB_list(n)
    pau_bas_Q = PB_list(n)
    ##Construct PTM
    PTMatrix = [ real( tr(pau_bas_P[i]*G'*pau_bas_Q[j]*G) ) for i = 1:length(pau_bas_P), j = 1:length(pau_bas_Q) ]
    return PTMatrix
end


unorm(A) = max(abs.(A .- mean(A,dims=2))...)
function mnorm(ρ,n)
    μ = 0.0
    for P in 1:(4^n)-1
        Pst = String_from_number(P,n,"IXYZ")
        PO = OP_from_string_n(Pst)
        μ = max(μ, abs(tr(PO*ρ)))
    end
    return μ
end
