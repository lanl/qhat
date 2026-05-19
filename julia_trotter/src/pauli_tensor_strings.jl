# pauli_tensor_strings.jl
using LinearAlgebra
using ITensors

include("quantum_utils.jl")  # expects PDict['I','X','Y','Z'] etc.

# ----------------------------------------------------------------------
# Build once: dense Complex Pauli matrices derived from quantum_utils.jl
# ----------------------------------------------------------------------
global PAULI_CHARS = ('I', 'X', 'Y', 'Z')

global PDict_dense = Dict{Char, Matrix{ComplexF64}}()
for p in PAULI_CHARS
  PDict_dense[p] = Matrix{ComplexF64}(PDict[p])
end

# ----------------------------------------------------------------------
# Build once: Pauli multiplication table derived from PDict_dense
# σ_a σ_b = phase * σ_c
# ----------------------------------------------------------------------
const PAULI_MUL = Dict{Tuple{Char,Char}, Tuple{ComplexF64,Char}}()

for a in PAULI_CHARS, b in PAULI_CHARS
  R = PDict_dense[a] * PDict_dense[b]

  best_err   = Inf
  best_phase = 1.0 + 0.0im
  best_c     = 'I'

  for c in PAULI_CHARS
    C = PDict_dense[c]
    phase = tr(C' * R) / tr(C' * C)
    err   = norm(R - phase * C)
    if err < best_err
      best_err   = err
      best_phase = phase
      best_c     = c
    end
  end

  best_err > 1e-10 && error("Failed to build Pauli multiplication table for ($a,$b), err=$best_err")
  PAULI_MUL[(a,b)] = (ComplexF64(best_phase), best_c)
end

@inline mul_pauli(a::Char, b::Char) = PAULI_MUL[(a,b)]

# ----------------------------------------------------------------------
# Self-contained local operator payload: stores the ITensor explicitly
# ----------------------------------------------------------------------
struct SitePauliOp
  site::Int       # 1..N
  p::Char         # 'X','Y','Z' (never 'I' stored)
  op::ITensor     # indices (prime(sites[site]), sites[site])
end

@inline function make_site_op(sites::Vector{<:Index}, site::Int, p::Char)::SitePauliOp
  s  = sites[site]
  sp = prime(s)
  return SitePauliOp(site, p, ITensor(PDict_dense[p], sp, s))
end

# ----------------------------------------------------------------------
# PauliTensorString (parametrized on Index type to avoid Vector invariance)
# Represents: coeff * ⊗_j σ_{p_j}, stored sparsely (skip identities)
# ----------------------------------------------------------------------
struct PauliTensorString{I<:Index}
  coeff::ComplexF64
  sites::Vector{I}
  ops::Vector{SitePauliOp}   # sorted by .site
end

# Construct from full Pauli string like "IZXY..."
function PauliTensorString(coeff::Number, sites::AbstractVector{I}, pstr::AbstractString) where {I<:Index}
  sitesv = collect(sites)  # ensure Vector{I}
  N = length(sitesv)
  length(pstr) == N || error("Pauli string length $(length(pstr)) must equal number of sites $N")

  ops = SitePauliOp[]
  for (j,ch) in enumerate(pstr)
    ch == 'I' && continue
    push!(ops, make_site_op(sitesv, j, ch))
  end
  # already increasing in j, but sort anyway for safety
  sort!(ops; by = o -> o.site)

  return PauliTensorString{I}(ComplexF64(coeff), sitesv, ops)
end

# Construct from sparse list (site => Pauli char)
function PauliTensorString(coeff::Number, sites::AbstractVector{I}, ops_in::AbstractVector{<:Pair{Int,Char}}) where {I<:Index}
  sitesv = collect(sites)

  ops = SitePauliOp[]
  for (j,p) in ops_in
    p == 'I' && continue
    push!(ops, make_site_op(sitesv, j, p))
  end
  sort!(ops; by = o -> o.site)

  return PauliTensorString{I}(ComplexF64(coeff), sitesv, ops)
end

# Pretty helper: return the full-length Pauli string
function paulistring(P::PauliTensorString)
  N = length(P.sites)
  s = fill('I', N)
  for o in P.ops
    s[o.site] = o.p
  end
  return String(s)
end

# ----------------------------------------------------------------------
# 1) Product with a state ITensor (state has unprimed site indices)
# ----------------------------------------------------------------------
function apply_pauli_string(P::PauliTensorString, ψ::ITensor)::ITensor
  x = ψ
  @inbounds for o in P.ops
    x = noprime(o.op * x)
  end
  return x
end

function Base.:*(P::PauliTensorString, ψ::ITensor)::ITensor
  return P.coeff * apply_pauli_string(P, ψ)
end

# scalar scaling
Base.:*(a::Number, P::PauliTensorString{I}) where {I<:Index} =
  PauliTensorString{I}(ComplexF64(a)*P.coeff, P.sites, P.ops)

Base.:*(P::PauliTensorString, a::Number) = a * P

# ----------------------------------------------------------------------
# 2) Product of two PauliTensorStrings (operator composition)
# (A*B)*ψ == A*(B*ψ)
# ----------------------------------------------------------------------
function Base.:*(A::PauliTensorString{I}, B::PauliTensorString{I}) where {I<:Index}
  A.sites == B.sites || error("PauliTensorStrings must have identical `sites` to multiply")

  opsA = A.ops
  opsB = B.ops
  i = 1
  j = 1

  out_ops = SitePauliOp[]
  phase = 1.0 + 0.0im
  sites = A.sites

  while i <= length(opsA) || j <= length(opsB)
    if j > length(opsB) || (i <= length(opsA) && opsA[i].site < opsB[j].site)
      push!(out_ops, opsA[i])
      i += 1
    elseif i > length(opsA) || opsB[j].site < opsA[i].site
      push!(out_ops, opsB[j])
      j += 1
    else
      site = opsA[i].site
      # local product order is A then B
      (ph, c) = mul_pauli(opsA[i].p, opsB[j].p)
      phase *= ph

      if c != 'I'
        push!(out_ops, make_site_op(sites, site, c))
      end

      i += 1
      j += 1
    end
  end

  # out_ops is already sorted due to merge process
  return PauliTensorString{I}(A.coeff * B.coeff * phase, sites, out_ops)
end

# ----------------------------------------------------------------------
# 3) Time evolution: exp(-i*m*P*t) * ψ
# For P = coeff * S, with S^2 = I (Pauli string):
# exp(-i a S) ψ = cos(a) ψ - i sin(a) (Sψ), where a = m*t*coeff
# ----------------------------------------------------------------------
function time_evolve(P::PauliTensorString, ψ::ITensor, t::Real; m::Real=1.0)::ITensor
  a  = (m * t) * P.coeff
  Sψ = apply_pauli_string(P, ψ)   # apply S only (no coeff)
  return (cos(a) * ψ) - (1im * sin(a) * Sψ)
end



# ----------------------------------------------------------------------
# 4) Commutator between Pauli Strings 
# Uses the fact that product of pauli strings is always another Pauli string. So we need to only look at the coefficients. And the ordering only changes the coefficients
# ----------------------------------------------------------------------



function commutator(P::PauliTensorString, Q::PauliTensorString)::PauliTensorString
    A = P*Q
    B = Q*P
    return PauliTensorString(A.coeff - B.coeff,A.sites, A.ops)
end






