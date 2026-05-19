
"""
Parse a file with:

  # key = value
  ...
  <coeff> <PauliString>

Returns:
  metadata::Dict{String,String}
  ham::Dict{String,ComplexF64}   (PauliString => coefficient)

Notes:
- Metadata values are kept as strings (no typing/inference).
- Coeff supports real scientific notation and complex like 1.2+0.3i (or ...j).
- Duplicate Pauli strings are summed.
"""
function parse_hamiltonian_file(path::AbstractString)
    metadata = Dict{String,String}()
    ham = Dict{String,ComplexF64}()

    open(path, "r") do io
        for (lineno, raw) in enumerate(eachline(io))
            line = strip(raw)
            isempty(line) && continue

            # Metadata lines
            if startswith(line, "#")
                body = strip(replace(line, "#" => ""; count=1))
                isempty(body) && continue

                if occursin("=", body)
                    parts = split(body, "="; limit=2)
                    key = strip(parts[1])
                    val = strip(parts[2])
                    metadata[key] = val
                else
                    # Comment without '=' — keep it as a key with empty value
                    metadata[strip(body)] = ""
                end
                continue
            end

            # Hamiltonian lines: "<coeff> <paulistring>"
            toks = split(line)
            length(toks) < 2 && error("Line $lineno: expected '<coef> <paulistring>', got: '$raw'")

            coeff_str = toks[1]
            pauli = join(toks[2:end], "")  # in case spacing is weird; usually it's one token

            coeff = parse_coeff(coeff_str; lineno=lineno)

            ham[pauli] = get(ham, pauli, 0.0 + 0.0im) + coeff
        end
    end

    return metadata, ham
end

"""
Parse numeric coefficient token into ComplexF64.

Supports:
- real: -3.180033e+00
- complex: 1.2+0.3im, 1.2+0.3i, 1.2+0.3j
"""
function parse_coeff(s::AbstractString; lineno::Int=0)::ComplexF64
    t = strip(String(s))

    # Allow 'i' or 'j' suffix for imaginary unit; normalize to Julia's "im" form.
    # If it's already got "im", keep it.
    if occursin("im", t)
        return ComplexF64(parse(Complex{Float64}, t))
    end

    # If token contains i/j, translate to im
    if occursin('i', t) || occursin('j', t)
        # Replace trailing 'i'/'j' with "im" in forms like "a+bi" or "bi"
        # Also handle "a+bi" by replacing all 'i'/'j' with "" then appending "im" appropriately.
        # Easiest robust approach: replace i/j with "im" and let Julia parse Complex.
        t2 = replace(t, "j" => "im", "i" => "im")
        try
            return ComplexF64(parse(Complex{Float64}, t2))
        catch
            lineno == 0 ? error("Could not parse coefficient: '$s'") :
                          error("Line $lineno: could not parse coefficient: '$s'")
        end
    end

    # Otherwise it's real
    try
        return ComplexF64(parse(Float64, t), 0.0)
    catch
        lineno == 0 ? error("Could not parse real coefficient: '$s'") :
                      error("Line $lineno: could not parse real coefficient: '$s'")
    end
end

# --- Optional: CLI usage ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia parser.jl <path>")
        exit(1)
    end
    path = ARGS[1]
    meta, ham = parse_hamiltonian_file(path)
    println("metadata keys: ", length(meta))
    println("hamiltonian terms: ", length(ham))
end
