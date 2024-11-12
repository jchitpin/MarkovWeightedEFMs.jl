## Import files in this subfolder
include("error-checking.jl")
include("helpers.jl") # includes structs

"""
    canonicalize_smiles(m::Vector{String})

Canonicalize input SMILES string `m` for RXNMapper mapping.
"""
function canonicalize_smiles(m::Vector{String})
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.smiles_utils
    return rxn_mapper.canonicalize_smi.(m)
end

"""
    canonicalize_smiles(m::String)

Canonicalize input SMILES string `m` for RXNMapper mapping.
"""
function canonicalize_smiles(m::String)
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.smiles_utils
    return rxn_mapper.canonicalize_smi(m)
end

# For printing summary results below in the following format
print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))

# Print summary to console
import Base.print
"""
    print(res::CHMCAtomicErrorSummary)

Write CHMCAtomicErrorSummary to `stdout`.
"""
function print(res::CHMCAtomicErrorSummary)
    summary = join([#
        "############################################################\n",
        "## ERROR CHECKING STOICHIOMETRY MATRIX AND FLUX VECTOR #####\n",
        "# (1)  SUM OF ABSOLUTE FLUX RECONSTRUCTION ERROR:\n",
        "#      %s\n",
        "#      %s\n",
        "# (2)  REACTIONS THAT ARE DUPLICATES:\n",
        "#      %s\n",
        "#      %s\n",
        "# (3)  REACTIONS WTIH ZERO FLUX:\n",
        "#      %s\n",
        "#      %s\n",
        "# (4)  REACTIONS WTIH NEGATIVE FLUX:\n",
        "#      %s\n",
        "#      %s\n",
        "# (5)  INTERNAL REACTIONS W/ NON-INTEGER STOICHIOMETRIES:\n",
        "#      %s\n",
        "#      %s\n",
        "# (6)  UNIMOLECULAR SOURCE REACTIONS W/ STOICH == 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (7)  UNIMOLECULAR SOURCE REACTIONS W/ STOICH != 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (8)  MULTIMOLECULAR SOURCE REACTIONS W/ STOICH == 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (9)  MULTIMOLECULAR SOURCE REACTIONS W/ STOICH != 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (10) UNIMOLECULAR SINK REACTIONS W/ STOICH == 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (11) UNIMOLECULAR SINK REACTIONS W/ STOICH != 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (12) MULTIMOLECULAR SINK REACTIONS W/ STOICH == 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (13) MULTIMOLECULAR SINK REACTIONS W/ STOICH != 1:\n",
        "#      %s\n",
        "#      %s\n",
        "# (14) REACTIONS W/ NO SUBSTRATES OR PRODUCTS:\n",
        "#      %s\n",
        "#      %s\n",
        "# (15) # METABOLITES PARTICIPATING IN NO REACTIONS:\n",
        "#      %s\n",
        "#      %s\n",
        "# STATUS:\n",
        "#      %s\n",
        "############################################################\n"
    ])

    outcome = pass_fail(res)
    status = ""
    if all(outcome .== "PASSED.")
        status = "PASSED. THESE INPUTS SATISFY ATOMIC CHMC REQUIREMENTS."
    else
        status = "FAILED. THESE INPUTS FAIL ATOMIC CHMC REQUIREMENTS."
    end
    interleaved = [#
        collect(Iterators.flatten(zip(string_errorsummary(res), outcome)));
        status
    ]
    print_formatted(summary, interleaved...)
end

# Return summary of errors in stoichiometry matrix/flux vector/smiles
"""
    find_atomic_chmc_input_errors(S::Matrix{<:Real}, v::Vector{<:Real})

Performs a series of checks on stoichiometry matrix `S` and steady state flux
vector `v` to determine whether they meet the atomic CHMC requirements. Returns
an `CHMCAtomicErrorSummary` structure that can be printed via [`print(res::CHMCAtomicErrorSummary)`](@ref).
"""
function find_atomic_chmc_input_errors(S::Matrix{<:Real}, v::Vector{<:Real})
    # Flux-specific checks
    res1 = find_atomic_chmc_input_errors_fluxes(S, v)

    # Stoichiometry-specific checks
    res2 = find_atomic_chmc_input_errors_rest(S)

    res = CHMCAtomicErrorSummary(#
        res1[1],
        res2[1],
        res1[2],
        res1[3],
        res2[2],
        res2[3],
        res2[4],
        res2[5],
        res2[6],
        res2[7],
        res2[8],
        res2[9],
        res2[10],
        res2[11],
        res2[12]
    )

    return res
end
"""
    find_atomic_chmc_input_errors(S::Matrix{<:Real})

Performs a series of checks on stoichiometry matrix `S` to determine whether
they meet the atomic CHMC requirements for atomic EFM enumeration. Returns
an `CHMCAtomicErrorSummary` structure that can be printed via [`print(res::CHMCAtomicErrorSummary)`](@ref).
"""
function find_atomic_chmc_input_errors(S::Matrix{<:Real})
    # Stoichiometry-specific checks
    res2 = find_atomic_chmc_input_errors_rest(S)

    res = CHMCAtomicErrorSummary(#
        0.0,
        res2[1],
        Vector{Int64}(),
        Vector{Int64}(),
        res2[2],
        res2[3],
        res2[4],
        res2[5],
        res2[6],
        res2[7],
        res2[8],
        res2[9],
        res2[10],
        res2[11],
        res2[12]
    )

    return res
end

# Fix errors identified from |find_atomic_chmc_input_errors|
"""
    correct_atomic_chmc_input_errors(#
        res::CHMCAtomicErrorSummary,
        S::Matrix{<:Real},
        v::Vector{<:Real},
        mets::Vector{String},
        rxns::Vector{String}
    )

Corrects all errors identified by [`find_atomic_chmc_input_errors`](@ref) except
for the steady state flux requirement. Returns a corrected version of each input.

`S` is the m by n stoichiometry matrix.
`v` is the steady state flux vector of length n.
`mets` is the vector of metabolite names of length m.
`rxns` is the vector of reaction names of length n.
"""
function correct_atomic_chmc_input_errors(#
    res::CHMCAtomicErrorSummary,
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String}
)
    if sanitize_correct_atomic_chmc_input_errors(res, S, mets, rxns) == true
        return nothing
    end

    return correct_atomic_chmc_input_errors_inner(res, S, v, mets, rxns)
end

"""
    correct_atomic_chmc_input_errors(#
        res::CHMCAtomicErrorSummary,
        S::Matrix{<:Real},
        mets::Vector{String},
        rxns::Vector{String}
    )

Corrects all errors identified by [`find_atomic_chmc_input_errors`](@ref) except
for the steady state flux requirement. Returns a corrected version of each input.

`S` is the m by n stoichiometry matrix.
`mets` is the vector of metabolite names of length m.
`rxns` is the vector of reaction names of length n.
"""
function correct_atomic_chmc_input_errors(#
    res::CHMCAtomicErrorSummary,
    S::Matrix{<:Real},
    mets::Vector{String},
    rxns::Vector{String}
)
    if sanitize_correct_atomic_chmc_input_errors(res, S, mets, rxns) == true
        return nothing
    end
    v = zeros(size(S,2))
    S, _, mets, rxns = correct_atomic_chmc_input_errors_inner(res, S, v, mets, rxns)

    return S, mets, rxns
end

# Remove reactions due to RXNMapper character limitations/pseudometabolites
"""
    correct_atomic_chmc_input_smiles(#
        S::Matrix{<:Real},
        v::Vector{<:Real},
        mets::Vector{String},
        rxns::Vector{String},
        smiles::Vector{String},
        H::Bool=false
    )
Updates the inputs by (1) removing pseudometabolites and (2) reactions whose
reaction strings exceed the RXNMapper character limit. Pseudometabolites are
defined as either metabolites with no defined chemical structure/SMILES
representation, or metabolites whose SMILES strings contain symbols absent
from the periodic table (including metabolites with an R-group). Reactions
that exceed the RXNMapper character limit are converted into unimolecular
sources and sinks to maintain steady state flux. Returns the corrected inputs
with a list of removed metabolites/reactions.

`S` is the m by n stoichiometry matrix.

`v` is the steady state flux vector of length n.

`mets` is the vector of metabolite names of length m.

`rxns` is the vector of reaction names of length n.

`smiles` is the vector of metabolite smiles matching `mets`.

`H = false` will not check for hydrogen conservation in each reaction.
"""
function correct_atomic_chmc_input_smiles(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    H::Bool = false
)
    return correct_atomic_chmc_input_smiles_inner(S, v, mets, rxns, smiles, H)
end

"""
    correct_atomic_chmc_input_smiles(#
        S::Matrix{<:Real},
        mets::Vector{String},
        rxns::Vector{String},
        smiles::Vector{String},
        H::Bool=false
    )
Updates the inputs by (1) removing pseudometabolites and (2) reactions whose
reaction strings exceed the RXNMapper character limit. Pseudometabolites are
defined as either metabolites with no defined chemical structure/SMILES
representation, or metabolites whose SMILES strings contain symbols absent
from the periodic table (including metabolites with an R-group). Reactions
that exceed the RXNMapper character limit are converted into unimolecular
sources and sinks to maintain steady state flux. Returns the corrected inputs
with a list of removed metabolites/reactions.

`S` is the m by n stoichiometry matrix.

`mets` is the vector of metabolite names of length m.

`rxns` is the vector of reaction names of length n.

`smiles` is the vector of metabolite smiles matching `mets`.

`H = false` will not check for hydrogen conservation in each reaction.
"""
function correct_atomic_chmc_input_smiles(#
    S::Matrix{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    H::Bool = false
)
    # Initialize flux vector to reuse existing functions
    v = zeros(length(rxns))

    S, _, mets, rxns, smiles, rm = correct_atomic_chmc_input_smiles_inner(#
        S, v, mets, rxns, smiles, H
    )

    return S, mets, rxns, smiles, rm
end

"""
    exchange_atomic_chmc_input_metabolites(#
        S::Matrix{<:Real},
        v::Vector{<:Real},
        mets::Vector{String},
        rxns::Vector{String},
        smiles::Vector{String},
        ii::Vector{Int64}
    )
Removes internal metabolites `ii` from the remaining inputs by converting
each metabolite i ∈ ii into a pair of sink and source metabolites. Returns an
updated copy of inputs `S`, `v`, `mets`, `rxns`, and `smiles`.

`S` is the m by n stoichiometry matrix.

`v` is the steady state flux vector of length n.

`mets` is the vector of metabolite names of length m.

`rxns` is the vector of reaction names of length n.

`smiles` is the vector of metabolite smiles matching `mets`.

`ii` is the vector of metabolite indices to exchange for pairs of sinks/sources.
"""
function exchange_atomic_chmc_input_metabolites(#
    S::Matrix{Int16},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    ii::Vector{Int64}
)
    # Error checking
    sanitize_exchange(S, ii)
    @assert(length(mets) == size(S,1))
    @assert(length(smiles) == size(S,1))
    @assert(length(rxns) == size(S,2))
    @assert(length(v) == size(S,2))

    return exchange_atomic_chmc_input_metabolites_inner(S, v, mets, rxns, smiles, ii)
end
"""
    exchange_atomic_chmc_input_metabolites(#
        S::Matrix{<:Real},
        mets::Vector{String},
        rxns::Vector{String},
        smiles::Vector{String},
        ii::Vector{Int64}
    )
Removes internal metabolites `ii` from the remaining inputs by converting
each metabolite i ∈ ii into a pair of sink and source metabolites. Returns an
updated copy of inputs `S`, `mets`, `rxns`, and `smiles`.

`S` is the m by n stoichiometry matrix.

`mets` is the vector of metabolite names of length m.

`rxns` is the vector of reaction names of length n.

`smiles` is the vector of metabolite smiles matching `mets`.

`ii` is the vector of metabolite indices to exchange for pairs of sinks/sources.
"""
function exchange_atomic_chmc_input_metabolites(#
    S::Matrix{Int16},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    ii::Vector{Int64}
)
    # Error checking
    sanitize_exchange(S, ii)
    @assert(length(mets) == size(S,1))
    @assert(length(smiles) == size(S,1))
    @assert(length(rxns) == size(S,2))

    v = zeros(size(S,2))
    S, _, mets, rxns, smiles = exchange_atomic_chmc_input_metabolites_inner(#
        S, v, mets, rxns, smiles, ii
    )

    return S, mets, rxns, smiles
end

# Get source metabolites for main function
"""
    get_source_metabolites(S::Matrix{Int16})

Returns a vector of all source metabolites in stoichiometry matrix `S`.

`S` is the m by n stoichiometry matrix.
"""
function get_source_metabolites(S::Matrix{Int16})
    source_mets = Vector{Int64}()
    for j in 1:size(S, 2)
        neg_stoich = S[findall(<(0), @view S[:, j]), j]
        pos_stoich = S[findall(>(0), @view S[:, j]), j]
        if isempty(neg_stoich) && length(pos_stoich) == 1 && pos_stoich[1] == 1
            push!(source_mets, findfirst(>(0), @view S[:, j]))
        end
        if isempty(neg_stoich) && length(pos_stoich) == 1 && pos_stoich[1] > 1
            @warn(join([#
                "Metabolite row index $(findfirst(>(0), @view S[:, j])) ",
                "should have a source stoichiometry of 1."
            ]))
        end
    end
    return sort(source_mets)
end

# Get maximum number of designated atom types in a vector of SMILES
function get_max_atoms_old(s::Vector{String}, a::Symbol)
    sanitize_atom(String(a))
    f(x) = Dict(parse_formula(molecularformula(smilestomol(x))))
    x = f.(s)
    y = zeros(Int64, length(s))
    for i in eachindex(x)
        if haskey(x[i], String(a))
            y[i] = x[i][String(a)]
        end
    end
    return y
end
"""
    get_max_atoms(s::Vector{String}, a::Symbol)

Returns the maximum number of atoms of type `a` for each SMILES string in `s`.
Returns a vector of all source metabolites in stoichiometry matrix `S`.

`s` is the vector of metabolite SMILES strings.

`a` is a periodic table element symbol.
"""
function get_max_atoms(s::Vector{String}, a::Symbol)
    a = String(a)
    sanitize_atom(a)
    y = zeros(Int64, length(s))
    for i in eachindex(s)
        y[i] = length(findall(==(a), parse_from_1(parse_from_2(s[i]))))
    end
    return y
end

"""
    map_reaction_strings(#
        S::Matrix{<:Real},
        smiles::Vector{String},
        rxns::Vector{String},
        H::Bool=false
    )

Constructs reaction strings from the `smiles` and then returns the atom-mapped
reaction strings via RXNMapper.

`S` is the stoichiometry matrix.

`smiles` is the vector of SMILES strings for metabolites in `S`.

`rxns` is the vector of reaction names in `S`.

`H = false` will not check for hydrogen conservation in each reaction.
"""
function map_reaction_strings(#
    S::Matrix{<:Real},
    smiles::Vector{String},
    rxns::Vector{String},
    H::Bool = false
)
    # Error checking
    sanitize_smiles(smiles)

    # Construct reaction smiles (checks if each reaction is atom balanced)
    rs = rxn_string(S, smiles, rxns, H)

    # Construct atom tracing strings
    ms = trace_rxn_string(rs)

    return rs, ms
end

# Precompute atom tracing dictionary
"""
    precompute_atom_tracing_dictionary(#
        S::Matrix{Int16},
        ms::Vector{String},
        amax::Vector{Int64},
        a::Symbol
    )

Precomputes a dictionary mapping each atom in a given substrate stoichiometric
copy to its product atom across every reaction. The dictionary is used in
[`steady_state_efm_distribution`](@ref).

`S` is the stoichiometry matrix.

`ms` is the vector of mapped reaction SMILES strings.

`amax` is the total number of atom `a` in each metabolite row of `S`.

`a` is a periodic table element symbol.
"""
function precompute_atom_tracing_dictionary(#
    S::Matrix{Int16},
    ms::Vector{String},
    amax::Vector{Int64},
    a::Symbol
)
    sanitize_atom(String(a))
    d = Dict{Tuple{Int64, Int64, Int64, Int64}, Tuple{Int64, Int64}}()
    for j in findall(.!isempty.(ms))
        for sub_idx in findall(<(0), @view S[:, j])
            for sub_foi in 1:amax[sub_idx]
                for sub_copy in 1:abs(S[sub_idx, j])
                    sub_rsi = get_srsi(S[:, j], Int16(sub_idx), sub_copy)
                    prod_rsi, prod_foi = trace_atoms(#
                        ms[j],
                        (Int16(sub_rsi), Int16(sub_foi), String(a))
                    )
                    prod_idx = get_prod_meti(S[:, j], prod_rsi)
                    d[(sub_idx, sub_foi, sub_copy, j)] = (prod_idx, prod_foi)
                end
            end
        end
    end

    return d
end

# CHMC (higher order; computes EFM probabilities and weights)
import MarkovWeightedEFMs.CHMC.Standard.steady_state_efm_distribution
"""
    steady_state_efm_distribution(#
        S::Matrix{<:Integer},
        v::Vector{<:Real},
        ms::Vector{String},
        I::Tuple{Int64,Int64,Symbol},
        D::Dict{#
          NTuple{4,Int64},
          Tuple{Int64,Int64}
        } = Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}();
        solver = nothing,
        verbose::Bool = true
    )

Enumerates the atomic EFMs from stoichiometry matrix `S` and compute the
steady state probabilities of each EFM according to the discrete-time,
cycle-history Markov chain.

`S` is a fully-connected, unimolecular, m by n stoichiometry matrix with m
metabolites and n reactions.

`v` is the n-length steady state flux vector associated with *S*.

`ms` is the vector of mapped reaction SMILES strings.

`I` is a triplet of the metabolite row index in `S`, the atom index, and the
atom index type. It is the initial starting state for rooting the cycle-history
Markov chain.

`D` is the precomputed dictionary from [`precompute_atom_tracing_dictionary`](@ref).

`solver` is the type used for eigenvector calculations. Default is `nothing`
and `LinearSolve` will pick the best solver. Consult `LinearSolve.jl` for 
specifying other solvers.
"""
function steady_state_efm_distribution(#
    S::Matrix{<:Integer},
    v::Vector{<:Real},
    ms::Vector{String},
    I::Tuple{Int64,Int64,Symbol},
    D::Dict{#
        NTuple{4,Int64},
        Tuple{Int64,Int64}
    } = Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}();
    solver = nothing,
    verbose::Bool = true
)
    stime = time()
    II = (I[1], I[2], String(I[3]))

    # Error checking
    sanitize_stoichiometry_fluxes(S, v)
    sanitize_stoichiometry_strings(S, ms)
    sanitize_initial(II, S, ms)

    # Check if network is open or closed
    if verbose == true
        @info("(1/4): Checking that input matrix S is open loop.")
    end
    status = check_open_closed(Int64.(S))
    @assert(status == "open")

    # Construct CHMC directly from stoichiometry matrix and reaction SMILES
    if verbose == true
        @info("(2/4): Constructing atomic CHMC.")
    end

    T, d, invD, R = construct_atomic_chmc(S, v, ms, II, D, verbose = verbose)

    # Enumerate all EFMs/simple cycles from the cycle-history matrix/prefix
    if verbose == true
        @info("(3/4): Enumerating EFMs and estimating their probabilities.")
    end
    ϕ = enumerate_efms(T, d, verbose = verbose)
    e = first.(ϕ)

    # Compute the steady state probabilities of each atomic EFM
    if verbose == true
        @info("       Computing CHMC stationary distribution.")
    end
   
    p = solve_efm_probabilities(T, ϕ)
    #π = linearsolve_eigenvector(T, solver)
    #@assert(all(.!isnan.(π)), "Eigenvector cannot contain NaNs.")
    #p = solve_probabilities(T, ϕ, π)
    #@info p
    #@assert(all(.!isnan.(p)), "Steady state CHMC probabilities cannot be NaN.")

    # Compute EFM weights from proportionality constant
    if verbose == true
        @info("(4/4): Computing EFM weights.")
    end
    v_ext = collect(keys(invD))[findfirst(==((0,0)), collect(values(invD)))]
    idx = findall(x -> v_ext ∈ x, e)
    α = v[findfirst(==(1), S[I[1], :])] / sum(p[idx])
    w = p * α

    if verbose == true
        etime = time()
        elapsed = round((etime - stime) / 60, digits = 2)
        @info("Total computation time: $elapsed min.")
    end

    return CHMCAtomicSummary(I, ϕ, p / sum(p), w, invD, d, T, R)
end

"""
    enumerate_atomic_efms(#
        S::Matrix{<:Integer},
        ms::Vector{String},
        I::Tuple{Int64,Int64,Symbol},
        D::Dict{#
          NTuple{4,Int64},
          Tuple{Int64,Int64}
        } = Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}();
        verbose::Bool = true
    )

Enumerates the atomic EFMs from stoichiometry matrix `S` via the
cycle-history Markov chain. Returns the EFMs as MC states and their corresponding
simple cycles as CHMC states, and the MC and CHMC dictionaries.

`S` is an open-loop m by n stoichiometry matrix with m metabolites and
n reactions.

`ms` is the vector of mapped reaction SMILES strings.

`I` is a triplet of the metabolite row index in `S`, the atom index, and the
atom index type. It is the initial starting state for rooting the cycle-history
Markov chain.

`D` is the precomputed dictionary from [`precompute_atom_tracing_dictionary`](@ref).
"""
function enumerate_atomic_efms(#
    S::Matrix{Int16},
    ms::Vector{String},
    I::Tuple{Int64,Int64,Symbol},
    D::Dict{NTuple{4,Int64},Tuple{Int64,Int64}};
    verbose::Bool = true
)
    stime = time()

    # Error checking
    sanitize_stoichiometry_strings(S, ms)
    sanitize_atom(String(I[3]))

    # Check if network is open or closed
    @assert(check_open_closed(Int64.(S)) == "open")

    # Construct CHMC directly from stoichiometry matrix and reaction SMILES
    if verbose == true
        @info("(1/2): Constructing atomic CHMC.")
    end
    II = (I[1], I[2], String(I[3]))
    v = ones(Float64, size(S,2))
    T, d, invD, R = construct_atomic_chmc(S, v, ms, II, D, verbose = verbose)

    # Enumerate all EFMs/simple cycles from the cycle-history matrix/prefix
    if verbose == true
        @info("(2/2): Enumerating EFMs.")
    end
    ϕ = enumerate_efms(T, d, verbose = verbose)
    #p = zeros(length(ϕ))
    #w = zeros(length(ϕ))
    p = nothing
    w = nothing

    if verbose == true
        etime = time()
        elapsed = round((etime - stime) / 60, digits = 2)
        @info("Total computation time: $elapsed min.")
    end

    return CHMCAtomicSummary(I, ϕ, p, w, invD, d, T, R)
end

# Get sequence of metabolites in atomic CHMC EFMs
"""
    get_efm_metabolite_atom_indices(res::CHMCAtomicSummary, i::Int64)

Converts the sequence of metabolite indices for atomic EFM index `i` into
metabolite names. The length of the vector of metabolite names is one element
shorter than `res.e[i]` because the pseudometabolite connecting sink and source
reactions is omitted.

`res` are the results from [`steady_state_efm_distribution`](@ref).

`i` is the index for a given EFM.
"""
function get_efm_metabolite_atom_indices(res::CHMCAtomicSummary, i::Int64)
    @assert(i <= length(res.e), "EFM index $i exceeds number of EFMs in res.")
    e = first.(res.e)
    x = [res.dmc[e[i][k]] for k in 1:length(e[i])]
    filter!(y -> y != (0, 0), x)
    return x
end

"""
    get_efm_reaction_indices(res::CHMCAtomicSummary, i::Int64)

Converts the sequence of metabolite indices for atomic EFM index `i` into
reaction names. 

`res` are the results from [`steady_state_efm_distribution`](@ref).

`i` is the index for a given EFM.
"""
function get_efm_reaction_atom_indices(res::CHMCAtomicSummary, i::Int64)
    ks = collect(keys(res.dchmc))
    vs = collect(values(res.dchmc))
    e = first.(res.e)

    function issubarray(needle, haystack)
        getView(vec, i, len) = view(vec, i:i+len-1)
        ithview(i) = getView(haystack, i, length(needle))
        return any(i -> ithview(i) == needle, 1:length(haystack)-length(needle)+1)
    end

    # Find the CHMC state containing the sequence of MC states for the given EFM
    c = 1
    while c > 0
        if issubarray(e[i][2:end], ks[c]) && ks[c][end] == e[i][end]
            break
        else
            c += 1
        end
    end

    # Find the CHMC state for the EFM with 'd' fewer MC states in the sequence
    x = Vector{Int64}(undef, length(e[i])-1)
    x[1] = parse(Int64, first(vs[c])[2:end])
    d = 0
    for k in 2:length(x)
        d += 1
        e = 1
        while e > 0
            if issubarray(ks[c][1:(end-d)], ks[e]) && ks[e][end] == ks[c][end-d]
                x[k] = parse(Int64, first(vs[e])[2:end])
                break
            else
                e += 1
            end
        end
    end
    reverse!(x)

    y = Vector{Int16}(undef, length(x))
    for k in 1:length(x)
        if k != length(x)
            y[k] = res.R[x[k], x[k+1]]
        else
            y[k] = res.R[x[k], x[1]]
        end
    end

    return y
end

"""
    chmc_to_mc_matrix(res::CHMCAtomicSummary, v::Vector{<:Real})

Converts CHMC transition matrix `res.T` to a Markov chain with probabilities
taken from steady state flux vector `v`.

`res` are the results from [`steady_state_efm_distribution`](@ref).

`v` is the steady state flux vector used as input from [`steady_state_efm_distribution`](@ref).
"""
function chmc_to_mc_matrix(res::CHMCAtomicSummary, v::Vector{<:Real})
    if res.T isa SparseMatrixCSC
        T′ = spzeros(length(res.dmc), length(res.dmc))
    else
        T′ = zeros(length(res.dmc), length(res.dmc))
    end
    ks = last.(collect(keys(res.dchmc)))
    vs = first.(collect(values(res.dchmc)))
    vs = [parse(Int64, v[2:end]) for v in vs]

    for i in 1:size(res.R, 1)
        for j in 1:size(res.R, 2)
            if !iszero(res.R[i, j])
                idx = minimum(findall(==(i), vs)) # CHMC state i
                jdx = minimum(findall(==(j), vs)) # CHMC state j
                T′[ks[idx], ks[jdx]] = v[res.R[i, j]]
            end
        end
    end
    function normalize_rows(A::SparseMatrixCSC)
        A = permutedims(A, [2, 1])
        sums = sum(A, dims = 1)
        I,J,V = findnz(A)
        for idx in 1:length(V)
            V[idx] /= sums[J[idx]]
        end
        return permutedims(sparse(I, J, V), [2, 1])
    end
    if res.T isa SparseMatrixCSC
        T′ = normalize_rows(T′)
    else
        T′ = T′ ./ sum(T′, dims = 2)
    end
    dropzeros!(T′)
    return T′
end

"""
    preprocess_all_for_atomic_chmc(#
        S::Matrix{<:Real},
        v::Vector{<:Real},
        mets::Vector{String},
        rxns::Vector{String},
        smiles::Vector{String},
        atom::Symbol
    )

Wrapper function to run pre-processing functions and return the updated
model, atomic info, and error/warning logs.

`S` is the m by n stoichiometry matrix.
`v` is the steady state flux vector of length n.
`mets` is the vector of metabolite names of length m.
`rxns` is the vector of reaction names of length n.
`smiles` is the vector of metabolite SMILES strings of length m.
`a` is a periodic table element symbol.
"""
function preprocess_all_for_atomic_chmc(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    atom::Symbol
)
    errors = find_atomic_chmc_input_errors(S, v)
    #print(errors) # summary of errors associated with S/v
    output = correct_atomic_chmc_input_errors(errors, S, mets, rxns)
    if !isnothing(output)
        S = output[1]
        mets = output[2]
        rxns = output[3]
    end
    S, v, mets, rxns, smiles, logs = correct_atomic_chmc_input_smiles(S, v, mets, rxns, smiles)
    smiles = canonicalize_smiles(smiles)
    rs, ms = map_reaction_strings(S, smiles, rxns, false)
    atom_max = get_max_atoms(smiles, atom)
    D_atom = precompute_atom_tracing_dictionary(S, ms, atom_max, atom)

    # Identify source metabolites
    src_mets = get_source_metabolites(S)
    max_src_met_atoms = atom_max[src_mets]

    model_data = (#
        S = S, v = v, mets = mets, rxns = rxns, smiles = smiles, rs = rs,
        ms = ms
    )
    atom_info = (#
        atom = atom, src_mets = src_mets,
        max_src_met_atoms = max_src_met_atoms, D = D_atom
    )
    logs = (model_errors = errors, smiles_warnings = logs)

    return model_data, atom_info, logs
end

