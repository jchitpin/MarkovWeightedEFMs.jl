## Periodic table elements
function string_chemical_element()
    r = [#
        "A[cglmrstu]",
        "B[aehikr]?",
        "C[adeflmnorsu]?",
        "D[bsy]",
        "E[rsu]",
        "F[elmr]?",
        "G[ade]",
        "H[efgos]?",
        "I[nr]?",
        "Kr?",
        "L[airuv]",
        "M[dgnot]",
        "N[abdeiop]?",
        "Os?",
        "P[abdmortu]?",
        "R[abefghnu]",
        "S[bcegimnr]?",
        "T[abcehilm]",
        "U(u[opst])?",
        "V",
        "W",
        "Xe",
        "Yb?",
        "Z[nr]"
    ]
    return r
end

function regex_chemical_element()
    return Regex(join(string_chemical_element(), "|"))
end

# Error checking
function sanitize_stoichiometry(S::Matrix{<:Real})
    for j in 1:size(S, 2)
        subs = findall(<(0), S[:, j])
        prods = findall(>(0), S[:, j])

        # Internal reactions must have integer stoichiometries
        if length(subs) > 0 && length(prods) > 0
            @assert(#
                all(isinteger.(S[subs, j])),
                join([#
                    "Internal reaction stoichiometries must be integers. ",
                    "Check reaction $j."
                ])
            )
            @assert(#
                all(isinteger.(S[prods, j])),
                join([#
                    "Internal reaction stoichiometries must be integers. ",
                    "Check reaction $j."
                ])
            )
        end
        # Source reactions must be unimolecular with a stoichiometry of 1
        if length(subs) == 0 && length(prods) > 0
            @assert(#
                length(prods) == 1,
                "Source reaction $j must be unimolecular."
            )
            @assert(#
                S[prods[1], j] == 1,
                "Source reaction $j must have a stoichiometry of 1."
            )
        end
        # Sink reactions must be unimolecular with a stoichiometry of 1
        if length(subs) > 0 && length(prods) == 0
            @assert(#
                length(subs) == 1,
                "Sink reaction $j must be unimolecular."
            )
            @assert(#
                S[subs[1], j] == -1,
                "Sink reaction $j must have a stoichiometry of -1."
            )
        end
    end
    Scols = mapslices(x -> [x], S, dims = 1)[:]
    @assert(#
        length(Scols) == length(unique(Scols)),
        "Reaction columns in the stoichiometry matrix must be unique."
    )
end

function sanitize_smiles(s::Vector{String})
    # SMILES strings must be canonicalized
    s′ = canonicalize_smiles(s)
    if s′ != s
        @assert(#
            s′ == s,
            join([#
                "The input SMILES strings must be canonicalized by function ",
                "canonicalize_smiles()."
            ])
        )
    end

    # SMILES strings cannot contain periods
    for i in eachindex(s)
        @assert(#
            isnothing(match(r"\.", s[i])),
            "Input SMILES string at index $i cannot contain a period."
        )
    end
end

function sanitize_stoichiometry_fluxes(S::Matrix{<:Real}, v::Vector{<:Real})
    @assert(#
        length(v) == size(S, 2),
        "The number of fluxes in v must match the number of reactions in S."
    )
end

function sanitize_stoichiometry_strings(S::Matrix{<:Integer}, ms::Vector{String})
    @assert(#
        length(ms) == size(S, 2),
        join([#
            "The number of mapped reaction smiles strings must match ",
            "the number of reactions in S."
        ])
    )
end

function sanitize_exchange(S::Matrix{<:Integer}, ii::Vector{Int64})
    @assert(#
        all(ii .<= size(S, 1)),
        "Elements in ii must match metabolite rows in S."
    )
end

function sanitize_atom(a::String)
    @assert(#
        !isnothing(match(regex_chemical_element(), a)),
        "Invalid periodic table element symbol."
    )
    @assert(#
        a != "[H+]",
        "Cannot construct atomic CHMC for hydrogen atoms."
    )
    @assert(#
        a != "H",
        "Cannot construct atomic CHMC for hydrogen atoms."
    )
end

function sanitize_atom_smiles(y::Vector{Char}, a::String)
    @assert(#
        haskey(countmemb(parse_from_1(y)), a),
        join([#
            "The specified root atom type is not present in the "
            "root metabolite."
        ])
    )
end

function sanitize_initial(#
    I::Tuple{Int64,Int64,String},
    S::Matrix{<:Integer},
    ms::Vector{String}
)
    # Check for valid source metabolite index
    srcs = get_source_metabolites(Int16.(S))
    @assert(I[1] ∈ srcs, "I[3] must correspond to a source metabolite index.")

    # Check the atom type
    sanitize_atom(I[3])

    # Find smiles string of metabolite in S
    j = findall(!=(0), S[I[1],:])
    if all(isempty.(ms[j]))
        @assert(#
            any(.!isempty.(ms[j])),
            join([#
                "Metabolite I[1] has no corresponding mapped reaction SMILES. ",
                "This metabolite participates in no internal reaction and ",
                "has no meaningful atomic CHMC."
            ])
        )
    end

    ij = (I[1], j[findfirst(!isempty, ms[j])])
    if S[ij[1], ij[2]] < 0
        mets, _  = split(ms[ij[2]], ">>")
        met_indices = findall(<(0), S[:, ij[2]])
    else
        _, mets  = split(ms[ij[2]], ">>")
        met_indices = findall(>(0), S[:, ij[2]])
    end
    met_stoichs = abs.(S[met_indices, ij[2]])
    idx = findfirst(==(I[1]), met_indices)
    mets = split(mets, ".")
    m = mets[sum(met_stoichs[1:idx])]

    amax = get_max_atoms([string(m)], Symbol(I[3]))[1]
    @assert(#
        I[2] <= amax,
        join([#
            "Atom I[3] at index I[2] does not exist in the SMILES string of ",
            "metabolite I[1]."
        ])
    )
end

function sanitize_correct_atomic_chmc_input_errors(#
    res::CHMCAtomicErrorSummary,
    S::Matrix{<:Real},
    mets::Vector{String},
    rxns::Vector{String}
)
    # Error checking
    @assert(size(S, 1) == length(mets), "length(mets) must equal rows in S.")
    @assert(size(S, 2) == length(rxns), "length(rxns) must equal cols in S.")
    outcome = pass_fail(res)
    if outcome[1] == "FAILED."
        error(#
            join([#
                "Fluxes are not at steady state.\n",
                "It is required that: sum(abs.(S * v)) <= sqrt(eps())"
            ])
        )
    elseif all(outcome .== "PASSED.")
        @info("Inputs do not require correction.")
        return true
    else
        return false
    end
end

function sanitize_correct_atomic_chmc_input_smiles(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
)
    res = find_atomic_chmc_input_errors(S)
    outcome = pass_fail(res)
    if !all(outcome .== "PASSED.")
        @assert(#
            all(outcome .== "PASSED."),
            join([#
                "Input stoichiometry matrix S and flux vector v fail ",
                "function find_atomic_chmc_input_errors()"
            ])
        )
    end
    @assert(#
        size(S, 1) == length(smiles),
        join([#
            "Ensure the array of SMILES strings matches the metabolites ",
            "ordered in the stoichiometry matrix."
        ])
    )
    @assert(size(S, 1) == length(mets), "length(mets) must equal rows in S.")
    @assert(size(S, 2) == length(rxns), "length(rxns) must equal cols in S.")
    @assert(size(S, 2) == length(v), "length(v) must equal cols in S.")
end

