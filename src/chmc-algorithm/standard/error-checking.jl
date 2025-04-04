# Error checking for steady_state_efm_distribution (unimolecular)
function sanitize_stoich(S::Matrix{<:Integer})
    subs = length.([findall(<(0), S[:, i]) for i in 1:size(S, 2)])
    prds = length.([findall(>(0), S[:, i]) for i in 1:size(S, 2)])

    for j in 1:size(S, 2)
        @assert(#
            !(all(S[:, j] .== 0)),
            join([#
                "Reaction column in S cannot contain all zeros. ",
                "(A reaction must involve at least one metabolite ",
                "being produced or consumed)"
            ])
        )
    end
    for i in 1:size(S,1)
        @assert(#
            !(all(S[i,:] .== 0)),
            join([#
                "Metabolite row in S cannot contain all zeros. ",
                "(A metabolite must participate in at least one reaction."
            ])
        )
    end
    @assert(#
        all(subs .∈ Ref([0, 1])),
        join([#
            "Only unimolecular reactions allowed with one ",
            "stoichiometric unit of substrate."
        ])
    )
    @assert(#
        all(prds .∈ Ref([0, 1])),
        join([#
            "Only unimolecular reactions allowed with one ",
            "stoichiometric unit of product."
        ])
    )
    @assert(#
        all(prds - subs .== 0),
        join([#
            "Only unimolecular reactions allowed with one stoichiometric unit ",
            "of substrate converting to one stoichiometric unit of product."
        ])
    )

    Scols = mapslices(x -> [x], S, dims = 1)[:]
    @assert(#
        length(Scols) == length(unique(Scols)),
        "Reaction columns in the stoichiometry matrix must be unique."
    )
end

function sanitize_flux(v::Vector{<:Real})
    @assert(all(v .>= 0), "Fluxes must be ≥ 0.")
end

function sanitize_stoich_flux(S::Matrix{<:Integer}, v::Vector{<:Real})
    d = round.(S * v, digits = 5)
    @assert(#
        all(d .== 0.0),
        join([#
            "Network must be fully connected and ",
            "fluxes must satisfy metabolic steady state."
        ])
    )
end

function sanitize_transition_matrix(T::Matrix{<:Real})
    @assert(size(T, 1) == size(T, 2), "T must be a square matrix.")
    @assert(#
        all([round(sum(T[i, :]), digits = 5) == 1 for i in 1:size(T, 1)]),
        join([#
            "T[i, :] must be a right stochastic matrix with rows summing to ",
            "one. Check that the stoichiometry matrix is fully-connected and ",
            "includes non-zero fluxes connected to each metabolite. "
        ])
    )
    @assert(#
        all([T[i, i] == 0 for i in size(T, 1)]),
        "The diagonal entries of the transition probability matrix must be zero."
    )
end

function sanitize_transition_matrix(T::SparseMatrixCSC{Float64, Int64})
    @assert(size(T, 1) == size(T, 2), "T must be a square matrix.")
    @assert(#
        all([round(sum(T[i, :]), digits = 5) == 1 for i in 1:size(T, 1)]),
        join([#
            "T[i, :] must be a right stochastic matrix with rows summing to ",
            "one. Check that the stoichiometry matrix is fully-connected and ",
            "includes non-zero fluxes connected to each metabolite. "
        ])
    )
    @assert(#
        all([T[i, i] == 0 for i in size(T, 1)]),
        "The diagonal entries of the transition probability matrix must be zero."
    )
end

function sanitize_initial(I::Int64, s::Int64)
    @assert(0 < I <= s, "I must belong in size(S, 1).")
end

# Error checking for reshape_efm_vector/matrix
function sanitize_efms(ϕ::Vector{Vector{Int64}}, S::Matrix{<:Integer})
    @assert(#
        maximum(vcat(ϕ...)) <= size(S, 1),
        "EFM indices must be in 1:size(S, 1)."
    )
end

function sanitize_efms(ϕ::Matrix{Int64}, S::Matrix{<:Integer})
    @assert(all(ϕ .∈ Ref([-1, 0, 1])), "Only binary EFMs accepted.")
    @assert(#
        size(ϕ, 1) == size(S, 2),
        "Rows in E must equal columns in S (same number of reactions)."
    )
end

