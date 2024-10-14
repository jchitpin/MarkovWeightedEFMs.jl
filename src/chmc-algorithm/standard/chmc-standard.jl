## Import files in this subfolder
include("error-checking.jl")
include("helpers.jl")

# CHMC: (unimolecular; computes EFM probabilities and weights)
"""
    steady_state_efm_distribution(#
        S::Matrix{<:Integer},
        v::Vector{<:Real},
        I::Int64 = 1;
        solver = nothing,
        issparse::Bool = false
    )

Enumerate the EFMs from stoichiometry matrix `S` and compute the
steady state probabilities of each EFM according to the discrete-time,
cycle-history Markov chain.

`S` is a fully-connected, unimolecular, m by n stoichiometry matrix with m
metabolites and n reactions.

`v` is the n-length steady state flux vector associated with *S*.

`I` is the initial starting state for rooting the cycle-history Markov chain.
The choice of initial starting state does not affect the steady state EFM
probabilities. The default is 1 and must be a whole number between 1:m.

`solver` is the type used for eigenvector calculations. Default is `nothing`
and `LinearSolve` will pick the best solver. Consult `LinearSolve.jl` for 
specifying other solvers.

`issparse` is true will use a sparse transition probability for the CHMC
to save memory.

# Example
```julia
julia> S = [#
 -1  0  0  0  0  0  0  0  0  0  1
  1 -1  1 -1  0  0  0  0  0  0  0
  0  1 -1  0 -1  1  0  0  0  0  0
  0  0  0  1  0  0 -1  0  0  0  0
  0  0  0  0  1 -1  1 -1  1 -1  0
  0  0  0  0  0  0  0  0  0  1 -1
  0  0  0  0  0  0  0  1 -1  0  0
];
julia> v = [3, 2, 1, 2, 3, 2, 2, 1, 1, 3, 3];
julia> res = steady_state_efm_distribution(S, v);
julia> res.e # EFM state sequences
6-element Vector{Vector{Int64}}:
 [3, 2, 3]
 [3, 2, 4, 5, 3]
 [3, 5, 3]
 [6, 1, 2, 4, 5, 6]
 [7, 5, 7]
 [6, 1, 2, 3, 5, 6]

julia> res.p # EFM probabilities
6-element Vector{Float64}:
 0.10638297872340426
 0.0425531914893617
 0.25531914893617025
 0.1914893617021277
 0.14893617021276595
 0.25531914893617025

julia> res.w # EFM weights
6-element Vector{Float64}:
 0.7142857142857142
 0.2857142857142857
 1.7142857142857144
 1.2857142857142858
 0.9999999999999999
 1.7142857142857144
```
"""
function steady_state_efm_distribution(#
    S::Matrix{<:Integer},
    v::Vector{<:Real},
    I::Int64 = 1;
    solver = nothing,
    issparse::Bool = false,
    verbose::Bool = true
)
    stime = time()

    # Error checking
    sanitize_initial(I, size(S, 1))

    # Pre-process open-loop stoichiometry matrix
    if verbose == true
        @info("(1/4): Checking whether input matrix S is open/closed loop.")
    end
    status = check_open_closed(S)
    if status == "open"
        if verbose == true
            @info("-> Open network detected.")
        end
        I = I + 1 # pseudo-metabolite row index is always 1
        S = close_network(S, v)
    end

    # Error check S and v and generate transition probability matrix
    if verbose == true
        @info("(2/4): Converting matrix S to transition probability matrix.")
    end
    T = stoichiometry_to_transition_matrix(S, v)
    if issparse == true
        T = sparse(T)
        dropzeros!(T)
    end

    # Construct CHMC, enumerate unimolecular EFMs, and compute EFM probabilities
    if verbose == true
        @info("(3/4): Enumerating EFMs and estimating their probabilities.")
    end
    res = steady_state_efm_distribution(T, I; solver = solver, verbose = verbose)
    e = res.e
    p = res.p

    # Compute EFM weights from proportionality constant
    if verbose == true
        @info("(4/4): Computing EFM weights.")
    end
    c = length.(e) .- 1
    α = sum(v) / sum(p .* c)
    w = p * α

    # Remove pseudo metabolite from EFMs in an open-network
    if status == "open"
        filter!.(i -> i != 1, e) # Remove '1' index for pseudo metabolite
        e = e .|> i -> i .- 1    # Subtract '1' from all EFM indices
    end

    if verbose == true
        etime = time()
        elapsed = round((etime - stime) / 60, digits = 2)
        @info("Total computation time: $elapsed min.")
    end

    return CHMCStandardSummary(e, p, w)
end

# CHMC (unimolecular; computes EFM probabilities)
"""
    steady_state_efm_distribution(#
        T::Union{Matrix{<:Real}, SparseMatrixCSC{Float64, Int64}},
        I::Int64 = 1;
        solver = nothing,
        verbose::Bool = true
    )

Enumerate the EFMs from (right) transition probability matrix `T` whose rows
sum to one, and compute the steady state probabilities of each EFM according to
the discrete-time, cycle-history Markov chain.

`T` is the discrete-time transition probability matrix with probabilities
proportional to the outgoing fluxes.

`I` is the initial starting state for rooting the cycle-history Markov chain.
The choice of initial starting state does not affect the steady state EFM
probabilities. The default is 1 and must be a whole number between 1:m.

`solver` is the type used for eigenvector calculations. Default is `nothing`
and `LinearSolve` will pick the best solver. Consult `LinearSolve.jl` for 
specifying other solvers.

# Example
```julia
julia> T = [#
  0.0  1.0   0.0       0.0  0.0   0.0  0.0
  0.0  0.0   0.5       0.5  0.0   0.0  0.0
  0.0  0.25  0.0       0.0  0.75  0.0  0.0
  0.0  0.0   0.0       0.0  1.0   0.0  0.0
  0.0  0.0   0.333333  0.0  0.0   0.5  0.166667
  1.0  0.0   0.0       0.0  0.0   0.0  0.0
  0.0  0.0   0.0       0.0  1.0   0.0  0.0
];
julia> res = steady_state_efm_distribution(T);
julia> res.e # EFM state sequences
6-element Vector{Vector{Int64}}:
 [3, 2, 3]
 [3, 2, 4, 5, 3]
 [3, 5, 3]
 [6, 1, 2, 3, 5, 6]
 [7, 5, 7]
 [6, 1, 2, 4, 5, 6]

julia> res.p # EFM probabilities
6-element Vector{Float64}:
 0.13723110896294213
 0.035797649943428565
 0.26989203316869914
 0.20302350628312624
 0.13926980198123254
 0.2147858996605714
```
"""
function steady_state_efm_distribution(#
    T::Union{Matrix{<:Real}, SparseMatrixCSC{Float64, Int64}},
    I::Int64 = 1;
    solver = nothing,
    verbose::Bool = true
)
    # Error checking
    sanitize_transition_matrix(T)

    # Construct cycle-history transition probability matrix and prefix dictionary
    T′, d = trie_matrix(T, I; verbose = verbose)

    # Enumerate all EFMs/simple cycles from the cycle-history matrix/prefix
    ϕ = enumerate_efms(T′, d; verbose = verbose)

    # Compute the steady state probabilities of being at a given state
    if verbose == true
        @info("       (7/7): Computing CHMC stationary distribution.")
    end
    p = solve_efm_probabilities(T′, ϕ, solver)

    return CHMCStandardSummary(first.(ϕ), p/sum(p), zeros(Int64,length(first.(ϕ))))
end

# Stoichiometry to transition probability matrix
"""
    stoichiometry_to_transition_matrix(S::Matrix{<:Integer}, v::Vector{<:Real})

Convert stoichiometry matrix `S` with vector of steady state fluxes to a
right stochastic transition probability matrix with rows summing to one.

`S` is the m by n stoichiometry matrix with m metabolites and n reactions.

`v` is the steady state flux vector with length n.

# Examples
```julia
julia> S = [#
  -1  0  0  0  0  0  0  0  0  0  1
   1 -1  1 -1  0  0  0  0  0  0  0
   0  1 -1  0 -1  1  0  0  0  0  0
   0  0  0  1  0  0 -1  0  0  0  0
   0  0  0  0  1 -1  1 -1 -1  1  0
   0  0  0  0  0  0  0  1  0  0 -1
   0  0  0  0  0  0  0  0  1 -1  0
]
julia> v = [2, 2, 2, 2, 2, 2, 2, 2, 4]
julia> stoich_to_transition(S, v)
7x7 Matrix{Float64}:
 0.0  1.0   0.0       0.0  0.0   0.0  0.0
 0.0  0.0   0.5       0.5  0.0   0.0  0.0
 0.0  0.25  0.0       0.0  0.75  0.0  0.0
 0.0  0.0   0.0       0.0  1.0   0.0  0.0
 0.0  0.0   0.333333  0.0  0.0   0.5 0.166667
 1.0  0.0   0.0       0.0  0.0   0.0  0.0
 0.0  0.0   0.0       0.0  1.0   0.0  0.0
```
"""
function stoichiometry_to_transition_matrix(S::Matrix{<:Integer}, v::Vector{<:Real})
    # Error-checking for a fully-connected, unimolecular network at steady state
    sanitize_stoich(S)
    sanitize_flux(v)
    sanitize_stoich_flux(S, v)

    # Initialize (right stochastic) transition probability matrix
    T = zeros(size(S, 1), size(S, 1))
 
    # Irreversible unimolecular reactions of index i --> index products
    idx_products = [findall(<(0), @view S[i, :]) for i in 1:size(S, 1)]

    # Fill probabilities of transition probability matrix
    for i in 1:size(T, 1)
        if !isempty(idx_products[i])
            # Indices of products in row i
            idx_j = vcat([findall(>(0), @view S[:,j]) for j in idx_products[i]]...)

            # Assign normalized flux probabilities to transition matrix row i
            T[i, idx_j] = v[idx_products[i]] / sum(v[idx_products[i]])
        end
    end

    # Error-checking that flux arcs are not all zero for a given metabolite state
    sanitize_transition_matrix(T)

    return T
end

# Convert matrix of EFMs to nested vector (unimolecular)
"""
    reshape_efm_matrix(ϕ::Matrix{Int64}, S::Matrix{<:Real})

Convert a matrix of EFMs `ϕ` to a nested vector of EFMs from a stoichiometry
matrix `S`. Stoichiometry matrix may only contain unimolecular reactions.

`ϕ` is the n by k EFM matrix with n reactions (rows) and k EFMs (cols).

`S` is the m by n stoichiometry matrix with m metabolites (rows) and n
reactions (cols).

# Examples
```julia
julia> ϕ = [#
  1  1  0  0  0  0
  1  0  1  0  0  0
  0  0  1  0  0  1
  0  1  0  0  0  1
  1  0  0  1  0  0
  0  0  0  1  0  1
  0  1  0  0  0  1
  1  1  0  0  0  0
  0  0  0  0  1  0
  0  0  0  0  1  0
  1  1  0  0  0  0
]
julia> S = [#
  -1  0  0  0  0  0  0  0  0  0  1
   1 -1  1 -1  0  0  0  0  0  0  0
   0  1 -1  0 -1  1  0  0  0  0  0
   0  0  0  1  0  0 -1  0  0  0  0
   0  0  0  0  1 -1  1 -1 -1  1  0
   0  0  0  0  0  0  0  1  0  0 -1
   0  0  0  0  0  0  0  0  1 -1  0
]
julia> efm_vector = reshape_efm_matrix(ϕ, S)
6-element Vector{Vector{Int64}}:
 [1, 2, 3, 5, 6, 1]
 [1, 2, 4, 5, 6, 1]
 [2, 3, 2]
 [3, 5, 3]
 [5, 7, 5]
 [3, 2, 4, 5, 3]
```
"""
function reshape_efm_matrix(ϕ::Matrix{Int64}, S::Matrix{<:Real})
    # Error checking
    sanitize_stoich(S)
    sanitize_efms(ϕ, S)

    # Convert EFM matrix to nested vector of EFMs
    nested_ϕ = [findall(!=(0), ϕ[:, j]) for j in 1:size(ϕ, 2)]

    # Ensure EFM reactions are correctly ordered according to S
    for e in nested_ϕ
        c = 1
        while c < length(e)
            if findfirst(>(0), S[:,e[c]]) == findfirst(<(0), S[:,e[c+1]])
                c += 1
            else
                push!(e, e[c+1])
                deleteat!(e, c+1)
            end
        end
    end

    # Convert reaction indices to metabolite state indices
    for i in eachindex(nested_ϕ)
        nested_ϕ[i] = [#
            findfirst(<(0), S[:,nested_ϕ[i][j]]) for j in eachindex(nested_ϕ[i])
        ]
    end

    # Ensure all EFMs loop back to their starting state
    for e in nested_ϕ
        if e[1] != e[end]
            push!(e, e[1])
        end
    end

    return nested_ϕ
end

# Convert nested vector of EFMs to matrix (unimolecular)
"""
    reshape_efm_vector(ϕ::Vector{Vector{Int64}}, S::Matrix{<:Real})

Convert nested vector of EFM indices `ϕ` with length k to an n by k matrix
of EFMs based on m by n stoichiometry matrix `S`. Stoichiometry matrix may
only contain unimolecular reactions.

`ϕ` is the nested vector of EFMs with length k and elements corresponding to
EFM metabolite indices in `S`.

`S` is the m by n stoichiometry matrix with m metabolites (rows) and n
reactions (cols).

# Examples
```julia
julia> ϕ = [#
  [1, 2, 3, 5, 6, 1],
  [1, 2, 4, 5, 6, 1],
  [2, 3, 2], [3, 5, 3], [5, 7, 5], [2, 4, 5, 3, 2]
]
julia> S = [#
  -1  0  0  0  0  0  0  0  0  0  1
   1 -1  1 -1  0  0  0  0  0  0  0
   0  1 -1  0 -1  1  0  0  0  0  0
   0  0  0  1  0  0 -1  0  0  0  0
   0  0  0  0  1 -1  1 -1 -1  1  0
   0  0  0  0  0  0  0  1  0  0 -1
   0  0  0  0  0  0  0  0  1 -1  0
]
julia> efm_matrix = reshape_efm_vector(ϕ, S)
11x6 Matrix{Int64}:
 1  1  0  0  0  0
 1  0  1  0  0  0
 0  0  1  0  0  1
 0  1  0  0  0  1
 1  0  0  1  0  0
 0  0  0  1  0  1
 0  1  0  0  0  1
 1  1  0  0  0  0
 0  0  0  0  1  0
 0  0  0  0  1  0
 1  1  0  0  0  0
```
"""
function reshape_efm_vector(ϕ::Vector{Vector{Int64}}, S::Matrix{<:Real})
    # Error checking
    sanitize_stoich(S)
    sanitize_efms(ϕ, S)

    # Ensure all EFMs loop back to their starting state
    for i in eachindex(ϕ)
        if ϕ[i][1] != ϕ[i][end]
            push!(ϕ[i], ϕ[i][1])
        end
    end

    # Vector of reactions and (i, j) coordinates of metabolite substrates/products
    subs = [findfirst(<(0), S[:, j]) for j in 1:size(S, 2)]
    prods = [findfirst(>(0), S[:, j]) for j in 1:size(S, 2)]
    reaction_pairs = tuple.(subs, prods)

    # Initialize and fill EFM matrix by converting EFM metabolite indices to rxns
    E = zeros(Int64, length(reaction_pairs), length(ϕ))
    for j in eachindex(ϕ)
        for i in 1:(length(ϕ[j])-1)
            k = (ϕ[j][i], ϕ[j][i+1])
            ii = findfirst(#
                [k == reaction_pairs[l] for l in eachindex(reaction_pairs)]
            )
            E[ii, j] = 1
        end
    end
    return E
end

