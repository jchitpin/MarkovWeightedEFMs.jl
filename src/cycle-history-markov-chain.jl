# Error checking helpers for steady_state_(c/d)tmc_efm_distribution
function sanitize_stoich(S::Matrix{<:Int64})
  subs = length.([findall(<(0), S[:,i]) for i in 1:size(S,2)])
  prods = length.([findall(>(0), S[:,i]) for i in 1:size(S,2)])
  @assert(#
    all(subs .∈ Ref([0, 1])),
    join([#
      "Only unimolecular reactions allowed with one stoichiometric unit of ",
      "substrate."
    ])
  )
  @assert(#
    all(prods .∈ Ref([0, 1])),
    join([#
      "Only unimolecular reactions allowed with one stoichiometric unit of ",
      "product."
    ])
  )
  @assert(#
    all(prods - subs .== 0),
    join([#
      "Only unimolecular reactions allowed with one stoichiometric unit of ",
      "substrate converting to one stoichiometric unit of product."
    ])
  )
end
function sanitize_flux(v::Vector{<:Real})
  @assert(all(v .>= 0), "Fluxes must be ≥ 0.")
end
function sanitize_stoich_flux(S::Matrix{<:Int64}, v::Vector{<:Real})
  d = S * v
  d = round.(d, digits=5)
  @assert(#
    all(d .== 0.0),
    join([#
      "Network must be fully connected and ",
      "fluxes must satisfy metabolic steady state."
    ])
  )
end
function sanitize_transition_matrix(T::Matrix{<:AbstractFloat})
  @assert(size(T,1) == size(T,2), "T must be a square matrix.")
  @assert(#
    all([round(sum(T[i,:]), digits = 5) == 1 for i in 1:size(T,1)]),
    join([#
      "T[i,:] must be a right stochastic matrix with rows summing to one. ",
      "Check that the stoichiometry matrix is fully-connected and ",
      "includes non-zero fluxes connected to each metabolite. "
    ])
  )
end
function sanitize_initial(I::Int64, s::Int64)
  @assert(0 < I <= s, "I must belong in size(T,1).")
end
function sanitize_efms(ϕ::Vector{Vector{Int64}}, S::Matrix{<:Real})
  @assert(#
    maximum(vcat(ϕ...)) <= size(S,1),
    "EFM indices must be 1:size(S,1)."
  )
end
function sanitize_efms(ϕ::Matrix{Int64}, S::Matrix{<:Real})
  @assert(all(ϕ .∈ Ref([-1, 0, 1])), "Only binary EFMs accepted.")
  @assert(#
    size(ϕ,1) == size(S,2),
    "Rows in E must equal columns in S (same number of reactions)."
  )
end
function sanitize_efms(ϕ::Vector{Vector{Int64}}, s::Int64)
  for e in ϕ
    @assert(all(e .∈ Ref(1:s)), "EFM $e does not match transition matrix T.")
  end
end

# Stoichiometry to transition probability matrix
"""
    stoich_to_transition(S::Matrix{<:Real}, v::Vector{<:Real})

Convert a stoichiometry matrix with vector of steady state fluxes to a
right stochastic transition probability matrix with rows summing to one.

*S* is the m by n stoichiometry matrix with m metabolites and n reactions.

*v* is the steady state flux vector with length n.

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
function stoich_to_transition(S::Matrix{<:Real}, v::Vector{<:Real})

  # Error-checking for a fully-connected, unimolecular network at steady state
  sanitize_stoich(S)
  sanitize_flux(v)
  sanitize_stoich_flux(S, v)

  # Initialize (right stochastic) transition probability matrix
  T = zeros(size(S,1), size(S,1))
 
  # Irreversible unimolecular reactions of index i --> index products
  idx_products = [findall(<(0), S[i,:]) for i in 1:size(S,1)]

  # Fill probabilities of transition probability matrix
  for i in 1:size(T,1)
    if !isempty(idx_products[i])
      # Indices of products in row i
      idx_j = vcat([findall(>(0), S[:,j]) for j in idx_products[i]]...)

      # Assign normalized flux probabilities to transition matrix row i
      T[i, idx_j] = v[idx_products[i]] / sum(v[idx_products[i]])
    end
  end

  # Error-checking that flux arcs are not all zero for a given metabolite state
  sanitize_transition_matrix(T)

  return T
end

# Construct cycle-history Markov chain
function trie_matrix(#
  T::Matrix{<:AbstractFloat},
  I::Int64=1,
)

# Error checking I
  sanitize_initial(I, size(T, 1))

  # Construct trie to get number of nodes
  t = trie(T, I)

  # Initialize new transition probability matrix
  A = zeros(length(t), length(t))

  # Traverse tree and populate new transition probability matrix with trie nodes
  function traverse_trie(#
    prefix::Vector{Int64},
    t::Dict{#
      Vector{Int64},
      NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
    },
    T::Matrix{<:AbstractFloat},
    A::Matrix{Float64}
  )
    i = prefix[end] # old state
    ii = parse(Int64, t[prefix].id[2:end]) # new state

    # Fill matrix with non-children (upstream) transitions
    childs = t[prefix].children
    parents = findall(>(0), T[i,:])
    parents = parents[parents .∉ Ref(childs)]

    for j in parents
      r = T[i,j]
      l = prefix[1:findfirst(x -> x == j, prefix)]
      jj = parse(Int64, t[l].id[2:end])
      A[ii,jj] = r
    end

    # Fill matrix with children transitions
    for j in childs
      r = T[i,j]
      jj = parse(Int64, t[[prefix; j]].id[2:end])
      A[ii,jj] = r
      traverse_trie([prefix; j], t, T, A)
    end
  end
  traverse_trie([I], t, T, A)

  return A, t
end

# Enumerate EFMs via counting simple cycles
function enumerate_efms(#
  T′::Matrix{Float64},
  d::Dict{#
    Vector{Int64},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
  }
)

  # Prefixes with cycle-history states
  prefixes = collect(keys(d))
  prefixes_transformed = Vector{Vector{Int64}}(undef, length(prefixes))
  for i in 1:length(prefixes)
    temp = Vector{Int64}()
    for j in 1:length(prefixes[i])
      push!(#
        temp,
        parse(Int64, d[prefixes[i][1:j]].id[2:end])
      )
    end
    prefixes_transformed[i] = temp
  end

  # For each prefix sequence, check if the last state can complete a simple
  # cycle by revisiting a previous state in the sequence
  # These are the EFM simple cycles
  simple_cycles_transformed = Vector{Vector{Int64}}()
  for i in 1:length(prefixes)
    upstream = findall(>(0), T′[prefixes_transformed[i][end],:])
    for j in 1:length(upstream)
      idx = findfirst(prefixes_transformed[i] .== upstream[j])
      if ~isnothing(idx)
        cycle = [prefixes_transformed[i][end]; prefixes_transformed[i][idx:end]]
        push!(simple_cycles_transformed, cycle)
      end
    end
  end

  # Convert simple cycle history states back to regular states
  e = Dict(parse(Int64, d[k].id[2:end]) => k[end] for k in keys(d))
  simple_cycles_original = Vector{Vector{Int64}}()
  for i in 1:length(simple_cycles_transformed)
    temp = Vector{Int64}()
    for j in 1:length(simple_cycles_transformed[i])
      push!(#
        temp,
        e[simple_cycles_transformed[i][j]]
      )
    end
    push!(simple_cycles_original, temp)
  end

  # Aggregate cycle history cycles over the same EFM
  res = NamedTuple{#
    (:EFM, :TransformedCycles),
    Tuple{Vector{Int64}, Vector{Vector{Int64}}}
  }[]
  while !isempty(simple_cycles_original)
    idx = findall(#
      in(Ref(Set(simple_cycles_original[1][2:end]))),
      [Set(i[2:end]) for i in simple_cycles_original]
    )
    push!(res, (EFM=simple_cycles_original[1], TransformedCycles=simple_cycles_transformed[idx]))
    splice!(simple_cycles_original, idx)
    splice!(simple_cycles_transformed, idx)
  end
  return res
end

# Main function (count EFMs and compute steady state probabilities/weights)
"""
    steady_state_efm_distribution(#
        S::Matrix{<:Int64},
        v::Vector{<:Real},
        I::Int64=1
    )

Enumerates the EFMs from stoichiometry matrix S and compute the
steady state probabilities of each EFM according to the cycle-history,
discrete-time Markov chain model.

*S* is a fully-connected, unimolecular, m by n stoichiometry matrix with m
metabolites and n reactions.

*v* is the n-length steady state flux vector associated with *S*.

*I* is the initial starting state for rooting the cycle-history Markov chain.
The choice of initial starting state does not affect the steady state EFM
probabilities. The default is 1 and must be a whole number between 1:m.

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
 [5, 3, 5]
 [7, 5, 7]
 [6, 1, 2, 4, 5, 6]
 [3, 2, 4, 5, 3]
 [6, 1, 2, 3, 5, 6]

julia> res.p # EFM probabilities
6-element Vector{Float64}:
 0.10638297872340426
 0.25531914893617025
 0.14893617021276595
 0.25531914893617025
 0.0425531914893617
 0.1914893617021277

julia> res.w # EFM weights
6-element Vector{Float64}:
 0.7142857142857142
 1.7142857142857144
 0.9999999999999999
 1.7142857142857144
 0.2857142857142857
 1.2857142857142858
```
"""
function steady_state_efm_distribution(#
  S::Matrix{<:Int64},
  v::Vector{<:Real},
  I::Int64=1
)

  # Check S and v and generate transition probability matrix
  T = stoich_to_transition(S, v)

  # Construct cycle-history transition probability matrix and prefix dictionary
  T′, d = trie_matrix(T, I)

  # Enumerate all EFMs/simple cycles from the cycle-history matrix/prefix
  ϕ = enumerate_efms(T′, d)
  e = [ϕ[i].EFM for i in 1:length(ϕ)]

  # Compute the steady state probabilities of being at a given state
  mc = MarkovChain(T′)
  π = stationary_distributions(mc)[1]

  # Compute steady state edge probabilities and aggregate for each EFM
  p = Vector{Float64}(undef, length(ϕ))
  for i in 1:length(ϕ)
    for j in ϕ[i].TransformedCycles
      p[i] += π[j[1]] * T′[j[1], j[2]]
    end
  end
  p = p ./ sum(p)

  # Compute EFM weights from proportionality constant
  c = [length(ϕ[i].EFM)-1 for i in 1:length(ϕ)]
  α = sum(v) / sum(p .* c)
  w = p * α

  return (e=e, p=p, w=w)
end

# Convert matrix of EFMs to nested vector (rows are reactions; cols are EFMs)
"""
    reshape_efm_matrix(#
        ϕ::Matrix{Int64},
        S::Matrix{<:Real}
    )

Convert a matrix of EFMs ϕ to a nested vector of EFMs from a stoichiometry
matrix S. Stoichiometry matrix may only contain
unimolecular reactions.

*ϕ* is the n by k EFM matrix with n reactions (rows) and k EFMs (cols).

*S* is the m by n stoichiometry matrix with m metabolites (rows) and n
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
function reshape_efm_matrix(#
  ϕ::Matrix{Int64},
  S::Matrix{<:Real}
)

  # Error checking
  sanitize_stoich(S)
  sanitize_efms(ϕ, S)

  # Convert EFM matrix to nested vector of EFMs
  nested_ϕ = [findall(!=(0), ϕ[:,j]) for j in 1:size(ϕ,2)]

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
  for i in 1:length(nested_ϕ)
    nested_ϕ[i] = [findfirst(<(0), S[:,nested_ϕ[i][j]]) for j in 1:length(nested_ϕ[i])]
  end

  # Ensure all EFMs loop back to their starting state
  for e in nested_ϕ
    if e[1] != e[end]
      push!(e, e[1])
    end
  end

  return nested_ϕ
end

# Convert nested vector of EFMs to matrix (rows are reactions; cols are EFMs)
"""
    reshape_efm_vector(#
        ϕ::Vector{Vector{Int64}},
        S::Matrix{<:Real}
    )

Convert nested vector of EFM indices ϕ with length k to an n by k matrix
of EFMs based on m by n stoichiometry matrix S. Stoichiometry matrix may
only contain unimolecular reactions.

*ϕ* is the nested vector of EFMs with length k and elements corresponding to
EFM metabolite indices in S.

*S* is the m by n stoichiometry matrix with m metabolites (rows) and n
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
function reshape_efm_vector(#
  ϕ::Vector{Vector{Int64}},
  S::Matrix{<:Real}
)

  ## Error checking
  sanitize_stoich(S)
  sanitize_efms(ϕ, S)

  # Ensure all EFMs loop back to their starting state
  for i in 1:length(ϕ)
    if ϕ[i][1] != ϕ[i][end]
      push!(ϕ[i], ϕ[i][1])
    end
  end

  # Vector of reactions and (i, j) coordinates of metabolite substrates/products
  subs = [findfirst(<(0), S[:,j]) for j in 1:size(S,2)]
  prods = [findfirst(>(0), S[:,j]) for j in 1:size(S,2)]
  reaction_pairs = tuple.(subs, prods)

  # Initialize and fill EFM matrix by converting EFM metabolite indices to rxns
  E = zeros(Int64, length(reaction_pairs), length(ϕ))
  for j in 1:length(ϕ)
    for i in 1:(length(ϕ[j])-1)
      k = (ϕ[j][i], ϕ[j][i+1])
      ii = findfirst([k == reaction_pairs[l] for l in 1:length(reaction_pairs)])
      E[ii, j] = 1
    end
  end
  return E
end

