#
using BenchmarkTools
using MarkovWeightedEFMs
using QuantEcon
using MATLAB
m = "/home/jchitpin/Documents/PhD/Code/MATLAB/efmtool/"
#

#= Toy network 39 (inverse of network 38)
S = [#
     1 -1 -1  0  0  0  0  0  0  0
     0  1  0 -1  1  0  0  0  0  0
     0  0  1 -1  1 -1  0  0  0  0
     0  0  0  1 -1  0 -1  0  0  0
     0  0  0  1 -1  0  0 -1  0  0
     0  0  0  1 -1  0  0  0 -1  0
     0  0  0  0  0  1  0  0  0 -1
]
v = [3, 1, 2, 2, 1, 1, 1, 1, 1, 1]


using MarkovWeightedEFMs
using GLMakie
T = [#
  0   1/3 2/3 0   0   0   0
  0   0   0   1/3 1/3 1/3 0
  0   0   0   2/9 2/9 2/9 1/3
  1/2 1/4 1/4 0   0   0   0
  1/2 1/4 1/4 0   0   0   0
  1/2 1/4 1/4 0   0   0   0
  1   0   0   0   0   0   0
]
tree_plot(T,1,show_all=false)
=#

#= Toy network (TP) this one leads to rxn_map and efm_map containing multiple elements
S = [#
     1 -1 -1 -1  0  0  0  0  0  0  0
     0  1  0  0 -1  1  0  0  0  0  0
     0  0  1  0 -1  1 -1  1  0  0  0
     0  0  0  1  0  0 -1  1  0  0  0
     0  0  0  0  1 -1  1 -1 -1  0  0
     0  0  0  0  1 -1  1 -1  0 -1  0
     0  0  0  0  1 -1  1 -1  0  0 -1
]
v = [4, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2]
S * v

res = steady_state_efm_distribution_higher(S, v, m, 1)
=#

# Error-checking helpers
function sanitize_stoich_higher(S::Matrix{<:Real})
  subs = length.([findall(<(0), S[:,i]) for i in 1:size(S,2)])
  prods = length.([findall(>(0), S[:,i]) for i in 1:size(S,2)])
  [@assert(#
    !(all(S[:,j] .== 0)),
    "Reaction column cannot be empty (all zeroes)."
  ) for j in 1:size(S,2)]
  external_cols = findall(==(1), [sum(S[:,j] .!= 0) for j in 1:size(S,2)])
  if !isempty(external_cols)
    src_cols = external_cols[findall([sum(S[:,j]) > 0 for j in external_cols])]
    sink_cols = external_cols[findall([sum(S[:,j]) < 0 for j in external_cols])]
    msg = "An open network must contain both source and sink reactions."
    @assert(!isempty(src_cols), msg)
    @assert(!isempty(sink_cols), msg)
  end
end
function sanitize_flux_higher(v::Vector{<:Real})
  @assert(all(v .>= 0), "Fluxes must be ≥ 0.")
end
function sanitize_stoich_flux_higher(S::Matrix{<:Real}, v::Vector{<:Real})
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
function sanitize_efms_higher(ϕ::Matrix{Float64}, S::Matrix{<:Real})
  @assert(all(ϕ .∈ Ref([-1, 0, 1])), "Only binary EFMs accepted.")
  @assert(#
    size(ϕ,1) == size(S,2),
    "Rows in E must equal columns in S (same number of reactions)."
  )
end
function sanitize_efms_higher(ϕ::Vector{Vector{Int64}}, S::Matrix{<:Real})
  @assert(#
    maximum(vcat(ϕ...)) <= size(S,1),
    "EFM indices must be 1:size(S,1)."
  )
end
function sanitize_transition_matrix_higher(T::Matrix{<:Real})
  @assert(size(T,1) == size(T,2), "T must be a square matrix.")
  @assert(#
    all([round(sum(T[i,:]), digits = 5) == 1 for i in 1:size(T,1)]),
    join([#
      "T[i,:] must be a right stochastic matrix with rows summing to one. ",
      "Check that the stoichiometry matrix is fully-connected and ",
      "includes non-zero fluxes connected to each metabolite. ",
      "Also ensure there are no duplicate reactions between metabolites."
    ])
  )
end
function sanitize_initial_higher(I::Int64, S::Matrix{<:Real})
  @assert(0 < I <= size(S,1), "I must belong in size(S,1).")
end

# Main function
function steady_state_efm_distribution_higher(#
  S::Matrix{<:Real},
  v::Vector{<:Real},
  m::String,
  I::Int64=1;
  test::Vector{<:Real}
)
  # Error-checking inputs
  sanitize_stoich_higher(S)
  sanitize_flux_higher(v)
  sanitize_stoich_flux_higher(S, v)
  sanitize_initial_higher(I, S)

  # Enumerate EFMs in original network
  ϕ = call_efmtool(Float64.(S), m) # S must be Float64 for MATLAB!

  # Transform into simplified unimolecular network with aggregated fluxes
  S′, v′ = transform_network(S, v)

  # Convert to closed, simplified network
  csnet = close_network(S, S′, v, v′) # 3 fields: internal/external rxns, fluxes

  # Generate transition probability matrix
  if isnothing(csnet.R)
    T = stoich_to_transition_higher(csnet.S′′, csnet.v′′)
  else
    T = stoich_to_transition_higher(hcat(csnet.S′′, csnet.R), csnet.v′′)
  end

  # Construct cycle-history transition probability matrix and prefix dictionary
  T′2, d2 = trie_matrix_higher(T, I)

  # Enumerate all EFMs/simple cycles from the cycle-history matrix/prefix
  ϕ′′2 = enumerate_efms(T′2, d2)
  e2 = [ϕ′′2[i].EFM for i in 1:length(ϕ′′2)]

  # Map reactions in closed, simplified network to those in original network
  if isnothing(csnet.R)
    rxn_map2 = map_rxns(S, csnet.S′′)
    A = reshape_efm_vector_higher(e2, csnet.S′′)
  else
    rxn_map2 = map_rxns(S, csnet.S′′, csnet.R)
    A = reshape_efm_vector_higher(e2, hcat(csnet.S′′, csnet.R))
  end

  # Map EFMs in closed, simplified network to those in original network
  efm_map2 = map_efms2(ϕ, A, rxn_map2) # ϕ′′ is now e in matrix form

  # Compute the steady state probabilities of closed, simplified network
  mc2 = MarkovChain(T′2)
  π′′2 = stationary_distributions(mc2)[1]

  # Compute steady state edge probabilities and aggregate for each EFM
  p2 = Vector{Float64}(undef, length(ϕ′′2))
  for i in 1:length(ϕ′′2)
    for j in ϕ′′2[i].TransformedCycles
      p2[i] += π′′2[j[1]] * T′2[j[1], j[2]]
    end
  end
  p2 = p2 ./ sum(p2)

  # Sum unimolecular EFMs belonging to each original EFM
  q = [#
    sum(p[findall([any(efm_map[i] .∈ Ref(j)) for i in 1:length(efm_map)])])
    for j in 1:size(ϕ,2)
  ]

  # p[29] = 0.19355 (index 1)
  # sum(p[[7, 14, 21, 24, 25, 33]]) = 0.38710 (index 2)
  # sum(p[[2, 4, 8, 10, 17, 26]]) = 0.29032 (index 3)

  # Divide EFM weights by the number of particles required to complete each EFM
  if isnothing(csnet.R)
    num_particles = test
    #@error("Not sure what to do for closed networks at the moment.")
  else
    external_cols = findall(==(1), [sum(S[:,j] .!= 0) for j in 1:size(S,2)])
    src_cols = external_cols[findall([sum(S[:,j]) > 0 for j in external_cols])]
    num_particles = vec(sum(ϕ[src_cols,:], dims=1))

    # What about internal loop particle counts?
    for j in size(ϕ,2)
      if num_particles[j] == 0
        rxns = S[:,findall(>(0), ϕ[:,j])]
        met_neg = findall(<(0), rxns[:,1])
        met_pos = findall(>(0), rxns[:,1])
        idx_neg = findfirst.(==(met_neg[1]), collect(keys(d)))
        idx_neg = maximum(idx_neg[.!isnothing.(idx_neg)])
        idx_pos = findfirst.(==(met_pos[1]), collect(keys(d)))
        idx_pos = maximum(idx_pos[.!isnothing.(idx_pos)])
        if idx_neg < idx_pos
          num_particles[j] = length(met_neg)
        else
          num_particles[j] = length(met_pos)
        end
        #num_particles[j] = maximum([id_pos; id_neg])
      end
    end
    #num_particles[num_particles .== 0] .= 1
  end
  qq = q ./ num_particles
  qq = qq ./ sum(qq)

  # Compute EFM weights from proportionality constant
  # (Now the sum of reaction coefficients in each EFM rather than length)
  c = [sum(ϕ[:,i]) for i in 1:size(ϕ,2)]
  α = sum(v) / sum(qq .* c)
  w = qq * α

  return (e=ϕ, p=qq, w=w)
end

# Main function helpers
function call_efmtool(S::Matrix{Float64}, efmtool_dir::String)
  stoich = mxarray(S)
  reversibilities = mxarray(zeros(size(S,2)))
  mat"cd($efmtool_dir)"
  mat"model.stoich = $stoich"
  mat"model.reversibilities = $reversibilities"
  mat"$res = CalculateFluxModes(model)"
  return res["efms"]
end



function transform_network(S::Matrix{<:Real}, v::Vector{<:Real})
  # Simplify higher order reactions into unimolecular ones
  # Distribute fluxes uniformly across simplified reactions
  cols = Vector{Vector{Float64}}()
  vs = Vector{Float64}()
  for j in 1:size(S,2)
    b1 = length(findall(==(-1), S[:,j])) > 1
    b2 = length(findall(==(+1), S[:,j])) > 1
    if b1 + b2 > 0
      idx_neg = findall(<(0), S[:,j])
      idx_pos = findall(>(0), S[:,j])
      if isempty(idx_neg)
        idx_neg = zeros(Int64, length(idx_pos))
      end
      if isempty(idx_pos)
        idx_pos = zeros(Int64, length(idx_neg))
      end
      combos = unique(vcat(collect(Iterators.product(idx_neg, idx_pos))...))
      for c in combos
        temp = zeros(Float64, size(S, 1))
        if c[1] != 0.0
          temp[c[1]] = -1.0
        end
        if c[2] != 0.0
          temp[c[2]] = 1.0
        end
        push!(cols, temp)
        push!(vs, (v[j] * length(idx_neg)) / length(combos))
      end
    else
      push!(cols, S[:,j])
      push!(vs, v[j])
    end
  end

  # Remove non-unique simplified reaction but aggregate fluxes
  i = 0
  while i <= length(cols)
    i += 1
    if i > length(cols)
      break
    end
    rm_idx = findall(cols .== Ref(cols[i]))[2:end]
    cols = cols[setdiff(1:end, rm_idx)]
    vs[i] = vs[i] * (1 + length(rm_idx))
    vs = vs[setdiff(1:end, rm_idx)]
  end

  return hcat(cols...), vs
end
function check_closed(S::Matrix{<:Real})
  if isempty(findall(==(1), [sum(S[:,j] .!= 0) for j in 1:size(S,2)]))
    @info("Closed network detected.")
    return true
  else
    @info("Open network detected.")
    return false
  end
end
function close_network(#
  S::Matrix{<:Real},
  S′::Matrix{<:Real},
  v::Vector{<:Real},
  v′::Vector{<:Real}
)
  status = check_closed(S)
  if status == true
    return (S′′=S′,R=nothing, v′′=v′)
  else
    # Internal reaction stoichiometry matrix (with pseudo metabolite,
    # pseudo reactions, indices of external reactions in S, and flux vector
    # for hcat(S′′, R)
    S′′, R, v′′ = close_network_inner(S′, v′)

    return (S′′=S′′,R=R, v′′=v′′)
  end
end
function close_network_inner(#
  S′::Matrix{<:Real},
  v′::Vector{<:Real}
)
  # Stoichiometry matrix with no source and sink reactions and pseudo metabolite
  src_cols = findall([all(S′[:,j] .>= 0) for j in 1:size(S′,2)])
  sink_cols = findall([all(S′[:,j] .<= 0) for j in 1:size(S′,2)])
  T = S′[:, setdiff(1:size(S′,2), [src_cols; sink_cols])]
  T = vcat(T, zeros(1, size(T,2)))
  src_rows = [findfirst(S′[:,j] .> 0) for j in src_cols]
  sink_rows = [findfirst(S′[:,j] .< 0) for j in sink_cols]

  # Create pseudo reactions linking sink to pseudo metabolite
  pseudo_sink = zeros(Float64, size(T,1), length(sink_rows))
  [pseudo_sink[[sink_rows[i],end],i] = [-1, 1] for i in 1:length(sink_rows)]

  # Create pseudo reactions linking pseudo metabolite to source
  pseudo_src = zeros(Float64, size(T,1), length(src_rows))
  [pseudo_src[[end,src_rows[i]],i] = [-1, 1] for i in 1:length(src_rows)]

  # Bind pseudo reactions
  pseudo_rxns = hcat(pseudo_sink, pseudo_src)

  # Internal fluxes first then source/sink pseudo reaction fluxes
  v′′ = v′[setdiff(1:size(S′,2), [src_cols; sink_cols])]
  pseudo_fluxes = v′[[sink_cols; src_cols]]
  v′′ = vcat(v′′, pseudo_fluxes)

  # Bind pseudo reactions to stoichiometry matrix
  return T, pseudo_rxns, v′′
end
function trie_matrix_higher(#
  T::Matrix{<:AbstractFloat},
  I::Int64=1,
)
# Error checking I
  #sanitize_initial(I, size(T, 1))

  # Construct trie to get number of nodes
  t = trie_higher(T, I)

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
    parents_new = copy(parents)
    cc = Tuple{Vector{Int64}, Vector{Int64}}

    for j in parents_new
      if typeof(cc) != DataType && any(Ref(j) .∈ cc[1])
        r = sum(T[i,cc[1]])
      else
        r = T[i,j]
      end
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

  # Re-normalize transition matrix rows to sum to probability of one
  A = A ./ sum(A, dims=2)

  return A, t
end

function trie_higher(#
  T::Matrix{<:Real},
  I::Int64=1
)
  # Initialize dictionary of prefix and children
  d = Dict{#
    Vector{Int64},
    Vector{Int64}
  }()
  function traverse_trie(#
    prefix::Vector{Int64},
    T::Matrix{<:Real},
    d::Dict{#
      Vector{Int64}, Vector{Int64}
    }
  )
    downstream = filter!(x -> x ∉ prefix, findall(>(0), T[prefix[end],:]))
    d[prefix] = downstream
    for p in downstream
      traverse_trie([prefix; p], T, d)
    end
  end

  # Construct dictionary of prefix and children and add unique IDs
  traverse_trie([I], T, d)
  v1 = ["X" * string(i) for i in 1:length(keys(d))]
  v2 = collect(values(d))

  # Ensure root node/prefix is variable "X1"
  root_idx = findfirst([collect(keys(d))[i][end] == I for i in 1:length(d)])
  switch_idx = findfirst(v1 .== "X1")
  v1[switch_idx] = v1[root_idx]
  v1[root_idx] = "X1"

  vs = [(id=v1[i], children=v2[i]) for i in 1:length(v1)]

  return Dict(zip(keys(d), vs))
end

# Enumerate EFMs via counting simple cycles
function enumerate_efms(#
  T′::Matrix{Float64},
  d::Dict{#
    Vector{Int64},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
  }
)

  # Prefixes with cycle-history Markov chain states
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

  # For each CHMC prefix, check if the last state transitions back to a CHMC
  # state already contained in the prefix. (Transitioning "up" the prefix tree)
  # These simple cycles represent the EFMs
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

  # Aggregate CHMC simple cycles over each EFM
  res = NamedTuple{#
    (:EFM, :TransformedCycles),
    Tuple{Vector{Int64}, Vector{Vector{Int64}}}
  }[]
  simple_cycles_original_2 = [i[2:end] for i in simple_cycles_original]
  while !isempty(simple_cycles_original)
    # Find all simple cycles corresponding to the same EFM
    tmp = simple_cycles_original_2[1]
    tmp = [tmp[[i:length(tmp); collect(1:(i-1))]] for i in 1:length(tmp)]
    ids = vcat(#
      [findall(simple_cycles_original_2 .== Ref(tmp[i])) for i in 1:length(tmp)]
    ...)
    push!(#
      res,
      (#
        EFM=simple_cycles_original[1],
        TransformedCycles=simple_cycles_transformed[ids]
      )
    )
    splice!(simple_cycles_original, sort(ids))
    splice!(simple_cycles_original_2, sort(ids))
    splice!(simple_cycles_transformed, sort(ids))
  end

  # Aggregate CHMC simple cycles over each EFM
  #res = NamedTuple{#
    #(:EFM, :TransformedCycles),
    #Tuple{Vector{Int64}, Vector{Vector{Int64}}}
  #}[]
  #while !isempty(simple_cycles_original)
    #idx = findall(#
      #in(Ref(Set(simple_cycles_original[1][2:end]))),
      #[Set(i[2:end]) for i in simple_cycles_original]
    #)
    #push!(res, (EFM=simple_cycles_original[1], TransformedCycles=simple_cycles_transformed[idx]))
    #splice!(simple_cycles_original, idx)
    #splice!(simple_cycles_transformed, idx)
  #end
  return res
end

function map_rxns(#
  S::Matrix{<:Real},
  S′′::Matrix{<:Real},
  R::Matrix{<:Real}
)
  # Match unimolecular to original EFMs based on (all) matching reactions
  rxn_map = Vector{Vector{Int64}}(undef, size(hcat(S′′, R),2))
  for k in 1:size(S′′,2)
    idx = Ref(S′′[S′′[:,k] .!= 0,k]) .==
      [S[findall(S′′[:,k] .!= 0),j] for j in 1:size(S,2)]
    rxn_map[k] = findall(idx)
  end

  # Pseudo reaction (source/sink) column indices
  src_cols = findall([all(S[:,j] .>= 0) for j in 1:size(S,2)])
  sink_cols = findall([all(S[:,j] .<= 0) for j in 1:size(S,2)])
  external_cols = [src_cols; sink_cols]

  # Match pseudo reactions to external reactions if open network
  for k in 1:size(R,2)
    tmp = Vector{Int64}()
    for j in external_cols
      if !isempty(findall(<(0), R[1:(end-1),k]))
        if any(findall(<(0), R[1:(end-1),k]) .∈ findall(<(0), S[:,j]))
          push!(tmp, j)
        end
      end
      if !isempty(findall(>(0), R[1:(end-1),k]))
        if any(findall(>(0), R[1:(end-1),k]) .∈ findall(>(0), S[:,j]))
          push!(tmp, j)
        end
      end
    end
    rxn_map[size(S′′,2) + k] = tmp
  end

  # Each reaction in simplified network must match to one or more in original
  @assert(#
    all([isassigned(rxn_map, i) for i in 1:length(rxn_map)]),
    "A reaction in the simplified network failed to match with one in the original"
  )
  return rxn_map
end
function map_rxns(#
  S::Matrix{<:Real},
  S′′::Matrix{<:Real}
)
  # Match unimolecular to original EFMs based on (all) matching reactions
  rxn_map = Vector{Vector{Int64}}(undef, size(S′′,2))
  for k in 1:size(S′′,2)
    idx = Ref(S′′[S′′[:,k] .!= 0,k]) .==
      [S[findall(S′′[:,k] .!= 0),j] for j in 1:size(S,2)]
    rxn_map[k] = findall(idx)
  end

  # Pseudo reaction (source/sink) column indices
  src_cols = findall([all(S[:,j] .>= 0) for j in 1:size(S,2)])
  sink_cols = findall([all(S[:,j] .<= 0) for j in 1:size(S,2)])
  external_cols = [src_cols; sink_cols]

  # Each reaction in simplified network must match to one or more in original
  @assert(#
    all([isassigned(rxn_map, i) for i in 1:length(rxn_map)]),
    "A reaction in the simplified network failed to match with one in the original"
  )
  return rxn_map
end


function map_efms2(ϕ::Matrix{<:Real}, A::Matrix{<:Real}, rxn_map::Vector{Vector{Int64}})

  # Initialize output
  #efm_map = Vector{Vector{Int64}}(undef, size(A,2))

  # Original reaction indices for each EFM
  F = [findall(>(0), ϕ[:,j]) for j in 1:size(ϕ,2)]

  # Simplified reaction indices for each EFM
  B = [findall(==(1), A[:,j]) for j in 1:size(A,2)]

  # Enumerate mapped reaction vectors for each simplified EFM
  mapr = [combinations_of_one_element_set(rxn_map[i]) for i in B]

  # Only keep reaction vectors with unique elements
  unique_nnarray!(mapr)

  # Output indices of original EFMs mapping to each simplified EFM
  return efm_indices(mapr, F)
end

function efm_indices(mapr::Vector{Vector{Vector{Int64}}}, F::Vector{Vector{Int64}})
  # Initialize output
  #efm_map = Vector{Vector{Int64}}(undef, length(mapr))
  efm_map = fill(Int64[], length(mapr))
  for i in 1:length(mapr)
    tmp = Vector{Int64}()
    if !isempty(mapr[i])
      for j in 1:length(mapr[i])
        if any(issubset.(Ref(mapr[i][j]), F))
          push!(tmp, findfirst(issubset.(Ref(mapr[i][j]), F)))
        end
      end
    end
    if !isempty(tmp)
      efm_map[i] = tmp
    end
  end
  return efm_map
end

function unique_nnarray!(x::Vector{Vector{Vector{Int64}}})
  for i in 1:length(x)
    j = 1
    while j <= length(x[i])
      c = false
      for k in 2:length(x[i][j])
        sort!(x[i][j])
        if isequal(x[i][j][k], x[i][j][k-1])
          c = true
          break
        end
      end
      if c == true
        deleteat!(x[i], j)
      else
        j += 1
      end
    end
  end
end

# Recursive algorithm to generate combinations of n elements from n sets of m elements
function combinations_of_one_element_set(x::Vector{Vector{Int64}})
  if length(x) == 1
    return [y for y in x[1]]
  end
  x′ = Vector{Vector{Int64}}()
  for x_el in x[1]
    for next_x in combinations_of_one_element_set(x[2:end])
      push!(x′, [[x_el]; next_x])
    end
  end
  return x′
end

function stoich_to_transition_higher(S::Matrix{<:Real}, v::Vector{<:Real})
  # Error-checking for a fully-connected, unimolecular network at steady state
  sanitize_stoich_higher(S)
  sanitize_flux_higher(v)

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
  sanitize_transition_matrix_higher(T)

  return T
end
function reshape_efm_vector_higher(#
  ϕ::Vector{Vector{Int64}},
  S::Matrix{<:Real}
)
  ## Error checking
  sanitize_stoich_higher(S)
  sanitize_efms_higher(ϕ, S)

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
  E = zeros(Float64, length(reaction_pairs), length(ϕ))
  for j in 1:length(ϕ)
    for i in 1:(length(ϕ[j])-1)
      k = (ϕ[j][i], ϕ[j][i+1])
      ii = findfirst([k == reaction_pairs[l] for l in 1:length(reaction_pairs)])
      E[ii, j] = 1.0
    end
  end
  return E
end



