## Helpers
# Int16 version of findall
function findall_int16(testf::Function, A)
  T = Int16
  gen = (first(p) for p in pairs(A) if testf(last(p)))
  isconcretetype(T) ? collect(T, gen) : collect(gen)
end

# Check if network is open or closed
function check_open_closed(S::Matrix{<:Int64})
  col_source = findall(==(1), vec(sum(S, dims=1)))
  col_sink = findall(==(-1), vec(sum(S, dims=1)))

  if isempty(col_source) && isempty(col_sink)
    return "closed"
  else
    sanitize_open_stoich(S)
    return "open"
  end
end

# Close network if open
function close_network(S::Matrix{<:Int64}, v::Vector{<:Real})
  # Error-checking for a unimolecular network with steady state fluxes
  sanitize_flux(v)
  sanitize_stoich_flux(S, v)

  # Indices of source/sink columns
  col_source = findall(==(1), vec(sum(S, dims=1)))
  col_sink = findall(==(-1), vec(sum(S, dims=1)))

  # Close off stoichiometry matrix
  rvec = zeros(Int64, size(S,2))
  rvec[col_source] .= -1
  rvec[col_sink] .= 1
  return vcat(rvec', S)
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
    prefix::Vector{Int16},
    t::Dict{#
      Vector{Int16},
      NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
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
  traverse_trie(Int16.([I]), t, T, A)

  return A, t
end

# Construct trie dictionary keyed by prefix with values equal to children
function trie(T::Matrix{<:Real}, I::Int64=1)

  # Initialize dictionary of prefix and children
  d = Dict{#
    Vector{Int16},
    Vector{Int16}
  }()
  function traverse_trie(#
    prefix::Vector{Int16},
    T::Matrix{<:Real},
    d::Dict{#
      Vector{Int16}, Vector{Int16}
    }
  )
    downstream = filter!(x -> x ∉ prefix, findall_int16(>(0), @view T[prefix[end],:]))
    d[prefix] = downstream
    for p in downstream
      traverse_trie([prefix; p], T, d)
    end
  end

  # Construct dictionary of prefix and children and add unique IDs
  traverse_trie(Int16.([I]), T, d)
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
    Vector{Int16},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
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
    upstream = findall(>(0), @view T′[prefixes_transformed[i][end],:])
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

  return res
end

