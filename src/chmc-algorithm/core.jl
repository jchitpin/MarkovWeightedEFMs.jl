# Structs
struct CHMCStandardSummary
  e::Vector{Vector{Int64}}
  p::Vector{Float64}
  w::Vector{Float64}
end

struct CHMCAtomicErrorSummary
  absolute_flux_error::Float64
  reactions_duplicated::Vector{Int64}
  reactions_with_zero_flux::Vector{Int64}
  reactions_with_negative_flux::Vector{Int64}
  reactions_with_non_integer_stoichiometries::Vector{Int64}
  reactions_unimolecular_source_stoichiometry_one::Vector{Int64}
  reactions_unimolecular_source_stoichiometry_not_one::Vector{Int64}
  reactions_multimolecular_source_stoichiometry_one::Vector{Int64}
  reactions_multimolecular_source_stoichiometry_not_one::Vector{Int64}
  reactions_unimolecular_sink_stoichiometry_one::Vector{Int64}
  reactions_unimolecular_sink_stoichiometry_not_one::Vector{Int64}
  reactions_multimolecular_sink_stoichiometry_one::Vector{Int64}
  reactions_multimolecular_sink_stoichiometry_not_one::Vector{Int64}
  reactions_empty::Vector{Int64}
  unused_metabolites::Vector{Int64}
end

struct CHMCAtomicSummary
  e::Vector{Vector{Int64}} # EFM sequences
  p::Vector{Float64} # EFM probabilities
  w::Vector{Float64} # EFM weights
  dmc::Dict{Int64, Tuple{Int16,Int16}} # Dictionary of MC states to (met idx/atom idx)
  i::Tuple{Int64, Int64, String} # Root MC state
  #T::Matrix{Float64} # CHMC transition matrix
  T::Union{Matrix{Float64},SparseMatrixCSC{Float64, Int64}} # CHMC transition matrix
  dchmc::Dict{# Dictionary of CHMC states based on MC states
    Vector{Int16},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
  }
  #R::Matrix{Int16} # CHMC reaction flux matrix
  R::Union{Matrix{Int16},SparseMatrixCSC{Int16, Int64}} # CHMC reaction flux matrix
end

# Int16 version of findall
function findall_int16(f::Function, A)
  T = Int16
  gen = (first(p) for p in pairs(A) if f(last(p)))
  isconcretetype(T) ? collect(T, gen) : collect(gen)
end

# Check if network is open or closed
function check_open_closed(S::Matrix{<:Int64})
  col_source = findall(==(1), vec(sum(S, dims=1)))
  col_sink = findall(==(-1), vec(sum(S, dims=1)))
  if isempty(col_source) && isempty(col_sink)
    return "closed"
  else
    col_source = findall(==(1), vec(sum(S, dims=1)))
    col_sink = findall(==(-1), vec(sum(S, dims=1)))
    bool = all([!isempty(col_source), !isempty(col_sink)])
    @assert(#
      bool,
      "An open network must contain at least one source and sink reaction."
    )
    return "open"
  end
end

# Construct cycle-history Markov chain
function trie_matrix(T::Matrix{<:Real}, I::Int64=1)
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
    T::Matrix{<:Real},
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
function trie_matrix(T::SparseMatrixCSC{Float64, Int64}, I::Int64=1)
  # Construct trie to get number of nodes
  t = trie(T, I)

  # Initialize new transition probability matrix
  A = spzeros(length(t), length(t))

  # Traverse tree and populate new transition probability matrix with trie nodes
  function traverse_trie(#
    prefix::Vector{Int16},
    t::Dict{#
      Vector{Int16},
      NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
    },
    T::SparseMatrixCSC{Float64, Int64},
    A::SparseMatrixCSC{Float64, Int64}
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
function trie(#
  T::Union{#
    Matrix{<:Real},
    SparseMatrixCSC{Int16, Int64},
    SparseMatrixCSC{Float64, Int64}},
  I::Int64=1
)
  # Initialize dictionary of prefix and children
  d = Dict{#
    Vector{Int16},
    Vector{Int16}
  }()
  function traverse_trie(#
    prefix::Vector{Int16},
    T::Union{Matrix{<:Real}, SparseMatrixCSC{Int16, Int64}, SparseMatrixCSC{Float64, Int64}},
    d::Dict{#
      Vector{Int16}, Vector{Int16}
    }
  )
    downstream = filter!(#
      x -> x ∉ prefix,
      findall_int16(>(0), @view T[prefix[end],:])
    )
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
function enumerate_efms(#
  T′::SparseMatrixCSC,
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
      @info("Remaining cycles: $(length(simple_cycles_original))")
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

# Import standard/atomic CHMC results
"""
    import_chmc(fname::String)

Import standard/atomic CHMC results from text file `fname`.

`fname` is the filename containing the CHMC results.
"""
function import_chmc(fname::String)
  # Load text file contents
  f = open(fname)
  lines = readlines(f)

  # Headers
  idx = findall(!isnothing, match.(r"## ", lines))

  # Extract data within each header and return CHMC results structure
  if length(idx) == 3 # standard CHMC
    efms =  eval.(Meta.parse.(lines[(idx[1]+1):(idx[2]-2)]))
    probs = parse.(Float64, lines[(idx[2]+1):(idx[3]-2)])
    weights = parse.(Float64, lines[(idx[3]+1):(length(lines)-1)])

    return CHMCStandardSummary(efms, probs, weights)
  elseif length(idx) == 8 # atomic CHMC
    init = eval(Meta.parse(lines[idx[1]+1]))
    aefms =  eval.(Meta.parse.(lines[(idx[2]+1):(idx[3]-2)]))
    probs = parse.(Float64, lines[(idx[3]+1):(idx[4]-2)])
    weights = parse.(Float64, lines[(idx[4]+1):(idx[5]-2)])
    d = split.(lines[(idx[5]+1):(idx[6]-2)], "\t")
    x = parse.(Int64, first.(d))
    y = eval.(Meta.parse.(last.(d)))
    dmc = Dict(zip(x, tuple.(Int16.(first.(y)), Int16.(last.(y)))))
    T = eval(Meta.parse(lines[idx[6]+1]))
    d = split.(lines[(idx[7]+1):(idx[8]-2)], "\t")
    x = eval.(Meta.parse.(first.(d)))
    y = eval.(Meta.parse.(last.(d)))
    dchmc = Dict(zip(x, y))
    R = eval(Meta.parse(lines[(idx[8]+1)]))
    return (e=aefms, p=probs, w=weights, dmc=dmc, i=init, T=T, dchmc=dchmc, R=R)
  else
    @warn("Cannot import the specified text file.")
  end
end

function solve_efm_probabilities(#
  T::Matrix{<:Real},
  T′::Matrix{<:Real},
  ϕ::Vector{NamedTuple{#
    (:EFM, :TransformedCycles),
    Tuple{Vector{Int64}, Vector{Vector{Int64}}}
  }},
  solver::Symbol
)
  π = Vector{Float64}(undef, size(T,1))
  if solver == :Direct
    π = [1; (LinearAlgebra.I - T′'[2:end,2:end]) \ Vector(T′'[2:end,1])]
    π = π / sum(π)
  elseif solver == :IterativeSolver_gmres
      x = gmres((LinearAlgebra.I - T′'[2:end,2:end]), Vector(T'[2:end,1]); log=true)
      @assert(string(x[2])[1:3] != "Not", "Failed to converge on eigenvector.")
      π = [1; x[1]]
      π = π ./ sum(π)
  elseif solver == :Arnoldi
  bool = true
    while bool
      decomp, history = partialschur(T′', nev=1, restarts=10000, which=LR())
      if history.converged == true
        if all(decomp.Q .>= 0)
          pii = vec((T′' * decomp.Q) / sum(T′' * decomp.Q))
          π = pii
          bool = false
        end
      end
    end
  end

  # Compute steady state edge probabilities and aggregate for each EFM
  p = Vector{Float64}(undef, length(ϕ))
  for i in 1:length(ϕ)
    for j in ϕ[i].TransformedCycles
      p[i] += π[j[1]] * T′[j[1], j[2]]
    end
  end
  return p
end
function solve_efm_probabilities(#
  T::SparseMatrixCSC,
  T′::SparseMatrixCSC,
  ϕ::Vector{NamedTuple{#
    (:EFM, :TransformedCycles),
    Tuple{Vector{Int64}, Vector{Vector{Int64}}}
  }},
  solver::Symbol
)
  π = Vector{Float64}(undef, size(T,1))
  if solver == :Direct
    π = [1; (LinearAlgebra.I - T′'[2:end,2:end]) \ Vector(T′'[2:end,1])]
    π = π / sum(π)
  elseif solver == :IterativeSolver_gmres
      x = gmres((LinearAlgebra.I - T′'[2:end,2:end]), Vector(T'[2:end,1]); log=true)
      @assert(string(x[2])[1:3] != "Not", "Failed to converge on eigenvector.")
      π = [1; x[1]]
      π = π ./ sum(π)
  elseif solver == :Arnoldi
  bool = true
    while bool
      decomp, history = partialschur(T′', nev=1, restarts=10000, which=LR())
      if history.converged == true
        if all(decomp.Q .>= 0)
          pii = vec((T′' * decomp.Q) / sum(T′' * decomp.Q))
          π = pii
          bool = false
        end
      end
    end
  end

  # Compute steady state edge probabilities and aggregate for each EFM
  p = Vector{Float64}(undef, length(ϕ))
  for i in 1:length(ϕ)
    for j in ϕ[i].TransformedCycles
      p[i] += π[j[1]] * T′[j[1], j[2]]
    end
  end
  return p
end

