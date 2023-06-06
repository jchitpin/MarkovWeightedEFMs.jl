function trie_adj_at_root(#
  T::Matrix{<:Real},
  show_all::Bool,
  I::Int64=1,
)
  # Construct trie to get number of nodes
  t = trie(T, I)

  # Initialize adjacency matrix
  A = zeros(length(t), length(t))
  E = Vector{Tuple{Int64, Int64}}()

  function traverse_trie(#
    prefix::Vector{Int64},
    t::Dict{#
      Vector{Int16},
      NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
    },
    T::Matrix{<:Real},
    A::Matrix{Float64},
    E::Vector{Tuple{Int64, Int64}},
    show_all::Bool
  )
    i = prefix[end] # old state
    ii = parse(Int64, t[prefix].id[2:end]) # new state

    # Fill adjacency matrix with non-children (upstream) transitions
    childs = t[prefix].children
    parents = findall(>(0), T[i,:])
    if show_all == false
      parents = parents[parents .!= prefix[1]]
    end
    parents = parents[parents .∉ Ref(childs)]

    for j in parents
      p = T[i,j]
      l = prefix[1:findfirst(x -> x == j, prefix)]
      jj = parse(Int64, t[l].id[2:end])
      A[ii,jj] = p
      push!(E, (ii, jj))
    end

    # Fill adjacency matrix with children transitions
    for j in childs
      p = T[i,j]
      jj = parse(Int64, t[[prefix; j]].id[2:end])
      A[ii,jj] = p
      traverse_trie([prefix; j], t, T, A, E, show_all)
    end
  end
  traverse_trie([I], t, T, A, E, show_all)
  return A, E, t
end

function label_nodes(#
  t::Dict{#
    Vector{Int16},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
  },
)
  labels = [collect(keys(t))[i][end] for i in 1:length(t)]
  idx = parse.(Int64, [collect(values(t))[i].id[2:end] for i in 1:length(t)])
  return string.(labels[sortperm(idx)])
end

function color_nodes(
  t::Dict{#
    Vector{Int16},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
  },
  T::Matrix{<:Real},
  I::Int64;
  colors::Vector{Symbol} = [:blue, :green, :black]
)

  # Get indices of root, "leaf", and internal nodes
  root = t[[I]].id
  root = replace(root, "X"  => "")

  leaf_id = findall(>(0), T[:,I])
  idx = findall([collect(keys(t))[i][end] for i in 1:length(t)] .∈ Ref(leaf_id))
  leaf = [collect(values(t))[i].id for i in idx]
  leaf = replace.(leaf, "X" => "")

  # Convert to int
  root = parse(Int64, root)
  leaf = parse.(Int64, leaf)

  # Fill colour node vector
  node_color = Vector{Symbol}()
  for n in 1:length(t)
    if n .∈ Ref(root)
      push!(node_color, colors[1])
    elseif n .∈ Ref(leaf)
      push!(node_color, colors[2])
    else
      push!(node_color, colors[3])
    end
  end
  return node_color
end

function color_edges(A::Matrix{Float64}, E::Vector{Tuple{Int64, Int64}})
  edge_color = Vector{Symbol}()
  cis = findall(>(0), A)
  cis = getindex.(cis, [1 2])
  sorted_cis = cis[sortperm(cis[:,1]),:]
  for e in eachrow(sorted_cis)
    if Tuple(e) ∈ E
      push!(edge_color, :red)
    else
      push!(edge_color, :gray)
    end
  end
  return edge_color
end

function label_edges(A::Matrix{Float64}, E::Vector{Tuple{Int64, Int64}})
  edge_label = Vector{String}()
  cis = findall(>(0), A)
  cis = getindex.(cis, [1 2])
  sorted_cis = cis[sortperm(cis[:,1]),:]
  for e in eachrow(sorted_cis)
    push!(edge_label, string(round(A[e[1],e[2]], digits=2)))
  end
  return edge_label
end

function tangent_edges(#
  A::Matrix{Float64},
  E::Vector{Tuple{Int64, Int64}},
  tangents::Tuple{Tuple{<:Real, <:Real}, Tuple{<:Real, <:Real}} = ((0,1),(0,1))
)
  tans = Dict{#
    Int64,
    Union{Nothing, Tuple{Tuple{<:Real, <:Real}, Tuple{<:Real, <:Real}}}
  }()
  cis = findall(>(0), A)
  cis = getindex.(cis, [1 2])
  sorted_cis = cis[sortperm(cis[:,1]),:]
  i = 0
  for e in eachrow(sorted_cis)
    i += 1
    if Tuple(e) ∈ E
      tans[i] = tangents
    else
      tans[i] = (nothing)
    end
  end
  return tans
end

function tree_layout(#
  p::Combined{graphplot, Tuple{SimpleDiGraph{Int64}}},
  t::Dict{#
    Vector{Int64},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}},
  I::Int64=1;
  x_pad::Real=0.75,
  y_pad::Real=0.75
)
  # Get tuple corresponding to each coordinate in graph
  coords = Vector{NamedTuple{(:nid, :row, :pid), Tuple{Int64, Int64, Int64}}}()
  function get_coords(#
    prefix::Vector{Int64},
    t::Dict{#
      Vector{Int64}, NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
    },
    coords::Vector{NamedTuple{(:nid, :row, :pid), Tuple{Int64, Int64, Int64}}}
  )
    nid = parse(Int64, replace(t[prefix].id, "X" => ""))
    pid = length(prefix) == 1 ? 0 : parse(Int64, replace(t[prefix[1:(end-1)]].id, "X" => ""))

    push!(coords, (nid=nid, row=length(prefix), pid=pid))
    for c in t[prefix].children
      get_coords([prefix; c], t, coords)
    end
  end
  get_coords([I], t, coords)

  # Determine number of nodes in each row
  rows = [coords[i].row for i in 1:length(coords)]
  uniqs = [count(==(element), rows) for element in unique(rows)]

  # Set root coordinates
  xl₀ = 0.0
  xr₀ = 0.0
  y₀ = y_pad

  # Iterating over each row
  for r in 1:length(uniqs)
    all_coords = coords[rows .== r]
    y₀ -= y_pad
    if length(all_coords) == 1
      #p[:node_pos][][all_coords[].nid] = Point2f(mean([xl₀, xr₀]), y₀)
      p[:node_pos][][all_coords[].nid] = Point{2, Float32}(mean([xl₀, xr₀]), y₀)
    else
      xl₀ -= x_pad
      xr₀ += x_pad

      # First coord goes leftmost and last coord goes rightmost
      #p[:node_pos][][all_coords[1].nid] = Point2f(xl₀, y₀)
      #p[:node_pos][][all_coords[end].nid] = Point2f(xr₀, y₀)
      p[:node_pos][][all_coords[1].nid] = Point{2, Float32}(xl₀, y₀)
      p[:node_pos][][all_coords[end].nid] = Point{2, Float32}(xr₀, y₀)

      # Remaining coords spread in between
      interval = (-xl₀ + xr₀) / (length(all_coords) - 1)
      int_idx = 0
      for j in all_coords[2:(end-1)]
        int_idx += 1
        #p[:node_pos][][j.nid] = Point2f(xl₀ + (int_idx * interval), y₀)
        p[:node_pos][][j.nid] = Point{2, Float32}(xl₀ + (int_idx * interval), y₀)
      end
    end
  end
  # Update Makie
  notify(p[:node_pos])

  return xl₀, xr₀, 0, y₀
end

function buchheim_layout(
  p::Combined{graphplot, Tuple{SimpleDiGraph{Int64}}},
  T::Matrix{<:Real},
  I::Int64
)

  A, E, t = trie_adj_at_root(T, false, I)

  # Remove reversible edges from transition matrix A
  for e in E
    A[e[1], e[2]] = 0
  end

  G = SimpleDiGraph(A)
  f2, ax2, p2 = graphplot(G, layout = Buchheim())

  for j in 1:length(p[:node_pos][])
    p[:node_pos][][j] = p2[:node_pos][][j]
  end

  # Update Makie
  notify(p[:node_pos])

  xl₀ = minimum([p[:node_pos][][i][1] for i in 1:length(p[:node_pos][])])
  xr₀ = maximum([p[:node_pos][][i][1] for i in 1:length(p[:node_pos][])])
  yu₀ = maximum([p[:node_pos][][i][2] for i in 1:length(p[:node_pos][])])
  yd₀ = minimum([p[:node_pos][][i][2] for i in 1:length(p[:node_pos][])])
  return xl₀, xr₀, yu₀, yd₀
end


