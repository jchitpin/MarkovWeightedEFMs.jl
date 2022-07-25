# Error checking helpers
function sanitize_transition(T::Matrix{<:Real}, I::Int64)
  @assert(#
    size(T,1) == size(T,2),
    "Transition probability matrix T must be square."
  )
  @assert(#
    I <= size(T,1),
    "Not a valid starting state."
  )
end


function trie(T::Matrix{<:Real}, I::Int64=1)

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

"""
    function tree_plot(#
        T::Matrix{<:Real},
        I::Int64=1;
        node_label_textsize::Real=15,
        edge_label_textsize::Real=12,
        arrow_shift::Real=0.85,
        x_pad::Real=0.75,
        y_pad::Real=0.75,
        tfactor::Real=0.15,
        tangents::Tuple{#
            Tuple{<:Real, <:Real},
            Tuple{<:Real, <:Real}
        } = ((1,0),(0,1)),
        show_all::Bool=false
)

Plot augmented prefix tree from transition probability matrix T rooted
on state/node I.

*T* is a right stochastic transition probability matrix. Alternatively, one
could specify a right generator matrix (with rows summing to zero).

*I* is a state in 1:size(T,1).

*node_label_textsize* is the text size of the node labels indexed from T.

*edge_label_textsize* is the text size of the edge labels taken from T.

*arrow_shift* is the percentage shift of the arrow head from src to dst.

*x_pad* is the left/right x coordinate padding of the plotting box.

*y_pad* is the up/down y coordinate padding of the plotting box.

*tfactor* scales the distance of the bezier control point relative to the
distance of the src and dst nodes.

*tangents* is the tangent of the src vertex and dst vertex.

*show_all = true* explicitly plots the upstream transition from all EFMs that
pass through the initial state/node *I*. By default, these arrows stemming
from the green nodes are omitted for visual clarity.
"""
function tree_plot(#
  T::Matrix{<:Real},
  I::Int64=1;
  node_label_textsize::Real=15,
  edge_label_textsize::Real=12,
  arrow_shift::Real=0.85,
  tfactor::Real=0.15,
  x_pad::Real=0.75,
  y_pad::Real=0.75,
  tangents::Tuple{Tuple{<:Real, <:Real}, Tuple{<:Real, <:Real}} = ((1,0),(0,1)),
  show_all::Bool=false
)

  # Error checking
  sanitize_transition(T, I)

  # New Markov chain transition probability matrix
  A, E, t = trie_adj_at_root(T, show_all, I)

  # Vector of node and edge colours
  node_color = color_nodes(t, T, I)
  edge_color = color_edges(A, E)
  node_labels = label_nodes(t)
  edge_labels = label_edges(A, E)

  # Tangents specify curved edges for looped EFM transitions
  tangents = tangent_edges(A, E, tangents)

  # Construct simple digraph
  G = SimpleDiGraph(A)

  # Construct Makie plot
  f, ax, p = graphplot(#
    G,
    nlabels = node_labels,
    elabels = edge_labels,
    labels_textsize = node_label_textsize,
    elabels_textsize = edge_label_textsize,
    arrow_show = true,
    arrow_shift = arrow_shift,
    node_color = node_color,
    tangents = tangents,
    tfactor = repeat([tfactor], length(tangents)),
    edge_color = edge_color
  )

  # Update current plot with custom tree layout
  x1, x2, y1, y2 = buchheim_layout(p, T, I)

  # Padding boundary box
  x1 -= x_pad
  x2 += x_pad
  y1 += y_pad
  y2 -= y_pad

  # Makie node interactivity
  function node_drag_action(state, idx, event, axis)
      p[:node_pos][][idx] = event.data
      p[:node_pos][] = p[:node_pos][]
  end
  ndrag = NodeDragHandler(node_drag_action)
  deregister_interaction!(ax, :rectanglezoom)
  register_interaction!(ax, :ndrag, ndrag)

  # Pretty-print Makie plot by removing axes/maintaining aspect/centering
  hidedecorations!(ax); hidespines!(ax)
  ax.aspect = DataAspect()
  autolimits!(ax)
  limits!(ax, x1, x2, y2, y1)

  # Don't forget to return the Makie plot
  f
end

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
      Vector{Int64},
      NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
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
    Vector{Int64},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
  },
)
  labels = [collect(keys(t))[i][end] for i in 1:length(t)]
  idx = parse.(Int64, [collect(values(t))[i].id[2:end] for i in 1:length(t)])
  return string.(labels[sortperm(idx)])
end

function color_nodes(
  t::Dict{#
    Vector{Int64},
    NamedTuple{(:id, :children), Tuple{String, Vector{Int64}}}
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

