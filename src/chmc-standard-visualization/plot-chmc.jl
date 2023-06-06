# Main function for plotting
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

*T* is a right stochastic transition probability matrix with rows summing to one.

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

