module MarkovWeightedEFMs

  using Test
  using Statistics: mean
  using GeometryBasics: Point
  using Makie: Combined, deregister_interaction!, autolimits!, limits!
  using Makie: register_interaction!, hidedecorations!, hidespines!, DataAspect
  using Graphs: SimpleDiGraph
  using NetworkLayout: Buchheim
  using GraphMakie: graphplot, NodeDragHandler

  include("plot-cycle-history-markov-chain.jl")
  export tree_plot

  using QuantEcon: MarkovChain, stationary_distributions

  include("cycle-history-markov-chain.jl")
  export steady_state_efm_distribution
  export stoich_to_transition
  export reshape_efm_matrix, reshape_efm_vector

  using MATLAB

  include("higher-order-generalization.jl")
  export steady_state_efm_distribution_higher

end # module
