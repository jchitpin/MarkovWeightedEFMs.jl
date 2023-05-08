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

  # Note: Rust compiler required for RXNMapper
  # https://docs.alliancecan.ca/wiki/Julia
  using PyCall
  using PubChemCrawler: get_for_cids, parse_formula
  using ArnoldiMethod: partialschur
  using MolecularGraph: molecularformula, smilestomol
  using DataFrames
  using CSV

  include("cycle-history-markov-chain.jl")
  export steady_state_efm_distribution
  export stoich_to_transition
  export enumerate_efms
  export reshape_efm_matrix, reshape_efm_vector



  include("atom-tracing.jl")
  export get_cid_info, translate_id, rxn_string


  #using MATLAB
  #include("higher-order-generalization.jl")
  #export steady_state_efm_distribution_higher

end # module
