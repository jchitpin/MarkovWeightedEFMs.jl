module MarkovWeightedEFMs

  using Test
  using Statistics: mean
  using GeometryBasics: Point
  using Makie: Combined, deregister_interaction!, autolimits!, limits!
  using Makie: register_interaction!, hidedecorations!, hidespines!, DataAspect
  using Graphs: SimpleDiGraph
  using NetworkLayout: Buchheim
  using GraphMakie: graphplot, NodeDragHandler

  include("chmc-standard-visualization/plot-chmc.jl")
  export tree_plot

  # Note: Rust compiler required for RXNMapper
  # https://docs.alliancecan.ca/wiki/Julia for installing Python packages

  using ArnoldiMethod: partialschur, LR

  include("chmc-standard/chmc-standard.jl")
  export steady_state_efm_distribution
  export stoich_to_transition # for tree_plot
  export reshape_efm_matrix, reshape_efm_vector

  include("chmc-standard/io.jl")
  export export_chmc
  #export import_chmc # combined with chmc-atomic/io.jl function

  #using PyCall
  #using PubChemCrawler: get_for_cids, parse_formula
  #using MolecularGraph: molecularformula, smilestomol
  #using Printf: @printf
  #using Tables
  #using CSV

  #include("chmc-atomic/chmc-atomic.jl")
  #export steady_state_efm_distribution
  #export rxn_string, trace_rxn_string
  #export get_source_metabolites, get_max_atoms
  #export precompute_atom_tracing_dictionary
  #export count_atomic_chmc_state_space
  #export get_efm_met_sequence
  #export reconstruct_atomic_fluxes
  #export summary_atomic_chmc_inputs, print
  #export correct_atomic_chmc_inputs

  #include("chmc-atomic/io.jl")
  #export export_chmc
  #export import_chmc
  #export export_atom_tracing_dictionary
  #export import_atom_tracing_dictionary

  #using Catalyst
  #using Dates

  #include("chmc-atomic-visualization/plot-catalyst-network.jl")
  #export construct_catalyst_network
  #export catalyst_to_graphviz

  ##using DataFrames
  ##export get_cid_info, translate_id

end # module
