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

  using PyCall
  using PubChemCrawler: get_for_cids, parse_formula
  using MolecularGraph: molecularformula, smilestomol
  using Tables
  using CSV

  include("chmc-atomic/chmc-atomic.jl")
  export steady_state_efm_distribution
  export rxn_string, trace_rxn_string
  export get_source_metabolites, get_max_atoms
  export precompute_atom_tracing_dictionary
  export count_atomic_chmc_state_space
  export get_efm_met_sequence
  export reconstruct_atomic_fluxes

  include("chmc-atomic/io.jl")
  export export_atomic_chmc
  export import_atomic_chmc
  export export_atom_tracing_dictionary
  export import_atom_tracing_dictionary



  ##using DataFrames
  ##export get_cid_info, translate_id

end # module
