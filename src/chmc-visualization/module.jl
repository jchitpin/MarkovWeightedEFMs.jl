module Plots

  using Reexport: @reexport

  # For standard and atomic
  using SparseArrays
  using MarkovWeightedEFMs.CHMC: trie
  using Graphs: SimpleDiGraph
  using NetworkLayout: Buchheim
  using Makie: Combined, deregister_interaction!, autolimits!, limits!
  using Makie: register_interaction!, hidedecorations!, hidespines!, DataAspect
  using GraphMakie: graphplot, NodeDragHandler

  # For atomic
  using MarkovWeightedEFMs.CHMC.Atomic: get_srsi, get_prsi
  using MarkovWeightedEFMs.CHMC.Atomic: canonicalize_and_atom_map
  using MarkovWeightedEFMs.CHMC: CHMCAtomicSummary
  using RDKitMinimalLib: get_rxn, get_rxn_svg, get_svg, get_mol
  using Luxor: Drawing, origin, readsvg, placeimage, preview, finish, arrow
  using Luxor: background, sethue, fontsize, Point, text, svgstring
  using Luxor: image_as_matrix
  using Makie: RGBf, Figure, GridLayout, Axis, colgap!, rowgap!
  using Makie: Observable, image!
  using GraphMakie: graphplot!, EdgeClickHandler

  # Core functions
  include("core.jl")

  # Standard CHMC plotting
  include("standard/standard.jl")
  export plot_chmc

  # Atomic CHMC plotting
  include("atomic/atomic.jl")
  export plot_mapped_reaction
  export plot_atomic_chmc
end # module

