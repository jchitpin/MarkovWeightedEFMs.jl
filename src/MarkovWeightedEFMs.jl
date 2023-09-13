module MarkovWeightedEFMs

  using Reexport
  using Test

  # Algorithm
  include("chmc-algorithm/module.jl")
  @reexport using .CHMC

  # Plotting
  include("chmc-visualization/module.jl")
  @reexport using .Plots

  #=
  using Catalyst
  using Dates
  include("plot-catalyst-network.jl")
  export construct_catalyst_network
  export catalyst_to_graphviz
  =#
end # module
