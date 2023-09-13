module CHMC

  using Reexport: @reexport
  using SparseArrays
  using ExtendableSparse # unsure if needed
  using LinearAlgebra
  using IterativeSolvers
  using ArnoldiMethod: partialschur, LR

  # Core functions
  include("core.jl")
  export CHMCStandardSummary
  export CHMCAtomicSummary
  export CHMCAtomicErrorSummary
  export import_chmc

  include("standard/submodule.jl") # standard CHMC
  @reexport using .Standard

  #include("atomic/submodule.jl") # atomic CHMC
  #@reexport using .Atomic
end # module

