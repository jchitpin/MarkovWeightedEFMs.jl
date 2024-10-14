module CHMC
    using Reexport: @reexport
    using ProgressMeter
    using ChunkSplitters: getchunk
    using SparseArrays
    using LinearAlgebra
    using LinearSolve 

    # Core functions
    include("core.jl")
    export CHMCStandardSummary
    export CHMCAtomicSummary
    export CHMCAtomicErrorSummary

    include("standard/submodule.jl") # standard CHMC
    @reexport using .Standard

    include("atomic/submodule.jl") # atomic CHMC
    @reexport using .Atomic
end # module

