module Atomic
    using SparseArrays
    using ExtendableSparse
    using LinearAlgebra
    using LinearSolve

    using MarkovWeightedEFMs.CHMC: CHMCAtomicSummary, CHMCAtomicErrorSummary
    using MarkovWeightedEFMs.CHMC: findall_int16, enumerate_efms
    using MarkovWeightedEFMs.CHMC: solve_efm_probabilities
    using MarkovWeightedEFMs.CHMC: check_open_closed
    using PyCall # https://docs.alliancecan.ca/wiki/Julia for Python installation
    using PubChemCrawler: get_for_cids, parse_formula
    using MolecularGraph: molecularformula, smilestomol
    using Printf: @printf
    using Tables # for tmp_dir
    using CSV # for tmp_dir

    include("chmc-atomic.jl")
    export canonicalize_smiles
    export find_atomic_chmc_input_errors
    export print
    export correct_atomic_chmc_input_errors
    export correct_atomic_chmc_input_smiles
    export exchange_atomic_chmc_input_metabolites
    export get_source_metabolites
    export get_max_atoms
    export map_reaction_strings
    export precompute_atom_tracing_dictionary
    export steady_state_efm_distribution # imported from Standard
    export enumerate_atomic_efms
    export get_efm_metabolite_atom_indices
    export get_efm_reaction_atom_indices
    export chmc_to_mc_matrix
    #export compress

    include("io.jl")
    export export_atom_tracing_dictionary
    export import_atom_tracing_dictionary
end

