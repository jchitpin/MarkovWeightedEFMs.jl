  module Standard
    using SparseArrays
    using MarkovWeightedEFMs.CHMC: CHMCStandardSummary
    using MarkovWeightedEFMs.CHMC: findall_int16, check_open_closed
    using MarkovWeightedEFMs.CHMC: trie_matrix, trie, enumerate_efms
    using MarkovWeightedEFMs.CHMC: solve_efm_probabilities
    include("chmc-standard.jl")

    export steady_state_efm_distribution
    export stoichiometry_to_transition_matrix
    export reshape_efm_matrix
    export reshape_efm_vector

    include("io.jl")
    export export_chmc
  end

