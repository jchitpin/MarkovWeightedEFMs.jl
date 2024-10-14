using MarkovWeightedEFMs
using ProgressMeter
using SparseArrays
using BenchmarkTools
using Random
rng = MersenneTwister(0)

#T = sprand(rng, 15, 15, 0.3) # 20,349 CHMC states
T = sprand(rng, 15, 15, 0.45) # 1,117,312 CHMC states
#T = sprand(rng, 16, 16, 0.4) # 4,350,843 CHMC states
T = sprand(rng, 17, 17, 0.4) # 353,249,957 CHMC states

[T[i,i] = 0 for i in 1:size(T,1)]
T = T ./ sum(T, dims = 2)
I = Int16(1)
#T = Matrix(T)
verbose = true

@elapsed d = MarkovWeightedEFMs.CHMC.trie(T)
@elapsed T′, d = MarkovWeightedEFMs.CHMC.trie_matrix(T, 1, verbose = true)
MarkovWeightedEFMs.CHMC.enumerate_efms(T′, d; verbose = true)
MarkovWeightedEFMs.CHMC.enumerate_efms(T′, d; verbose = false)

module test
import Random
    using Random
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60.0
@benchmark MarkovWeightedEFMs.CHMC.trie_matrix(T, 1)

@btime MarkovWeightedEFMs.CHMC.trie_matrix(T, 1)
@elapsed MarkovWeightedEFMs.CHMC.trie_matrix(T, 1)


S = [#
 -1  0  0  0  0  0  0  0  0  0  1
  1 -1  1 -1  0  0  0  0  0  0  0
  0  1 -1  0 -1  1  0  0  0  0  0
  0  0  0  1  0  0 -1  0  0  0  0
  0  0  0  0  1 -1  1 -1  1 -1  0
  0  0  0  0  0  0  0  0  0  1 -1
  0  0  0  0  0  0  0  1 -1  0  0
];
v = [3, 2, 1, 2, 3, 2, 2, 1, 1, 3, 3];
res = steady_state_efm_distribution(S, v);

res = steady_state_efm_distribution(T);

