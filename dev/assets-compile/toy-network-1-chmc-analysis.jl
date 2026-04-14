using MarkovWeightedEFMs
using CairoMakie

# Toy network 1
S = [#
    -1  0  0  0  0  0  0  0  0  0  1
     1 -1  1 -1  0  0  0  0  0  0  0
     0  1 -1  0 -1  1  0  0  0  0  0
     0  0  0  1  0  0 -1  0  0  0  0
     0  0  0  0  1 -1  1 -1 -1  1  0
     0  0  0  0  0  0  0  1  0  0 -1
     0  0  0  0  0  0  0  0  1 -1  0
]
v = [3, 2, 1, 2, 3, 2, 2, 3, 1, 1, 3]
res = steady_state_efm_distribution(S, v)
efms = [#
  [3, 2, 3],
  [3, 2, 4, 5, 3],
  [3, 5, 3],
  [6, 1, 2, 3, 5, 6],
  [7, 5, 7],
  [6, 1, 2, 4, 5, 6]
]
res.e == efms # true

T = stoichiometry_to_transition_matrix(S, v)

# Static backend for plotting
CairoMakie.activate!(type = "png")

# Plot and save
fig = plot_chmc(T, 1)
save("toy-network-1-chmc-makie.png", fig)

