using MarkovWeightedEFMs
using CairoMakie

# Figure 1b) in the manuscript; stoichiometry matrix, fluxes, efms,
# probability transition matrix
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
efms = [#
  [1, 2, 3, 5, 6, 1],
  [1, 2, 4, 5, 6, 1],
  [2, 3, 2],
  [3, 5, 3],
  [5, 7, 5],
  [2, 4, 5, 3, 2]
]
T = stoich_to_transition(S, v)

# Static backend for plotting
CairoMakie.activate!(type = "png")

# Plot and save
fig = tree_plot(T, 1)
save("prefix.png", fig)


