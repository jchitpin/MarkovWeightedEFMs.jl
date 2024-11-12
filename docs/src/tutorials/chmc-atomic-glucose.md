# ACHMC (for AEFMs; quickstart)

This section shows how to use the wrapper function in
MarkovWeightedEFMs.jl to enumerate and assign AEFMs weights in a simple,
multispecies reaction network.

![Toy multispecies network](../assets/toy-network-1-achmc.png)

```@setup required
# Bodge until I figure out how to install Python and RXMapper on Git workflows
using MarkovWeightedEFMs
using SparseArrays

S = [#
  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 # Glc
  0 -1  0  0  0 -1  0  0  0  0  0  0  1  0  0  0 # ATP
  0  1 -1  0 -1  0  0  0  0  0  0  0  0  0  0  0 # G6P
  0  1  0  0  0  1  0  0  0  0  0  0  0 -1  0  0 # ADP
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 # 6PG
  0  0  0  0  1 -1  1  0  0  0  0  0  0  0  0  0 # F6P
  0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1  0 # Pi
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1 # H2O
  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0 # FDP
  0  0  0  0  0  0  0  1 -1  1 -1  0  0  0  0  0 # G3P
  0  0  0  0  0  0  0  1  1 -1  0 -1  0  0  0  0 # DHAP
]

v = [10, 10, 3, 3, 7, 8, 1, 7, 1, 1, 7, 7, 18, 18, 1, 1]

mets = [#
  "Glc",
  "ATP",
  "G6P",
  "ADP",
  "6PG",
  "F6P",
  "Pi",
  "H2O",
  "FDP",
  "G3P",
  "DHAP"
]

rxns = [#
  "Source Glc",
  "Hexokinase",
  "G6P dehydrogenase",
  "Sink 6PG",
  "Phosphoglucose isomerase",
  "6-phosphofructo-1-kinase",
  "Fructose 1,6-bisphosphatase",
  "Fructose-bisphosphate aldolase",
  "Triose phosphate isomerase",
  "Triose phosphate isomerase",
  "Sink G3P",
  "Sink DHAP",
  "Source ATP",
  "Sink ADP",
  "Sink Pi",
  "Source H2O"
]

smiles = [#
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)OP(=O)(O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
  "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O",
  "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)OP(=O)(O)O",
  "[O-]P(=O)([O-])[O-]",
  "O",
  "C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O",
  "C([C@H](C=O)O)OP(=O)(O)O",
  "C(C(=O)COP(=O)(O)O)O"
]

atom = :C # carbon atom type for AEFMs

mdl, atom_info, logs = preprocess_all_for_atomic_chmc(S, v, mets, rxns, smiles, atom)

I = (1, 1, atom)
ks = Int16.([5, 4, 6, 7, 2, 8, 3, 1])
vs = [(Int16(6), Int16(6)), (Int16(0), Int16(0)), (Int16(9), Int16(4)), (Int16(11), Int16(2)), (Int16(3), Int16(1)), (Int16(10), Int16(1)), (Int16(5), Int16(3)), (Int16(1), Int16(1))]
dmc = Dict(zip(ks, vs))

ks = [#
    Int16.([1, 2, 5, 6, 7]),
    Int16.([1, 2, 5, 6]),
    Int16.([1, 2, 5, 6, 7, 4]),
    Int16.([1]),
    Int16.([1, 2, 3]),
    Int16.([1, 2, 3, 4]),
    Int16.([1, 2]),
    Int16.([1, 2, 5]),
    Int16.([1, 2, 5, 6, 7, 8]),
    Int16.([1, 2, 5, 6, 7, 8, 4])
]
vs = [#
    (id = 7, children = Int16.([8])),
    (id = 6, children = Int16.([7])),
    (id = 8, children = Int16.([])),
    (id = 1, children = Int16.([2])),
    (id = 3, children = Int16.([])),
    (id = 4, children = Int16.([])),
    (id = 2, children = Int16.([3, 5])),
    (id = 5, children = Int16.([6])),
    (id = 9, children = Int16.([])),
    (id = 10, children = Int16.([]))
]
dchmc = Dict(zip(ks, vs))

T = SparseMatrixCSC([
    0 1.0 0 0 0 0 0 0 0 0
    0 0 0.3 0 0.7 0 0 0 0 0
    0 0 0 1.0 0 0 0 0 0 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1.0 0 0 0 0
    0 0 0 0 0.125 0 0.875 0 0 0
    0 0 0 0 0 0 0 0.875 0.125 0
    1.0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0.125 0 0 0.875
    1.0 0 0 0 0 0 0 0 0 0
])

R = [#
    (i = 1, j = 2, k = Int16(2)),
    (i = 2, j = 3, k = Int16(3)),
    (i = 2, j = 5, k = Int16(5)),
    (i = 3, j = 4, k = Int16(4)),
    (i = 4, j = 1, k = Int16(4)),
    (i = 5, j = 6, k = Int16(6)),
    (i = 6, j = 5, k = Int16(7)),
    (i = 6, j = 7, k = Int16(8)),
    (i = 7, j = 8, k = Int16(12)),
    (i = 8, j = 1, k = Int16(12)),
    (i = 7, j = 9, k = Int16(10)),
    (i = 9, j = 10, k = Int16(11)),
    (i = 10, j = 1, k = Int16(11)),
    (i = 9, j = 7, k = Int16(9))
]

res = CHMCAtomicSummary(#
    (1, 1, :C),
    [#
        (EFM = [6, 5, 6], Closures = [(6, 5)]),
        (EFM = [4, 1, 2, 5, 6, 7, 4], Closures = [(8, 1)]),
        (EFM = [4, 1, 2, 3, 4], Closures = [(4, 1)]),
        (EFM = [8, 7, 8], Closures = [(9, 7)]),
        (EFM = [4, 1, 2, 5, 6, 7, 8, 4], Closures = [(10, 1)])
    ],
    [0.09000000000000001, 0.56, 0.27, 0.01, 0.07],
    [1.0, 6.222222222222222, 3.0, 0.11111111111111112, 0.7777777777777778],
    dmc,
    dchmc,
    T,
    R
)

res_enum = CHMCAtomicSummary(#
    (1, 1, :C),
    [#
        (EFM = [6, 5, 6], Closures = [(6, 5)]),
        (EFM = [4, 1, 2, 5, 6, 7, 4], Closures = [(8, 1)]),
        (EFM = [4, 1, 2, 3, 4], Closures = [(4, 1)]),
        (EFM = [8, 7, 8], Closures = [(9, 7)]),
        (EFM = [4, 1, 2, 5, 6, 7, 8, 4], Closures = [(10, 1)])
    ],
    nothing,
    nothing,
    dmc,
    dchmc,
    T,
    R
)

efm_seq_1 = ["FDP", "F6P", "FDP"]
efm_seq_2 = ["Glc", "G6P", "F6P", "FDP", "DHAP"]
```

## Inputs

```julia
using MarkovWeightedEFMs
S = [#
  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 # Glc
  0 -1  0  0  0 -1  0  0  0  0  0  0  1  0  0  0 # ATP
  0  1 -1  0 -1  0  0  0  0  0  0  0  0  0  0  0 # G6P
  0  1  0  0  0  1  0  0  0  0  0  0  0 -1  0  0 # ADP
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 # 6PG
  0  0  0  0  1 -1  1  0  0  0  0  0  0  0  0  0 # F6P
  0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1  0 # Pi
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1 # H2O
  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0 # FDP
  0  0  0  0  0  0  0  1 -1  1 -1  0  0  0  0  0 # G3P
  0  0  0  0  0  0  0  1  1 -1  0 -1  0  0  0  0 # DHAP
]

v = [10, 10, 3, 3, 7, 8, 1, 7, 1, 1, 7, 7, 18, 18, 1, 1]

mets = [#
  "Glc",
  "ATP",
  "G6P",
  "ADP",
  "6PG",
  "F6P",
  "Pi",
  "H2O",
  "FDP",
  "G3P",
  "DHAP"
]

rxns = [#
  "Source Glc",
  "Hexokinase",
  "G6P dehydrogenase",
  "Sink 6PG",
  "Phosphoglucose isomerase",
  "6-phosphofructo-1-kinase",
  "Fructose 1,6-bisphosphatase",
  "Fructose-bisphosphate aldolase",
  "Triose phosphate isomerase",
  "Triose phosphate isomerase",
  "Sink G3P",
  "Sink DHAP",
  "Source ATP",
  "Sink ADP",
  "Sink Pi",
  "Source H2O"
]

smiles = [#
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
  "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)OP(=O)(O)O",
  "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
  "O=C1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O",
  "C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)OP(=O)(O)O",
  "[O-]P(=O)([O-])[O-]",
  "O",
  "C(C1C(C(C(O1)(COP(=O)(O)O)O)O)O)OP(=O)(O)O",
  "C([C@H](C=O)O)OP(=O)(O)O",
  "C(C(=O)COP(=O)(O)O)O"
]

atom = :C # carbon atom type for AEFMs
```

We can check that the flux vector satisfies the steady state requirements.

```@example required
all(S * v .== 0) # should evaluate as true
```



## Pre-processing data (wrapper)

The following function pre-processes the input metabolic network for
computing the AEFM weights for the specified atom type.

```julia
mdl, atom_info, logs = preprocess_all_for_atomic_chmc(S, v, mets, rxns, smiles, atom)
```

The variable `mdl` is a NamedTuple containing the updated stoichiometry
matrix, flux vector, metabolite/reaction names, metabolite SMILES strings,
reaction SMILES strings, and mapped reaction SMILES strings.

```@example required
keys(mdl)
```

The variable `logs` contains details about input metabolic network in
addition to listing the pseudometabolites and pseudoreactions dropped from
the network based on the input SMILES strings.

```@example required
keys(logs)
```

```@example required
print(logs.model_errors)
```

```@example required
logs.smiles_warnings
```

The variable `atom_info` contains the indices of all source metabolites in
the updated network and the number of occurrences for the input atom of
interest. It also contains an atom-mapping dictionary relating substrate-atom
positions to product-atom positions in each reaction. These are useful
for programmatically computing AEFMs across all source metabolite-atom
combinations of interest.

```@example required
keys(atom_info)
```

```@example required
atom_info.src_mets # source metabolite indices
```

```@example required
atom_info.max_src_met_atoms # counts of specified atom in each source metabolite
```

```@example required
atom_info.D # atom-mapping dictionary
```

## Computing ACHMC for a given metabolite/carbon atom state

The following atomic CHMC is rooted on the first carbon atom of the first
source metabolite in the stoichiometry.

```julia
I = (atom_info.src_mets[1], 1, atom) # initial state is 1st carbon of canonicalized glucose
res = steady_state_efm_distribution(mdl.S, mdl.v, mdl.ms, I, atom_info.D; verbose = false) # S must be Matrix{Int16}
```
```@example required
res
```

If we only wanted to enumerate the AEFMs, we would run:

```julia
res_enum = enumerate_atomic_efms(mdl.S, mdl.ms, I, atom_info.D, verbose = false)
```

Both functions produce the same output structure `res`, except that the
AEFM flux decomposition fields will be empty. The transition matrix will
also default to uniformly distributed probabilities along each row.

### Output

The output `res` is an immutable struct with 8 fields:

`res.i` is a tuple storing (i) the source metabolite index, (ii) source
metabolite atom index (based on canonicalized SMILES string), and (iii)
the atom type. This is a copy of the variable `I`.

```@example required
res.i
```

`res.e` is an array of AEFMs with all corresponding simple cycle closures.

```@example required
res.e
```

`res.p` is an array of AEFM probabilities normalized to one.

```@example required
res.p
```

`res.w` is an array of AEFM weights normalized by the (unimolecular)
reaction flux of the source metabolite.

```@example required
res.w
```

`res.dchmc` is a dictionary storing the ACHMC. The keys are the ACHMC
states (composed of Markov chain states in `res.dmc`). The values are the
ACHMC state and the Markov chain state children.

```@example required
res.dchmc
```

`res.dmc` is a dictionary converting Markov chain states to
metabolite-atom positions. The value `(0, 0)` always corresponds to the
external environment sink node (which connects back to the source
metabolite-atom state).

```@example required
res.dmc
```

`res.T` is a sparse array storing the ACHMC transition probability matrix.

```@example required
res.T
```

`res.R` is an array of tuples storing the reaction index/indices mapped to
each ACHMC transition matrix element.

```@example required
res.R
```

## Converting AEFM to sequence of metabolites

The corresponding AEFMs correspond to the movement of
metabolite/atom states through the reaction network. We can convert these
states into metabolites using `get_efm_metabolite_atom_indices`.  Note
that there is one fewer metabolite name than AEFM metabolite indices
because the pseudometabolite `(0, 0)` linking sink and source reactions is
omitted.

```julia
# First AEFM
efm_seq_1 = mets[first.(get_efm_metabolite_atom_indices(res, 1))]
```
```@example required
efm_seq_1
```

```julia
# Second AEFM
efm_seq_2 = mets[first.(get_efm_metabolite_atom_indices(res, 2))]
```
```@example required
efm_seq_2
```

## Visualizing the CHMC and mapped reactions

The following plotting function visualizes the ACHMC rooted on state `I`.
This is only recommended for exploring ACHMCs of small networks.

```julia
using GLMakie # Makie backend
GLMakie.activate!()

plot_atomic_chmc(res, S, mets, rs)
```

Each node in the main panel corresponds to a CHMC state
(metabolite and atomic index).

![ACHMC main panel](../assets/toy-network-1-chmc-makie-1.png)

Clicking on a CHMC transition will highlight
that transition and display the corresponding metabolic reaction on the upper
panel. The pair of purple highlighted atoms correspond to the movement of the
same atom from the LHS to RHS of the reaction.

![ACHMC main and upper panel](../assets/toy-network-1-chmc-makie-2.png)

Finally, the reaction and mapped reaction SMILES strings can also be plotted as
an SVG and previewed using a package like ElectronDisplay. If `fname != ""`,
the SVG is also saved to file. By default, `fname == ""` and the SVG is
not saved. The default canvas width and height are 1420 by 580 (pixels)
but these can be changed. If using ElectronDisplay and the image is cut
off, try resizing the plotting window or reducing the canvas dimensions.

```julia
using ElectronDisplay

# Reaction string
plot_mapped_reaction(rs[2], view=true, canvas_width = 1420, canvas_height = 580)
#plot_mapped_reaction(rs[2], "\path\to\save\name.svg", view=true)
```

![Reaction SMILES string](../assets/rs-2.svg)

```julia
# Mapped reaction string
plot_mapped_reaction(ms[2], view = true, canvas_width = 1420, canvas_height = 580)
#plot_mapped_reaction(ms[2], "\path\to\save\name.svg", view = true)
```

![Mapped reaction SMILES string](../assets/ms-2.svg)

