# ACHMC (for AEFMs)

This section shows how to use the functions in MarkovWeightedEFMs.jl to
enumerate and assign AEFMs weights to the following multispecies
reaction network.

![Toy multispecies network](../assets/toy-network-2-achmc.png)

```@setup required
using MarkovWeightedEFMs
```

## Inputs

```@example required
S = [#
  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  1  0  0  0  0  0 -1  0  0  0  0  0  0  0 -1  0
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  1  0
  0  0  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  1  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1
]

v = [10, 10, 10, 13, 10, 10, 10, 7, 7, 7, 7, 7, 7, 7, 7, 3, 3, 3]

mets = [#
  "6-phospho-D-gluconate", # 6PG
  "NADP+",
  "CO2",
  "H+",
  "NADPH",
  "ribulose-5-phosphate", # R5P
  "ribose-5-phosphate", # 5RP
  "AIR",
  "5-phosphoribosyl-4-carboxy-5-aminoimidazole", # CAIR
  "aspartate",
  "ATP",
  "ADP",
  "O4P3-",
  "SAICAR",
  "H2O",
  "HCO3-"
]

rxns = [#
  "source 6PG",
  "source NADP+",
  "6PG dehydrogenase",
  "sink H+",
  "sink NADPH",
  "ribose-5-phosphate isomerase",
  "sink ribose-5-phosphate",
  "source AIR",
  "sink CO2",
  "source aspartate",
  "source ATP",
  "SAICAR synthetase",
  "sink ADP",
  "sink O4P3-",
  "sink SAICAR",
  "source H2O",
  "HCO3- formation",
  "sink HCO3-"
]

smiles = [#
  "O=C(O)C(O)C(O)C(O)C(O)COP(=O)(O)O"
  "NC(=O)c1ccc[n+]([C@]2O[C@](COP(=O)(O)OP(=O)(O)OC[C@]3O[C@](n4cnc5c(N)ncnc54)[C@@](OP(=O)(O)O)[C@@]3O)[C@@](O)[C@@]2O)c1"
  "O=C=O"
  "[H+]"
  "NC(=O)C1=CN([C@]2O[C@](COP(=O)(O)OP(=O)(O)OC[C@]3O[C@](n4cnc5c(N)ncnc54)[C@@](OP(=O)(O)O)[C@@]3O)[C@@](O)[C@@]2O)C=CC1"
  "O=C(CO)[C](O)[C](O)COP(=O)(O)O"
  "O=P(O)(O)OC[C@]1O[C@@](O)[C@@](O)[C@@]1O"
  "Nc1cncn1C1O[C@](COP(=O)(O)O)[C@@](O)[C@@]1O"
  "Nc1c(C(=O)O)ncn1[C@]1O[C@](COP(=O)(O)O)[C@@](O)[C@@]1O"
  "[N+][C](CC(=O)O)C(=O)[O-]"
  "Nc1ncnc2c1ncn2[C@]1O[C@](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@](O)[C@@]1O"
  "Nc1ncnc2c1ncn2[C@]1O[C@](COP(=O)(O)OP(=O)(O)O)[C@@](O)[C@@]1O"
  "O=P(O)(O)O"
  "Nc1c(C(=O)NC(CC(=O)O)C(=O)O)ncn1[C@]1O[C@](COP(=O)(O)O)[C@@](O)[C@@]1O"
  "O"
  "O=C(O)O"
]
```

We can check that the flux vector satisfies the steady state requirements.

```@example required
all(S * v .== 0) # should evaluate as true
```

## Pre-processing data

The following functions check for issues with the inputs. The first function
`find_atomic_chmc_input_errors` identifies possible problems with the
stoichiometry matrix and flux vector. These problems, except for the steady
state flux requirement, can be addressed via `correct_atomic_chmc_input_errors`.
Finally, the last function `correct_atomic_chmc_input_smiles` checks and fixes
problems relating to the SMILES strings.

```@example required
# Confirm there are no issues with stoichiometry matrix 
errors = find_atomic_chmc_input_errors(S, v)
print(errors) # summary of errors associated with S/v

# S and v have no errors so the inputs are returned
correct_atomic_chmc_input_errors(errors, S, mets, rxns)
# S, mets, rxns = correct_atomic_chmc_input_errors(errors, S, mets, rxns) # otherwise

# Correct issues associated with RXNMapper character limit and pseudometabolites
S, v, mets, rxns, smiles, logs = correct_atomic_chmc_input_smiles(S, v, mets, rxns, smiles)
```

At this point, the SMILES strings (matching the updated `mets` if there were
errors in the initial inputs) should be canonicalized.

```@example required
smiles = canonicalize_smiles(smiles)
```

## Atom mapping reactions

The reaction SMILES strings are next constructed from the metabolite SMILES and
the atom mapping is performed via RXNMapper. In this tutorial, we will be
constructing an atomic CHMC rooted on a particular source metabolite carbon.
We precompute an atom tracing dictionary mapping the (carbon) atom in the
stoichiometric copy of a substrate to its product across each reaction.

```@example required
# Construct atom traced SMILES strings
rs, ms = map_reaction_strings(S, smiles, rxns, false)

# Precompute atom tracing dictionary
atom = :C # carbon
atom_max = get_max_atoms(smiles, atom)
D_C = precompute_atom_tracing_dictionary(S, ms, atom_max, atom)

# Identify source metabolites
src_mets = get_source_metabolites(S)
max_src_met_carbons = atom_max[src_mets]
nothing # hide
```

## Computing ACHMC for a given metabolite/carbon atom state

The following atomic CHMC is rooted on the first carbon atom of the first
source metabolite in the stoichiometry matrix 6-phospho-D-gluconate.

```@example required
I = (src_mets[1], 1, atom) # initial state is 1st carbon of 6-phospho-D-gluconate
res = steady_state_efm_distribution(S, v, ms, I, D_C; verbose = false)
```

If we only wanted to enumerate the AEFMs, we would run:

```@example required
enumerate_atomic_efms(S, ms, I, D_C, verbose = false)
```

Both functions produce the same output structure `res`, except that the
AEFM flux decomposition field will be a vector of zeros.


## Converting AEFM to sequence of metabolites

The corresponding AEFMs correspond to the movement of
metabolite/atom states through the reaction network. We can convert these
states into metabolites using `get_efm_metabolite_atom_indices`.  Note
that there is one fewer metabolite name than AEFM metabolite indices
because the pseudometabolite `(0, 0)` linking sink and source reactions is
omitted.

```@example required
# First AEFM
mets[first.(get_efm_metabolite_atom_indices(res, 1))]
```

```@example required
# Second AEFM
mets[first.(get_efm_metabolite_atom_indices(res, 2))]
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

![ACHMC main panel](../assets/toy-network-2-chmc-makie-1.png)

Clicking on a CHMC transition will highlight
that transition and display the corresponding metabolic reaction on the upper
panel. The pair of purple highlighted atoms correspond to the movement of the
same atom from the LHS to RHS of the reaction.

![ACHMC main and upper panel](../assets/toy-network-2-chmc-makie-2.png)

Finally, the reaction and mapped reaction SMILES strings can also be plotted as
an SVG and previewed using a package like ElectronDisplay. If `fname != ""`,
the SVG is also saved to file. By default, `fname == ""` and the SVG is not
saved.

```julia
using ElectronDisplay

# Reaction string
plot_mapped_reaction(rs[3], view=true)
#plot_mapped_reaction(rs[3], "\path\to\save\name.svg", view=true)
```

![Reaction SMILES string](../assets/rs-3.svg)

```julia
# Mapped reaction string
plot_mapped_reaction(ms[3], view=true)
```

![Mapped reaction SMILES string](../assets/ms-3.svg)

