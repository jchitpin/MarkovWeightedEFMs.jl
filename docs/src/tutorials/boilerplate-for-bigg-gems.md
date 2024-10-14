# Boilerplate code for ACHMC analysis of BiGG models

The following code is provided to quickly construct ACHMC models of BiGG
metabolic models and others stored in the SBML file format. Each code
block is designed to be standalone with calculations saved to file and
re-loaded in a subsequent code block. Unfortunately, BiGG does not store
SMILES strings so these will need to be manually compiled by the user.

Note we recommend starting with relatively small networks (<500
metabolites and reactions in the original metabolic model) for
computational feasibility.

## Inputs (from BiGG)

This section shows how to extract relevant information from a BiGG
metabolic model called `e_coli_core.xml`

```julia
# Import file names
im_sbml = "e_coli_core.xml"

# Export file names
ex_stoich                   = "stoich.csv"
ex_metabolites_compartments = "metabolites-compartments.csv"
ex_reactions                = "reactions.csv"
ex_formulas                 = "metabolite-formulas.csv"

# Packages
using SBML, CSV, Tables

# Load SBML model
mdl = readSBML(im_sbml)
metabolites, reactions, S = stoichiometry_matrix(mdl)

# Dense stoichiometry matrix
S = Array(S)

# Metabolite names
mets = [mdl.species[m].name for m in metabolites]
mets = replace.(mets, " " => "_")

# Metabolite names concatenated with compartment
metsc = [#
    join([mdl.species[m].name, mdl.species[m].compartment], "_")
    for m in metabolites
]
metsc = replace.(metsc, " " => "_")

# Metabolite formulas
formulas = [mdl.species[m].formula for m in metabolites]
formulas[isnothing.(formulas)] .= ""

# Reaction names
rxns = [mdl.reactions[r].name for r in reactions]

# Export to text file
CSV.write(ex_stoich, Tables.table(S), header = false)
CSV.write(ex_metabolites, Tables.table(mets), header = false)
CSV.write(ex_metabolites_compartments, Tables.table(metsc), header = false)
CSV.write(ex_reactions, Tables.table(rxns), header = false)
CSV.write(ex_formulas, Tables.table(formulas), header = false)
```

## Pre-processing inputs (for AEFM enumeration)

This section pre-processes the BiGG inputs to meet the ACHMC requirements.

```julia
# Import file names
im_stoich = "stoich.csv"
im_mets   = "metabolites-compartments.csv"
im_rxns   = "reactions.csv"
im_smiles = "smiles-isomeric.csv"

# Export file names
ex_stoich = "stoichiometry-matrix-processed.csv"
ex_mets   = "metabolites-processed.csv"
ex_rxns   = "reactions-processed.csv"
ex_smiles = "smiles-isomeric-processed.csv"
ex_reaction_smiles        = "reaction-smiles-processed.csv"
ex_mapped_reaction_smiles = "mapped-reaction-smiles-strings-processed.csv"
ex_dict_carbon            = "dictionary-atom-tracing-carbon.csv"

# Packages
using CSV, Tables, MarkovWeightedEFMs, BenchmarkTools

## Load previous inputs
# Stoichiometry matrix
S = CSV.read(im_stoich, Tables.matrix, header = false)

# Metabolites
mets = vec(CSV.read(im_mets, Tables.matrix, header = false))

# Reactions
rxns = String.(vec(CSV.read(im_rxns, Tables.matrix, header = false)))

## Pre-processing
# (1) Identify problems with S/v inputs
errors = find_atomic_chmc_input_errors(S)
print(errors) # summary of errors associated with S

# (2) Clean S inputs
S2, mets2, rxns2 = correct_atomic_chmc_input_errors(errors, S, mets, rxns)
print(find_atomic_chmc_input_errors(S2)) # confirm errors have been fixed

# (3) Construct vector of smiles corresponding to the remaining metabolites in S
# The SMILES strings for pseudometabolites with no defined chemical structure
# are given an arbitrary SMILES of 'R' (or character that does not represent
# a periodic table element)
# SMILES strings matching S2
smiles3 = vec(CSV.read(im_smiles, Tables.matrix, header = false))

# (4) Remove pseudometabolites and reactions exceeding RXNMapper character limit
S4, mets4, rxns4, smiles4, i4 = correct_atomic_chmc_input_smiles(#
  S2, mets2, rxns2, smiles3
)
i4.dropped_rows_pseudometabolites # pseudometabolite rows removed from S2
i4.dropped_cols_pseudometabolites # pseudometabolite reactions removed from S2
i4.dropped_cols_rxnmapper_limit # reactions in S2 removed bc of RXNMapper limit
print(find_atomic_chmc_input_errors(S4)) # confirm no errors

# (5) Construct the reaction strings and map atoms via RXNMAPPER
smiles5 = canonicalize_smiles(smiles4) # smiles strings must be canonicalized!
rs5, ms5 = map_reaction_strings(S4, smiles5, rxns4, false)

# (6) Precompute atom tracing dictionary (for carbons)
amax = get_max_atoms(smiles5, :C)
D_C = precompute_atom_tracing_dictionary(S4, ms5, amax, :C)

# Export to text file
CSV.write(ex_stoich, Tables.table(S4), header = false)
CSV.write(ex_mets, Tables.table(mets4), header = false)
CSV.write(ex_rxns, Tables.table(rxns4), header = false, quotestrings = true)
CSV.write(ex_smiles, Tables.table(smiles5), header = false)
CSV.write(ex_reaction_smiles, Tables.table(rs5), header = false)
CSV.write(ex_mapped_reaction_smiles, Tables.table(ms5), header = false)
CSV.write(ex_dict_carbon, D_C, header = false)
```

## Pre-processing inputs (for AEFM weight assignment)

Assuming there is a steady state flux vector `v`, the pre-processing steps
are slightly different:

```julia
## Pre-processing
# (1) Identify problems with S/v inputs
errors = find_atomic_chmc_input_errors(S, v)
print(errors) # summary of errors associated with S/v

# (2) Clean S/v inputs
S2, v2, mets2, rxns2 = correct_atomic_chmc_input_errors(errors, S, v, mets, rxns)
print(find_atomic_chmc_input_errors(S2, v2)) # confirm errors have been fixed

# (3) Construct vector of smiles corresponding to the remaining metabolites in S
smiles3 = vec(CSV.read(im_smiles, Tables.matrix, header = false))

# (4) Remove pseudometabolites and reactions exceeding RXNMapper character limit
S4, v4, mets4, rxns4, smiles4, i4 = correct_atomic_chmc_input_smiles(#
  S2, v2, mets2, rxns2, smiles3
)
i4.dropped_rows_pseudometabolites # 33 pseudometabolite rows removed from S2
i4.dropped_cols_pseudometabolites # 46 pseudometabolite reactions removed from S2
i4.dropped_cols_rxnmapper_limit # 3 reactions in S2 removed bc of RXNMapper limit
print(find_atomic_chmc_input_errors(S4, v4)) # confirm no errors

# (5) Construct the reaction strings and map atoms via RXNMAPPER
smiles5 = canonicalize_smiles(smiles4) # smiles strings must be canonicalized!
rs5, ms5 = map_reaction_strings(S4, smiles5, rxns, false)

# (6) Precompute atom tracing dictionary (for carbons)
amax = get_max_atoms(smiles5, :C)
D_C = precompute_atom_tracing_dictionary(S4, ms5, amax, :C)
```

## Enumerating AEFMs across all source metabolite carbons

```julia
# Verbosity of AEFM enumeration function (provides a progress meter)
verbose = true

# File names
stoich_loc = "stoichiometry-matrix-processed.csv"
smiles_loc = "smiles-isomeric-processed.csv"
mets_loc   = "metabolites-processed.csv"
D_C_loc    = "dictionary-atom-tracing-carbon.csv"
D_N_loc    = "dictionary-atom-tracing-nitrogen.csv"
ms_loc     = "mapped-reaction-smiles-strings-processed.csv"
rxns_loc   = "reactions-processed.csv"

# Packages
using CSV, Tables, MarkovWeightedEFMs, JLD2, Dates

## Load final data
# Load stoichiometry matrix
S = Int16.(CSV.read(stoich_loc, Tables.matrix, header = false))

# Load SMILES strings matching the stoichiometry rows
smiles = vec(CSV.read(smiles_loc, Tables.matrix, header = false))

# Load metabolites
mets = vec(CSV.read(mets_loc, Tables.matrix, header = false))

# Load atom tracing dictionary
D_C = import_atom_tracing_dictionary(D_C_loc)

# Load mapped reaction smiles strings
ms = vec(CSV.read(#
    ms_loc, Tables.matrix, delim = ';', ignoreemptyrows = false, header = false
))
g(x) = ismissing(x) ? "" : x
ms = g.(ms)

# Load reactions
rxns = String.(vec(CSV.read(rxns_loc, Tables.matrix, header = false)))

## Enumerate atomic efms
# Identify indices of all source metabolites and number of carbons/nitrogens
srcs = get_source_metabolites(Int16.(S))
amax_C = get_max_atoms(smiles, :C)

start_time = Dates.now()
for k in eachindex(srcs)
    res_dir = "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    for l in 1:amax_C[srcs[k]]
        I = (srcs[k], l, :C)
        res = enumerate_atomic_efms(S, ms, I, D_C; verbose = verbose)
    end
end
end_time = Dates.now()
ttime = end_time - start_time

@info "It took $ttime to enumerate all carbon AEFMs."
```

## Computing AEFM weights across all source metabolite carbons

Decomposing fluxes onto AEFMs uses the function
`steady_state_efm_distribution()` requires the steady state fluxes `v`:

```julia
for k in eachindex(srcs)
    res_dir = "$k-" * mets[srcs[k]]
    isdir(res_dir) || mkdir(res_dir)
    for l in 1:amax_C[srcs[k]]
        I = (srcs[k], l, :C)
        res = steady_state_efm_distribution(S, v, ms, I, D_C; verbose = verbose)
    end
end
```



