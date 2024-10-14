# Atomic CHMC

## Public functions

```@docs
canonicalize_smiles
find_atomic_chmc_input_errors
correct_atomic_chmc_input_errors
print(res::CHMCAtomicErrorSummary)
correct_atomic_chmc_input_smiles
exchange_atomic_chmc_input_metabolites
get_source_metabolites
get_max_atoms
map_reaction_strings
precompute_atom_tracing_dictionary
MarkovWeightedEFMs.CHMC.Atomic.steady_state_efm_distribution(#
    S::Matrix{<:Integer},
    v::Vector{<:Real},
    ms::Vector{String},
    I::Tuple{Int64,Int64,Symbol};
    D::Dict{#
      NTuple{4,Int64},
      Tuple{Int64,Int64}
    }=Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}(),
    tmp_dir::String=""
)
get_efm_metabolite_atom_indices
get_efm_reaction_atom_indices
chmc_to_mc_matrix
export_atom_tracing_dictionary
import_atom_tracing_dictionary
```

## Index

```@index
Pages = ["chmc-atomic.md"]
```


