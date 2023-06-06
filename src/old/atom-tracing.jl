











## Unused
# Convert array of chemical equations into reaction SMILES strings
function rxn_string(S::Matrix{<:Real}, m::Vector{String}, df::DataFrame, df_id::Symbol)
  # Error checking
  @assert(issubset(m, df[!,df_id]))

  # Construct SMILES reaction strings, expanding reversible reactions
  smiles = Vector{String}(undef, size(S,2))

  # Index of matching stoichiometry row molecule with corresponding dataframe row
  idx = [findfirst(==(m[i]), df[!,df_id]) for i in 1:length(m)]

  # Construct SMILES reaction string
  for j in 1:length(smiles)
    subs = df[idx[findall(<(0), S[:,j])], :CanonicalSMILES]
    prods = df[idx[findall(>(0), S[:,j])], :CanonicalSMILES]
    smiles[j] = join([join(subs, "."), ">>", join(prods, ".")])
  end
  return smiles
end

# Check if atom counts are conserved between rxnmapper and molecular formulas
function flatten_atoms(D::Dict{String, Dict{String, Dict{String, Int64}}})
  res = Vector{Set{String}}(undef, length(D))
  for i in 1:length(D)
    tmp = Vector{String}()
    for j in 1:length(values(collect(values(D))[i]))
      push!(#
        tmp,
        join([#
          collect(keys(collect(values(D))[i]))[j],
          string(sum(values(collect(values(collect(values(D))[i]))[j])))
        ])
      )
    end
    res[i] = Set(tmp)
  end
  return res
end
function flatten_atoms(s::String)
  P = parse_formula(s)
  filter!(x -> x[1] != "H", P)
  return Set([join(P[i]) for i in 1:length(P)])
end
function isconserved(#
  D::Dict{String, Dict{String, Dict{String, Int64}}},
  s::Vector{String}
)
  bool = [occursin("(", s[i]) for i in 1:length(s)]
  @assert(#
    all(bool .== false),
    "Ensure s is a molecular formula and not SMILES string."
  )
  x = flatten_atoms(D)
  y = flatten_atoms.(s)
  if length(x) == length(y) && all([x[i] ∈ y for i in 1:length(x)])
    return true
  else
    return false
  end
end

# Python wrapper to translate KEGG/PubChem IDs, etc.
"""
    translate_id(x::String, y::String, xz::Vector{String})

Translates list of compound identities z in the form of x into y.

*x* is the type of identifier.

*y* is the desired identifier translation.

*z* is the list of identifiers in the form of *x*.

# Example
```julia
julia> # Translate the following PubChem CID to a KEGG ID
julia> translate_id("PubChem CID", "KEGG", ["22247451"])
100%|█████████████████████████████████████████████| 1/1 [00:00<00:00,  1.57it/s]
Dict{Any, Any} with 1 entry:
  "KEGG" => Dict{Any, Any}("22247451"=>"C00001")
```
"""
function translate_id(x::String, y::String, xz::Vector{String})
  cts = pyimport("CTSgetPy.CTSgetPy")
  @assert(x ∈ cts.CTS_options(), "The specified x id does not exist.")
  @assert(y ∈ cts.CTS_options(), "The specified y id does not exist.")
  return cts.CTSget(x, y, xz)
end

# Query Pubchem database via REST api
function get_cid_info(cids::Vector{Int64})
  reg = r"-|\+"
  df = CSV.File(#
    get_for_cids(#
      cids;
      properties="CanonicalSMILES,MolecularFormula",
      output="CSV"
    )
  ) |> DataFrame
  replace!(#
    x -> occursin(reg, x) ? replace(x, reg => "") : x,
    df.MolecularFormula
  )
  @assert(#
    all(.!occursin.('.', df[!,3])),
    "Ionized metabolites with a period in Canonical SMILES not allowed. "
  )
  return df
end

