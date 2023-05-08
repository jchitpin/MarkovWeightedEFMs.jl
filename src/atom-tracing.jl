# Helpers
function string_chemical_element()
  r = [#
    "A[cglmrstu]",
    "B[aehikr]?",
    "C[adeflmnorsu]?",
    "D[bsy]",
    "E[rsu]",
    "F[elmr]?",
    "G[ade]",
    "H[efgos]?",
    "I[nr]?",
    "Kr?",
    "L[airuv]",
    "M[dgnot]",
    "N[abdeiop]?",
    "Os?",
    "P[abdmortu]?",
    "R[abefghnu]",
    "S[bcegimnr]?",
    "T[abcehilm]",
    "U(u[opst])?",
    "V",
    "W",
    "Xe",
    "Yb?",
    "Z[nr]"
  ]
  return r
end
function regex_chemical_element()
  return Regex(join(string_chemical_element(), "|"))
end

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

# Call rxnmapper to trace the atoms in a reaction SMILES string
function trace_rxn_string(s::Vector{String})
  rxnmap = pyimport("rxnmapper")
  rxn_mapper = rxnmap.RXNMapper()
  skip1 = findall(occursin.(r"R", s))
  skip2  = findall(occursin.(r"^>>", s))
  skip3  = findall(occursin.(r">>$", s))
  skip = [skip1; skip2; skip3]

  idx = setdiff(1:length(s), skip)
  m = Vector{String}(undef, length(s))
  for i in 1:length(s)
    if i ∈ idx
      try
        tmp = rxn_mapper.get_attention_guided_atom_maps(#
          s[[i]], canonicalize_rxns=false
        )
        m[i] = tmp[1]["mapped_rxn"]
      catch e
        @error(#
          "Reaction $i string length exceeded the 512 character limit."
        )
      end
    else
      m[i] = ""
    end
  end

  @assert(#
    all([isassigned(m, i) for i in 1:length(m)]),
    join([#
      "The reactions above must be within the 512 character limit for ",
      "RXNMapper. Either: (i) remove redundant characters in the SMILES ",
      "string (if possible), or (ii) remove the offending reaction and ",
      "replace with corresponding sink and source reactions carrying the ",
      "same metabolic flux and stoichiometric coefficients."
    ])
  )

  return m
end

# Canonical SMILES and atom mappings according to rxnmapper
function canonicalize_and_atom_map(m::Vector{String})
  rxnmap = pyimport("rxnmapper")
  rxn_mapper = rxnmap.smiles_utils
  m = [split.(m[i], r"\.|(>>)") for i in 1:length(m)]
  return m .|> x -> rxn_mapper.canonicalize_and_atom_map.(x)
end
function canonicalize_and_atom_map(m::String)
  rxnmap = pyimport("rxnmapper")
  rxn_mapper = rxnmap.smiles_utils
  m = split(m, r"\.|(>>)")
  return rxn_mapper.canonicalize_and_atom_map.(m)
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

# Main function to get mapping of atom a between substrate and product in mapped reaction m
function trace_atoms(m::String, a::Tuple{Int64, String, Int64})
  # Error checking
  rChem = regex_chemical_element()
  @assert(a[2] != "H")
  @assert(!isnothing(match(rChem, a[2])))

  # Atom map of mapped SMILES from rxnmapper
  x = canonicalize_and_atom_map(uppercase(m))

  # Get index of first product
  subs, prods = split(m, ">>")
  subs = split(subs, ".")
  @assert(subs[a[1]] != "[H+]")
  idx_prods = length(subs) + 1
  @assert(a[1] < idx_prods)
  
  # Substrate index of a[1] w.r.t. substrates in mapped SMILES
  atab = countmemb(collect(filter(isletter, x[a[1]][1])))
  @assert(#
    haskey(atab, a[2]),
    "Chemical element a[2] not present in substrate index a[1]."
  )

  # Find mapped SMILES substrate index of chosen atom index
  sidx = collect(x[a[1]][1])
  filter!(isletter, sidx)
  filter!(x -> x != 'H', sidx)
  @assert(length(sidx) == length(x[a[1]][2]))
  ms_idx = findall(==(a[2][1]), sidx)[a[3]] # Index of atoms corresponding to the a[3] occurrence of atom a[2].
  idx_1 = x[a[1]][2][ms_idx] # Mapped SMILES substrate index corresponding to the a[3] occurrence of atom a[2].

  # Find the chosen product atom index from mapped SMILES substrate index
  for i in idx_prods:length(x)
    if string(idx_1) ∈ x[i][2]
      midx = collect(x[i][1])
      filter!(isletter, midx)
      filter!(x -> x != 'H', midx)
      @assert(length(midx) == length(x[i][2]))
      idx_2 = length(findall(==(a[2][1]), midx[1:findfirst(==(idx_1), x[i][2])]))
      return (i-idx_prods+1, a[2], idx_2)
    end
  end
end

function backup_trace_atoms(m::String, a::Tuple{String, String, Int64})
  # Error checking
  rChem = regex_chemical_element()
  @assert(a[1] != "H")
  @assert(!isnothing(match(rChem, a[1])))

  # Atom map of mapped SMILES from rxnmapper
  x = canonicalize_and_atom_map(uppercase(m))

  # Substrate index of a[1] w.r.t. substrates in mapped SMILES
  idx_sub = Int64(0)
  a_met = parse_formula(replace(a[1], r"(H\d+|H)" => ""))
  a_met = Dict(zip(collect(join(first.(a_met))), last.(a_met)))
  @assert(#
    haskey(a_met, a[2][1]),
    "Chemical element a[2] not present in molecular formula a[1]."
  )
  @assert(#
    a[3] <= a_met[a[2][1]],
    "Element index a[3] exceeds the number of atoms in molecular formula a[1]."
  )
  idx_max_subs = findfirst([x[i][2][1] == "1" for i in 2:length(x)])
  for i in 1:idx_max_subs
    y = collect(uppercase(x[i][1]))
    filter!(isletter, y)
    filter!(x -> x != 'H', y)
    if a_met == countmemb(String(y))
      idx_sub = i
      break
    end
  end

  # Find mapped SMILES substrate index of chosen atom index
  sidx = collect(x[idx_sub][1])
  filter!(isletter, sidx)
  filter!(x -> x != 'H', sidx)
  @assert(length(sidx) == length(x[idx_sub][2]))
  ms_idx = findall(==(a[2][1]), sidx)[a[3]] # Index of atoms corresponding to the a[3] occurrence of atom a[2].
  idx_1 = x[idx_sub][2][ms_idx] # Mapped SMILES substrate index corresponding to the a[3] occurrence of atom a[2].

  # Find the chosen product atom index from mapped SMILES substrate index
  for i in (idx_max_subs+1):length(x)
    if string(idx_1) ∈ x[i][2]
      midx = collect(x[i][1])
      filter!(isletter, midx)
      mf = countmemb(midx) # molecular formula of product
      mf = join(permutedims(hcat(collect(keys(mf)), collect(values(mf)))))
      filter!(x -> x != 'H', midx)
      @assert(length(midx) == length(x[i][2]))
      idx_2 = length(findall(==(a[2][1]), midx[1:parse(Int64, idx_1)]))
      return (product_idx=i-idx_max_subs, a=(mf, a[2], idx_2))
    end
  end
end

# Dictionary of occurrences of unique elements in a vector or characters in a string
function countmemb(itr::Vector{Char})
    d = Dict{String, Int64}()
    for val in itr
      d[string(val)] = get(d, string(val), 0) + 1
    end
    return d
end

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
function rxn_string(S::Matrix{<:Real}, smiles::Vector{String})
  # Error checking
  @assert(size(S,1) == length(smiles))

  # Construct SMILES reaction strings
  rs = Vector{String}(undef, size(S,2))

  # Construct SMILES reaction string
  for j in 1:length(rs)
    subs_idx = findall(<(0), @view S[:,j])
    prods_idx = findall(>(0), @view S[:,j])
    subs_coeff = abs.(@view S[findall(<(0), @view S[:,j]), j])
    prods_coeff = @view S[findall(>(0), @view S[:,j]), j]

    # Internal reactions
    if !isempty(subs_idx) && !isempty(prods_idx)
      @assert(#
        all(isinteger.(subs_coeff)),
        join([#
          "Stoichiometry column $j must contain integer-valued substrate
          coefficients."
        ])
      )
      @assert(#
        all(isinteger.(prods_coeff)),
        join([#
          "Stoichiometry column $j must contain integer-valued product
          coefficients."
        ])
      )
      subs = [#
        join(repeat([smiles[subs_idx][i]], Int64(subs_coeff[i])), ".")
        for i in 1:length(subs_idx)
      ]
      prods = [#
               join(repeat([smiles[prods_idx][i]], Int64(prods_coeff[i])), ".")
        for i in 1:length(prods_idx)
      ]
      rs[j] = join([#
          join(subs, "."),
         ">>",
         join(prods, ".")
      ])

      # LHS-vs-RHS atom balancing for internal reactions (except hydrogens)
      formula(x) = Dict(#
        parse_formula(molecularformula(smilestomol(x)))
      )
      f((k,v)) = k => Float64(v)

      if !isempty(subs) && !isempty(prods)
        lhs_formula = merge(+, formula.(subs)...)
        rhs_formula = merge(+, formula.(prods)...)
        filter!(x -> x.first != "H", lhs_formula)
        filter!(x -> x.first != "H", rhs_formula)
        diff = merge(-, rhs_formula, lhs_formula)
        @assert(#
          all(values(diff) .== 0),
          "Non-hydrogen atoms are not balanced in reaction $j"
        )
      end
    else
      # Source reaction
      if isempty(subs_idx) && !isempty(prods_idx)
        @assert(#
          length(prods_idx) == 1,
          "Source reactions must contain a single source metabolite."
        )
        rs[j] = join([smiles[prods_idx[1]], ">>"])
      end
      # Sink reaction
      if !isempty(subs_idx) && isempty(prods_idx)
        @assert(#
          length(subs_idx) == 1,
          "Sink reactions must contain a single sink metabolite."
        )
        rs[j] = join([smiles[subs_idx[1]], ">>"])
      end
    end
  end
  return rs
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

