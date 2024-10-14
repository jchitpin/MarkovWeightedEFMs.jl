function dict_2_to_1_letter_chemical_element()
    d = Dict(#
        "Ac"  => 'å',
        "Ag"  => '⏦',
        "Al"  => '♾',
        "Am"  => '́',
        "Ar"  => '⋰',
        "As"  => 'æ',
        "At"  => 'ℵ',
        "Au"  => '≌',
        #"B"   => '϶',
        "Ba"  => '‷',
        "Be"  => '‶',
        "Bh"  => '‵',
        "Bi"  => '∽',
        "Bk"  => '⋍',
        "Br"  => '⋿',
        #"C"   => '⋒',
        "Ca"  => 'Χ',
        "Cd"  => '∷',
        "Ce"  => '⩴',
        "Cf"  => '⋓',
        "Cl"  => '⫤',
        "Cm"  => '⤋',
        "Co"  => '⟱',
        "Cr"  => 'Δ',
        "Cs"  => 'Ð',
        "Cu"  => 'Ϝ',
        "Db"  => '†',
        "Ds"  => 'ℸ',
        "Dy"  => '☡',
        "Er"  => '⪘',
        "Es"  => '♪',
        "Eu"  => '⏧',
        #"F"   => '≒',
        "Fe"  => '⤯',
        "Fm"  => '⤬',
        "Fr"  => '♀',
        "Ga"  => 'γ',
        "Gd"  => '≥',
        "Ge"  => '♊',
        #"H"   => '̂',
        "He"  => '⩯',
        "Hf"  => 'ħ',
        "Hg"  => '♡',
        "Ho"  => '⚥',
        "Hs"  => '⊹',
        #"I"   => '⟺',
        "In"  => '⨌',
        "Ir"  => '∭',
        #"K"   => 'Κ',
        "Kr"  => 'Ϟ',
        "La"  => 'Ł',
        "Li"  => 'ł',
        "Lr"  => '♂',
        "Lu"  => '✠',
        "Md"  => '↧',
        "Mg"  => '↤',
        "Mn"  => '↦',
        "Mo"  => '↥',
        "Mt"  => '♂',
        #"N"   => '∇',
        "Na"  => '⊼',
        "Nb"  => '≉',
        "Nd"  => '≭',
        "Ne"  => '♮',
        "Ni"  => '⋪',
        "No"  => '⇼',
        "Np"  => '≇',
        #"O"   => 'ø',
        "Os"  => 'Ø',
        #"P"   => '̡',
        "Pa"  => '∥',
        "Pb"  => '▱',
        "Pd"  => '▰',
        "Pm"  => '∂',
        "Po"  => '⪣',
        "Pr"  => 'ɤ',
        "Pt"  => '⬠',
        "Pu"  => '⬟',
        "Ra"  => '˔',
        "Rb"  => '⟩',
        "Re"  => '⥇',
        "Rf"  => 'ʼ',
        "Rg"  => '⌉',
        "Rh"  => '⤫',
        "Rn"  => '⤰',
        "Ru"  => '”',
        #"S"   => '♐',
        "Sb"  => '𝖠',
        "Sc"  => '𝖺',
        "Se"  => '𝖡',
        "Sg"  => '𝖻',
        "Si"  => '𝖢',
        "Sm"  => '𝖼',
        "Sn"  => '𝖣',
        "Sr"  => '𝖽',
        "Ta"  => 'τ',
        "Tb"  => '♉',
        "Tc"  => '⫶',
        "Te"  => 'ʧ',
        "Th"  => 'þ',
        "Ti"  => '∴',
        "Tl"  => 'θ',
        "Tm"  => '⟀',
        #"U"   => '˘',
        "Uub" => '⇑',
        "Uuh" => '⤒',
        "Uuo" => '⇕',
        "Uup" => '⥮',
        "Uuq" => 'Υ',
        "Uus" => '⤊',
        "Uut" => '⟰',
        #"V"   => '⌅',
        #"W"   => '∧',
        "Xe"  => 'ξ',
        #"Y"   => '¥',
        "Yb"  => 'ʒ',
        "Zn"  => 'ζ',
        "Zr"  => 'Ƶ',
    )
    return d
end
function dict_1_to_2_letter_chemical_element()
    d = dict_2_to_1_letter_chemical_element()
    return Dict(zip(collect(values(d)), collect(keys(d))))
end

function find_atomic_chmc_input_errors_fluxes(S::Matrix{<:Real}, v::Vector{<:Real})
    # Error checking
    sanitize_stoichiometry_fluxes(S, v)

    # Check for steady state flux
    abs_flux_error = sum(abs.(S * v))

    # Check there is some flux through all reactions
    zero_flux_rxns = findall(==(0), v)

    # Check that fluxes are non-negative
    neg_rxns = findall(<(0), v)

    return (abs_flux_error, zero_flux_rxns, neg_rxns)
end

function find_atomic_chmc_input_errors_rest(S::Matrix{<:Real})
    # Check for duplicate reactions
    Scols = mapslices(x -> [x], S, dims = 1)[:]
    uniq_rxns = Vector{Vector{<:Real}}()
    dupl_rxns = Vector{Int64}()
    j = 0
    while !isempty(Scols)
        j += 1
        if Scols[1] ∉ uniq_rxns
            push!(uniq_rxns, popfirst!(Scols))
        else
            popfirst!(Scols)
            push!(dupl_rxns, j)
        end
    end

    # Check for reaction with no species involved
    empty_reactions = Vector{Int64}()
    for j in 1:size(S, 2)
        if all(iszero(@view S[:, j]))
            push!(empty_reactions, j)
        end
    end

    # Check for metabolite not involved in any reaction
    unused_metabolites = Vector{Int64}()
    for i in 1:size(S, 1)
        if all(iszero(@view S[i, :]))
            push!(unused_metabolites, i)
        end
    end

    # Identify number of source/sink reactions
    src_a_rxns = Vector{Int64}() # (a) single product; abs(stoich) == 1
    src_b_rxns = Vector{Int64}() # (b) single product; abs(stoich) != 1
    src_c_rxns = Vector{Int64}() # (c) multi products; abs.(stoich) .== 1
    src_d_rxns = Vector{Int64}() # (d) multi products; abs.(stoich) .!= 1
    snc_a_rxns = Vector{Int64}() # (a) single substrate; abs(stoich) == 1
    snc_b_rxns = Vector{Int64}() # (b) single substrate; abs(stoich) != 1
    snc_c_rxns = Vector{Int64}() # (c) multi substrates; abs.(stoich) .== 1
    snc_d_rxns = Vector{Int64}() # (d) multi substrates; abs.(stoich) .!= 1
    for j in 1:size(S, 2)
        pidx = findall(>(0), @view S[:, j])
        nidx = findall(<(0), @view S[:, j])
        pstoich = S[pidx, j]
        nstoich = abs.(S[nidx, j])
        if isempty(nidx) && !isempty(pidx)
            if length(pidx) == 1 && pstoich[1] == 1
                push!(src_a_rxns, j)
            elseif length(pidx) == 1 && pstoich[1] != 1
                push!(src_b_rxns, j)
            elseif length(pidx) > 1 && all(pstoich .== 1)
                push!(src_c_rxns, j)
            elseif length(pidx) > 1 && all(pstoich .!= 1)
                push!(src_d_rxns, j)
            end
        elseif !isempty(nidx) && isempty(pidx)
            if length(nidx) == 1 && nstoich[1] == 1
                push!(snc_a_rxns, j)
            elseif length(nidx) == 1 && nstoich[1] != 1
                push!(snc_b_rxns, j)
            elseif length(nidx) > 1 && all(nstoich .== 1)
                push!(snc_c_rxns, j)
            elseif length(nidx) > 1 && all(nstoich .!= 1)
                push!(snc_d_rxns, j)
            end
        end
    end

    # Check for internal reactions with non-integer stoichiometries
    non_int_stoich_rxns = Vector{Int64}()
    for j in 1:size(S, 2)
        nidx = findall(<(0), @view S[:, j])
        pidx = findall(>(0), @view S[:, j])
        if !isempty(nidx) && !isempty(pidx)
            if !all(isinteger.(@view S[:, j]))
                push!(non_int_stoich_rxns, j)
            end
        end
    end

    res = (#
        dupl_rxns,
        non_int_stoich_rxns,
        src_a_rxns,
        src_b_rxns,
        src_c_rxns,
        src_d_rxns,
        snc_a_rxns,
        snc_b_rxns,
        snc_c_rxns,
        snc_d_rxns,
        empty_reactions,
        unused_metabolites
    )
    return res
end

function correct_atomic_chmc_input_errors_inner(#
    res::CHMCAtomicErrorSummary,
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String}
)
    # Flux vector must be a float
    v = Float64.(v)

    ## Fixing stoichiometry matrix and flux vector (order matters!)
    # (1) Remove metabolites not participating in any reaction
    S = S[setdiff(1:size(S, 1), res.unused_metabolites), :]
    mets = mets[setdiff(1:length(mets), res.unused_metabolites)]

    # (2) Flip stoichiometric coefficients such that fluxes are strictly positive
    for j in res.reactions_with_negative_flux
        to_pos = findall(<(0), S[:, j])
        to_neg = findall(>(0), S[:, j])
        S[to_pos,j] = S[to_pos, j] .* -1
        S[to_neg,j] = S[to_neg, j] .* -1
        v[j] = v[j] * -1
    end

    # (3) Ensure all source reactions are unimolecular with stoichiometry of 1
    for j in res.reactions_unimolecular_source_stoichiometry_not_one
        i = findfirst(>(0), @view S[:, j])
        v[j] = abs(S[i, j] * v[j])
        S[i, j] = 1
    end

    # (4) Ensure all sink reactions are unimolecular with stoichiometry of 1
    for j in res.reactions_unimolecular_sink_stoichiometry_not_one
        i = findfirst(<(0), @view S[:, j])
        v[j] = abs(S[i, j] * v[j])
        S[i, j] = -1
    end

    # (5) Create new reactions to split multimolecular sources
    S_sources = Vector{Matrix{Float64}}()
    v_sources = Vector{Vector{Float64}}()
    r_sources = Vector{Vector{String}}()
    jj = [#
        res.reactions_multimolecular_source_stoichiometry_one;
        res.reactions_multimolecular_source_stoichiometry_not_one
    ]
    for j in jj
        idx = findall(>(0), @view S[:, j])
        stoich = S[idx, j]
        S_tmp = zeros(Float64, size(S, 1), length(idx))
        v_tmp = zeros(Float64, length(idx))
        r_tmp = Vector{String}(undef, length(idx))
        for i in eachindex(idx)
            S_tmp[idx[i], i] = 1
            v_tmp[i] = v[j] * stoich[i]
            r_tmp[i] = join([rxns[j], "_", i])
        end
        push!(S_sources, S_tmp)
        push!(v_sources, v_tmp)
        push!(r_sources, r_tmp)
    end
    if !isempty(S_sources)
        S_sources = reduce(hcat, S_sources)
        v_sources = reduce(vcat, v_sources)
        r_sources = reduce(vcat, r_sources)
    else
        S_sources = Matrix{Float64}(undef, size(S, 1), 0)
    end

    # (6) Create new reactions to split multimolecular sinks
    S_sinks = Vector{Matrix{Float64}}()
    v_sinks = Vector{Vector{Float64}}()
    r_sinks = Vector{Vector{String}}()
    jj = [#
        res.reactions_multimolecular_sink_stoichiometry_one;
        res.reactions_multimolecular_sink_stoichiometry_not_one
    ]
    for j in jj
        idx = findall(<(0), @view S[:, j])
        stoich = S[idx, j]
        S_tmp = zeros(Float64, size(S, 1), length(idx))
        v_tmp = zeros(Float64, length(idx))
        r_tmp = Vector{String}(undef, length(idx))
        for i in eachindex(idx)
            S_tmp[idx[i], i] = -1
            v_tmp[i] = abs(v[j] * stoich[i])
            r_tmp[i] = join([rxns[j], "_", i])
        end
        push!(S_sinks, S_tmp)
        push!(v_sinks, v_tmp)
        push!(r_sinks, r_tmp)
    end
    if !isempty(S_sinks)
        S_sinks = reduce(hcat, S_sinks)
        v_sinks = reduce(vcat, v_sinks)
        r_sinks = reduce(vcat, r_sinks)
    else
        S_sinks = Matrix{Float64}(undef, size(S, 1), 0)
    end

    # (5) Create sources/sinks for internal reactions with non-integer stoichiometries
    S_int, v_int, rxns_int = convert_internal_rxn_to_sources_sinks(#
        S,
        v,
        rxns,
        res.reactions_with_non_integer_stoichiometries
    )

    # (6) Remove empty reactions, those with zero fluxes, multimolecular
    # sources/sinks, and internal reactions with non-integer stoichiometries
    rm_idx = [#
        res.reactions_empty;
        res.reactions_with_zero_flux;
        res.reactions_multimolecular_source_stoichiometry_one;
        res.reactions_multimolecular_source_stoichiometry_not_one;
        res.reactions_multimolecular_sink_stoichiometry_one;
        res.reactions_multimolecular_sink_stoichiometry_not_one;
        res.reactions_with_non_integer_stoichiometries
    ]
    S = S[:, setdiff(1:size(S, 2), rm_idx)]
    v = v[setdiff(1:length(v), rm_idx)]
    rxns = rxns[setdiff(1:length(rxns), rm_idx)]

    # (7) Combine remaining stoichiometry matrix/inputs with newly created ones
    if !isempty(hcat(S_sources, S_sinks, S_int))
        S = hcat(S, S_sources, S_sinks, S_int)
        v = [v; v_sources; v_sinks; v_int]
        rxns = [rxns; r_sources; r_sinks; rxns_int]
    end

    # (8) Aggregate duplicate reactions
    Scols = mapslices(x -> [x], S, dims = 1)[:]
    v2 = deepcopy(v)
    rxns2 = deepcopy(rxns)
    SS = Vector{Vector{Float64}}()
    vv = Vector{Float64}()
    rr = Vector{String}()
    while !isempty(Scols)
        if Scols[1] ∉ SS
            push!(SS, popfirst!(Scols))
            push!(rr, popfirst!(rxns2))
            push!(vv, popfirst!(v2))
        else
            idx = findfirst(==(Scols[1]), SS)
            popfirst!(Scols)
            rr[idx] = join([rr[idx], "_", popfirst!(rxns2)])
            vv[idx] = vv[idx] + popfirst!(v2)
        end
    end
    S = reduce(hcat, SS)
    v = vv
    rxns = rr

    # (9) Add missing sources/sinks
    src_rxns = Vector{Vector{Float64}}()
    snc_rxns = Vector{Vector{Float64}}()
    src_fluxes = Vector{Float64}()
    snc_fluxes = Vector{Float64}()
    src_rxn_names = Vector{String}()
    snc_rxn_names = Vector{String}()
    k = 0
    for i in 1:size(S, 1)
        nidx = findall(<(0), @view S[i, :])
        pidx = findall(>(0), @view S[i, :])
        if isempty(nidx) && !isempty(pidx)
            k += 1
            tmp = zeros(Float64, size(S, 1))
            tmp[i] = -1
            push!(snc_rxns, tmp)
            push!(snc_fluxes, abs(sum(S[i, pidx] .* v[pidx])))
            push!(snc_rxn_names, "Sink_rxn_$k")
        elseif !isempty(nidx) && isempty(pidx)
            k += 1
            tmp = zeros(Float64, size(S, 1))
            tmp[i] = +1
            push!(src_rxns, tmp)
            push!(src_fluxes, abs(sum(S[i, nidx] .* v[nidx])))
            push!(src_rxn_names, "Source_rxn_$k")
        end
    end
    if !all(isempty.([src_rxns, snc_rxns]))
        S = hcat(S, reduce(hcat, [src_rxns; snc_rxns]))
        v = [v; src_fluxes; snc_fluxes]
        rxns = [rxns; src_rxn_names; snc_rxn_names]
    end
    #return Int16.(S), v, mets, rxns
    return S, v, mets, rxns
end

function exchange_atomic_chmc_input_metabolites_inner(#
    S::Matrix{Int16},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    ii::Vector{Int64}
)
    # Looping over each metabolite to remove
    for i in ii
        # Identify all internal reactions that produce/consume metabolite i
        rxns_prod = findall(>(0), S[i, :])
        rxns_cons = findall(<(0), S[i, :])

        # Add 2 rows to stoichiometry matrix corresponding to i sink/source
        S = [S; zeros(2, size(S, 2))]
        push!(mets, mets[i] * "_sink")
        push!(mets, mets[i] * "_source")
        push!(smiles, smiles[i])
        push!(smiles, smiles[i])

        # For all internal reactions producing/consuming i, replace with sink/source
        S[end-1, rxns_prod] = S[i, rxns_prod]
        S[end, rxns_cons] = S[i, rxns_cons]

        # Summed sink/source fluxes
        vi_sink = sum(v[rxns_prod])
        vi_source = sum(v[rxns_cons])

        # Unimolecular sink reaction for metabolite i sink
        bools = [isempty(findall(<(0), S[:, j])) for j in rxns_prod]
        if !all(bools) # sink flux does not exist for metabolite i
            sink = zeros(size(S, 1))
            sink[end-1] = -1
            S = reduce(hcat, [S, sink])
            push!(v, vi_sink)
            push!(rxns, mets[i] * "_sink")
        else # sink flux exists so inject missing flux there
            v[rxns_prod[findfirst(bools)]] += vi_sink
        end

        # Unimolecular source reaction for metabolite i source
        bools = [isempty(findall(>(0), S[:, j])) for j in rxns_prod]
        if !all(bools) # source flux does not exist for metabolite i
            source = zeros(size(S, 1))
            source[end] = 1
            S = reduce(hcat, [S, source])
            push!(v, vi_source)
            push!(rxns, mets[i] * "_source")
        else # source flux exists so inject missing flux there
            v[rxns_cons[findfirst(bools)]] += vi_source
        end
    end

    # Remove rows ii from stoichiometry matrix
    S = S[setdiff(1:size(S, 1), ii), :]
    deleteat!(mets, ii)
    deleteat!(smiles, ii)

    return Int16.(S), v, mets, rxns, smiles
end

function reactions_with_pseudometabolites(S::Matrix{<:Real}, smiles::Vector{String})
    formula(x) = Dict(parse_formula(molecularformula(smilestomol(x))))
    idx = Vector{Int64}()
    for j in 1:size(S, 2)
        try
            formula.(smiles[findall(!=(0), S[:, j])])
        catch e
            if isa(e, ErrorException)
                push!(idx, j)
            end
        end
    end

    return idx
end

function test_rxn_string(s::Vector{String})
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.RXNMapper()
    skip1  = findall(occursin.(r"^>>", s))
    skip2  = findall(occursin.(r">>$", s))
    skip = [skip1; skip2]
    idx = Vector{Int64}()
    for i in setdiff(1:length(s), skip)
        try
            tmp = rxn_mapper.get_attention_guided_atom_maps(#
                s[[i]], canonicalize_rxns = false
            )
        catch e
            push!(idx, i)
        end
    end

    return idx
end

function isolated_pseudometabolites_inner(S::Matrix{<:Real}, smiles::Vector{String})
    formula(x) = Dict(parse_formula(molecularformula(smilestomol(x))))
    idx = Vector{Int64}()
    for i in 1:size(S, 1)
        try
            formula(smiles[i])
        catch e
            push!(idx, i)
        end
    end

    return idx
end

function isolated_pseudometabolites(S::Matrix{<:Real}, smiles::Vector{String})
    rxns_rm = Vector{Int64}()
    for i in isolated_pseudometabolites_inner(S, smiles)
        append!(rxns_rm, findall(!=(0), S[i, :]))
    end
    return rxns_rm
end

function update_model_by_reactions(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    rxns::Vector{String},
    jj::Vector{Int64}
)
    S_int, v_int, rxns_int = convert_internal_rxn_to_sources_sinks(S, v, rxns, jj)
    S = S[:, setdiff(1:size(S, 2), jj)]
    v = v[setdiff(1:length(v), jj)]
    rxns = rxns[setdiff(1:length(rxns), jj)]
    S = hcat(S, S_int)
    v = [v; v_int]
    rxns = [rxns; rxns_int]
    return S, v, rxns
end

function correct_atomic_chmc_input_smiles_inner(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    mets::Vector{String},
    rxns::Vector{String},
    smiles::Vector{String},
    H::Bool
)
    # Error checking
    sanitize_correct_atomic_chmc_input_smiles(S, v, mets, rxns, smiles)

    # (1) Identify internal reactions to remove because pseudometabolite SMILES
    # are undefined, add sources/sinks to balance fluxes, and remove reactions
    rm_idx_psd = reactions_with_pseudometabolites(S, smiles)
    S, v, rxns = update_model_by_reactions(S, v, rxns, rm_idx_psd)

    # (2) Update reaction SMILES strings for (3) below
    rs = rxn_string(S, smiles, rxns, H)

    # (3) Identify internal reactions to remove because they exceed the RXNMapper
    # character limit, add sources/sinks to balance fluxes, and remove reactions
    rm_idx_limit = test_rxn_string(rs)
    S, v, rxns = update_model_by_reactions(S, v, rxns, rm_idx_limit)

    # (4) Pseudometabolites are now disconnected from network with sink/sources.
    # Identify their corresponding sink/source reactions
    rxns_rm = isolated_pseudometabolites(S, smiles)

    # (5) Remove sinks/sources (6) from the stoichiometry matrix, reactions, fluxes
    S = S[:, setdiff(1:size(S, 2), rxns_rm)]
    v = v[setdiff(1:length(v), rxns_rm)]
    rxns = rxns[setdiff(1:length(rxns), rxns_rm)]

    # (6) Remove orphaned pseudometabolites participating in no reactions
    rm_idx_psd_mets = findall([iszero(S[i, :]) for i in 1:size(S, 1)])
    S = S[setdiff(1:size(S, 1), rm_idx_psd_mets), :]
    mets = mets[setdiff(1:length(mets), rm_idx_psd_mets)]
    smiles = smiles[setdiff(1:length(smiles), rm_idx_psd_mets)]

    # Aggregate duplicate source/sink reactions
    Scols = mapslices(x -> [x], S, dims = 1)[:]
    v2 = deepcopy(v)
    rxns2 = deepcopy(rxns)
    SS = Vector{Vector{Float64}}()
    vv = Vector{Float64}()
    rr = Vector{String}()
    while !isempty(Scols)
        if Scols[1] ∉ SS
            push!(SS, popfirst!(Scols))
            push!(rr, popfirst!(rxns2))
            push!(vv, popfirst!(v2))
        else
            idx = findfirst(==(Scols[1]), SS)
            popfirst!(Scols)
            rr[idx] = join([rr[idx], "_", popfirst!(rxns2)])
            vv[idx] = vv[idx] + popfirst!(v2)
        end
    end
    S = reduce(hcat, SS)
    v = vv
    rxns = rr

    # Aggregate removed rows/columns
    rm = (#
        dropped_rows_pseudometabolites = rm_idx_psd_mets,
        dropped_cols_pseudometabolites = rm_idx_psd,
        dropped_cols_rxnmapper_limit = rm_idx_limit
    )

  return Int16.(S), v, mets, rxns, smiles, rm
end

# Dictionary of occurrences of unique characters in an array of characters
function countmemb(itr::Vector{String})
    d = Dict{String, Int64}()
    for val in itr
        d[val] = get(d, val, 0) + 1
    end
    return d
end

function dict_regex()
    d = dict_2_to_1_letter_chemical_element()
    dd = Vector{Pair{Regex,Char}}()
    for k in keys(d)
        push!(dd, Regex(k) => d[k])
    end
    return dd
end

function parse_from_2(s::String)
    ss = filter(isletter, s)
    d = dict_2_to_1_letter_chemical_element()
    dd = dict_regex()
    return collect(uppercase(replace(ss, dd...)))
end

function parse_from_2_old(s::String)
    chars = Vector{Char}()
    ss = filter(isletter, s)
    d = dict_2_to_1_letter_chemical_element()
    i = 0
    while i < length(ss)
        i += 1
        if i < length(ss) && haskey(d, ss[i:i+1])
            push!(chars, only(d[ss[i:i+1]]))
            i += 1
        else
            push!(chars, only(uppercase(ss[i])))
        end

    end
    return chars
end

function parse_from_1(c::Vector{Char})
    strs = Vector{String}()
    d = dict_1_to_2_letter_chemical_element()
    for i in eachindex(c)
        if haskey(d, c[i])
            push!(strs, d[c[i]])
        else
            push!(strs, string(c[i]))
        end
    end

    return strs
end

# Convert array of chemical equations into reaction SMILES strings
function rxn_string(#
    S::Matrix{<:Real},
    smiles::Vector{String},
    rxns::Vector{String},
    H::Bool = false
)
    # Error checking
    @assert(size(S, 1) == length(smiles))
    for i in eachindex(smiles)
        @assert(#
            isnothing(match(r"\.", smiles[i])),
            "The input smiles STRING at index $i cannot contain a period."
        )
    end

    # Construct SMILES reaction strings
    rs = Vector{String}(undef, size(S, 2))

    # Construct SMILES reaction string
    for j in eachindex(rs)
        subs_idx = findall(<(0), @view S[:, j])
        prods_idx = findall(>(0), @view S[:, j])
        subs_coeff = abs.(@view S[findall(<(0), @view S[:, j]), j])
        prods_coeff = @view S[findall(>(0), @view S[:, j]), j]

        # Internal reactions
        if !isempty(subs_idx) && !isempty(prods_idx)
            @assert(#
                all(isinteger.(subs_coeff)),
                join([#
                    "Stoichiometry column of reaction $(rxns[j]) must contain ",
                    "integer-valued substrate coefficients."
                ])
            )
            @assert(#
                all(isinteger.(prods_coeff)),
                join([#
                    "Stoichiometry column of reaction $(rxns[j]) must contain ",
                    "integer-valued product coefficients."
                ])
            )
            subs = [#
                join(repeat([smiles[subs_idx][i]], Int64(subs_coeff[i])), ".")
                for i in eachindex(subs_idx)
            ]
            prods = [#
                join(repeat([smiles[prods_idx][i]], Int64(prods_coeff[i])), ".")
                for i in eachindex(prods_idx)
            ]
            rs[j] = join([#
                join(subs, "."),
                ">>",
                join(prods, ".")
            ])

            # LHS-vs-RHS atom balancing for internal reactions (except hydrogen)
            formula(x) = countmemb(parse_from_1(parse_from_2(x)))

            if !isempty(subs) && !isempty(prods)
                lhs_formula = merge(+, formula.(subs)...)
                rhs_formula = merge(+, formula.(prods)...)
                if H == false
                    filter!(x -> x.first != "H", lhs_formula)
                    filter!(x -> x.first != "H", rhs_formula)
                end
                diff = merge(-, rhs_formula, lhs_formula)
                @assert(#
                    all(iszero(collect(values(diff)))),
                    "Atoms are not balanced in reaction $(rxns[j])."
                )
            end
        else
            # Source reaction
            if isempty(subs_idx) && !isempty(prods_idx)
                @assert(#
                    length(prods_idx) == 1,
                    "Source reactions must contain a single source metabolite."
                )
                rs[j] = join([">>", smiles[prods_idx[1]]])
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

# Call rxnmapper to trace the atoms in a reaction SMILES string
function trace_rxn_string(s::Vector{String})
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.RXNMapper()
    skip1 = findall(occursin.(r"R", s))
    skip2  = findall(occursin.(r"^>>", s))
    skip3  = findall(occursin.(r">>$", s))
    skip = [skip1; skip2; skip3]

    # Find all pairs of reversible reactions
    splits = [split.(s[i], ">>") for i in eachindex(s)]
    rev_pairs = Vector{Tuple{Int16,Int16}}()
    for i in eachindex(splits)
        if !any(isempty.(splits[i]))
            if [splits[i][2], splits[i][1]] ∈ splits
                idx = findfirst(==([splits[i][2], splits[i][1]]), splits)
                if (idx, i) ∉ rev_pairs && idx != i
                    push!(rev_pairs, (i, idx))
                end
            end
        end
    end

    idx = setdiff(1:length(s), [skip; last.(rev_pairs)])
    m = Vector{String}(undef, length(s))
    for i in eachindex(s)
        if i ∈ idx
            try
                tmp = rxn_mapper.get_attention_guided_atom_maps(#
                    s[[i]], canonicalize_rxns=false
                )
                m[i] = tmp[1]["mapped_rxn"]
            catch e
                @error("Reaction $i string length exceeds 512 character limit.")
            end
        else
            m[i] = ""
        end
    end

    # Add reversible reaction mappings by flipping LHS/RHS.
    for i in eachindex(rev_pairs)
        lhs, rhs = split(m[first(rev_pairs[i])], ">>")
        m[last(rev_pairs[i])] = join([rhs, ">>", lhs])
    end

    @assert(#
        all([isassigned(m, i) for i in 1:length(m)]),
        join([#
            "The reactions above must be within the 512 character limit for ",
            "RXNMapper. To overcome this limitation, either: ",
            "(i) remove additional positional information in the SMILES ",
            "string (you may want to check whether this information loss ",
            "changes the atom mapping. ",
            "OR ",
            "(ii) remove this reaction and replace with corresponding sink ",
            "and source reactions carrying the metabolic flux with ",
            "stoichiometric coefficients equal to one."
        ])
    )

    return m
end

function get_first_product_index(m::String)
    subs, _ = split(m, ">>")
    subs = split(subs, ".")
    return length(subs) + 1
end

function trace_atoms(m::String, a::Tuple{Int16, Int16, String})
    # Atom map of mapped SMILES from RXNMapper
    x = canonicalize_and_atom_map(m)

    # Get index of first product
    idx_prods = get_first_product_index(m)

    # Chosen substrate with mapped atoms converted to 1-letter atom codes
    s_atoms = parse_from_2(x[a[1]][1])
    sanitize_atom_smiles(s_atoms, a[3])
    filter!(x -> x != 'H', s_atoms)

    # One letter chemical symbol code
    d = dict_2_to_1_letter_chemical_element()
    if haskey(d, a[3])
        chemical_symbol = d[a[3]]
    else
        chemical_symbol = a[3]
    end

    # Indices of the chosen atom within the chosen substrate
    s_a_idx = findall(==(chemical_symbol[1]), s_atoms)

    # Mapped SMILES substrate index of chosen atom index
    s_idx = x[a[1]][2][s_a_idx][a[2]]

    # Loop over product SMILES
    for i in idx_prods:length(x)
        # Find the product metabolite that the substrate atom moves to
        if s_idx ∈ x[i][2]
            # Mapped product atoms converted to 1-letter codes
            p_atoms = parse_from_2(x[i][1])
            filter!(x -> x != 'H', p_atoms)

            # Indices of chosen atom type in chosen product
            p_a_idx = findall(==(chemical_symbol[1]), p_atoms)

            # Product substrate atom index that maps from substrate
            idx = findfirst(==(s_idx), x[i][2][p_a_idx])

            return Int16(i-idx_prods+1), Int16(idx)
        end
    end
end

# Canonicalize SMILES and atom mappings according to RXNMAPPER
function canonicalize_and_atom_map(m::Vector{String})
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.smiles_utils
    m = [split.(m[i], r"\.|(>>)") for i in eachindex(m)]
    return m .|> x -> rxn_mapper.canonicalize_and_atom_map.(x)
end
function canonicalize_and_atom_map(m::String)
    rxnmap = pyimport("rxnmapper")
    rxn_mapper = rxnmap.smiles_utils
    m = split(m, r"\.|(>>)")
    return rxn_mapper.canonicalize_and_atom_map.(m)
end

# Substrate reaction smiles index from substrate index, stoichiometric copy
function get_srsi(Sj::Vector{Int16}, meti::Int16, copy::Int64)
    sub_indices = findall(<(0), Sj)
    idx = findfirst(==(meti), sub_indices)
    sub_stoichs = abs.(Sj[sub_indices])
    return Int16(sum(sub_stoichs[collect(1:(idx-1))]) + copy)
end

# Product reaction smiles index from product index, stoichiometric copy
function get_prsi(Sj::Vector{Int16}, meti::Int16, copy::Int64)
    prod_indices = findall(>(0), Sj)
    if isempty(prod_indices)
        return Int16(0)
    else
        idx = findfirst(==(meti), prod_indices)
        prod_stoichs = abs.(Sj[prod_indices])
        return Int16(sum(prod_stoichs[collect(1:(idx-1))]) + copy)
    end
end

# Product metabolite from product reaction smiles index
function get_prod_meti(Sj::Vector{Int16}, prod_rsi::Int16)
    prod_indices = findall(>(0), Sj)
    prod_stoichs = Int64.(abs.(Sj[prod_indices]))
    seq = vcat(fill.(prod_indices, prod_stoichs)...)
    prod_meti = seq[prod_rsi]
    return Int16(prod_meti)
end

# Main helper function for the higher order steady_state_distribution (with dictionary)
function construct_atomic_chmc(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    ms::Vector{String},
    I::Tuple{Int64,Int64,String},
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}};
    verbose::Bool = false
)
    # Construct trie to get number of nodes
    d, d′, invD = trie(S, v, ms, I, D)

    # First tumble down the tree to fill queue of prefixes
    if verbose == true
        @info("       Converting CHMC trie to CHMC probability matrix.")
    end
    rs, cs, vs, R, queue = trie_matrix_first_pass(d, [(I[1], I[2])])

    # Traverse remaining prefixes in the queue
    rs2 = Vector{Vector{Int64}}(undef, length(queue))
    cs2 = Vector{Vector{Int64}}(undef, length(queue))
    vs2 = Vector{Vector{Float64}}(undef, length(queue))
    R2 = Vector{#
        Vector{NamedTuple{(:i,:j,:k), Tuple{Int64,Int64,Int16}}}
    }(undef, length(queue))
    Threads.@threads for i in eachindex(queue)
        rs2[i], cs2[i], vs2[i], R2[i] = trie_matrix_standard(d, queue[i])
    end

    # Initialize and return sparse or dense matrix
    if !isempty(rs2)
        rs = [rs; reduce(vcat, rs2)]
        cs = [cs; reduce(vcat, cs2)]
        vs = [vs; reduce(vcat, vs2)]
        R = unique([R; reduce(vcat, R2)])
    end
    T = normalize_rows(sparse(rs, cs, vs, length(d), length(d)))

    return T, d′, invD, R
end

function normalize_rows(A::SparseMatrixCSC)
    A = permutedims(A, [2, 1])
    sums = sum(A, dims = 1)
    I, J, V = findnz(A)
    for idx in 1:length(V)
        V[idx] /= sums[J[idx]]
    end
    return permutedims(sparse(I,J,V), [2, 1])
end

function trie(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    ms::Vector{String},
    I::Tuple{Int64,Int64,String},
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
)
    # First tumble down the tree for key/value pairs
    p, queue = trie_first_pass(S, v, ms, I, D)

    # Traverse remaining prefixes in the queue
    p2 = Vector{Vector{Pair{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Vector{Tuple{Int16,Int16}},
            Vector{Float64},
            Vector{Tuple{Int16,Int16}},
            Vector{Float64},
            Vector{Float64},
            Vector{Int16},
            Vector{Int16},
            Vector{Int16}
        }
    }}}(undef, length(queue))
    Threads.@threads for i in eachindex(queue)
        p2[i] = trie_standard(S, v, ms, I[3], D, queue[i])
    end

    # Concatenate pairs and reshape to dictionary with CHMC state indices
    if !isempty(p2)
        p = [p; reduce(vcat, p2)]
    end

    # Dictionary of MC state to metabolite/atom index
    mc_ids = unique(reduce(vcat, first.(p)))
    invD = Dict(zip(eachindex(mc_ids), mc_ids))

    dd = Dict(zip(mc_ids, eachindex(mc_ids)))
    ks = Vector{Vector{Int16}}(undef, length(p))
    vs = Vector{Vector{Int16}}(undef, length(p))
    Threads.@threads for i in eachindex(ks)
        ks2 = Vector{Int16}(undef, length(first(p[i])))
        vs2 = Vector{Int16}(undef, length(first(last(p[i]))))
        for j in eachindex(first(p[i]))
            ks2[j] = dd[first(p[i])[j]]
        end
        for j in eachindex(first(last(p[i])))
            vs2[j] = dd[first(last(p[i]))[j]]
        end
        ks[i] = ks2
        vs[i] = vs2
    end
    de = Dict(zip(ks, NamedTuple{(:id, :children)}.(tuple.(eachindex(vs), vs))))

    # Add CHMC indices
    d = Dict(zip(#
        first.(p),
        tuple.(eachindex(p), last.(p))
    ))

    return d, de, invD
end

function trie_first_pass(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    ms::Vector{String},
    I::Tuple{Int64,Int64,String},
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
)
    # Initialize dictionary of prefixes with children and reaction indices
    p = Vector{Pair{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Vector{Tuple{Int16,Int16}}, # downstream states
            Vector{Float64},            # downstream fluxes
            Vector{Tuple{Int16,Int16}}, # upstream states
            Vector{Float64},            # upstream fluxes
            Vector{Float64},            # external fluxes
            Vector{Int16},              # downstream reaction indices
            Vector{Int16},              # upstream reaction indices
            Vector{Int16}               # external reaction indices
        }
    }}()
    pfxs = Vector{Vector{Tuple{Int16,Int16}}}()
    function traverse_trie(#
        pfx::Vector{Tuple{Int16,Int16}},
        pfxs::Vector{Vector{Tuple{Int16,Int16}}},
        S::Matrix{Int16},
        v::Vector{<:Real},
        ms::Vector{String},
        p::Vector{Pair{#
            Vector{Tuple{Int16,Int16}},
            Tuple{#
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Float64},
                Vector{Int16},
                Vector{Int16},
                Vector{Int16}
            }
        }},
        a::String,
        D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
    )
        # Get downstream/upstream, external states, fluxes, and reactions
        ds_s = Vector{Tuple{Int16,Int16}}()
        ds_f = Vector{Float64}()
        ex_f = Vector{Float64}()
        ds_r = Vector{Int16}()
        ex_r = Vector{Int16}()
        connecting_r = findall_int16(<(0), @view S[pfx[end][1], :])
        warning = join([#
            "The specified (metabolite index, atom index, ",
            "stoichiometric copy, reaction index) was not ",
            "found. Check that the right input dictionary (D) was specified ",
            "and that the initial starting state (I) is valid."
        ])
        for r in connecting_r
            if !isempty(ms[r])
                for c in 1:abs(S[pfx[end][1], r]) # substrate stoich copies
                    @assert(haskey(D, (pfx[end][1], pfx[end][2], c, r)), warning)
                    if D[(pfx[end][1], pfx[end][2], c, r)] ∉ ds_s
                        push!(ds_s, D[(pfx[end][1], pfx[end][2], c, r)])
                        push!(ds_f, v[r])
                        push!(ds_r, Int16(r))
                    else
                        idx = findfirst(#
                            ==(D[(pfx[end][1], pfx[end][2], c, r)]),
                            ds_s
                        )
                        ds_f[idx] += v[r]
                    end
                end
            else
                push!(ex_f, v[r])
                push!(ex_r, Int16(r))
            end
        end

        # Separate downstream/upstream reactions
        idx = findall(x -> x ∈ pfx, ds_s)
        us_s = ds_s[idx]
        us_f = ds_f[idx]
        deleteat!(ds_s, idx)
        deleteat!(ds_f, idx)
        us_r = ds_r[idx]
        deleteat!(ds_r, idx)

        @assert(#
            !(isempty(us_s) && isempty(ds_s) && isempty(ex_f)),
            join([#
                "Metabolite row index $(pfx[end][1]) does not participate ",
                "as a substrate in any reaction in S."
            ])
        )
        push!(p, pfx => (ds_s, ds_f, us_s, us_f, ex_f, ds_r, us_r, ex_r))
        add_ext_state!(p, pfx, ex_r)

        # Recursively enumerate CHMC and push downstream/upstream to dict
        if !isempty(ds_s)
            for d in ds_s[2:end]
                push!(pfxs, [pfx; d])
            end
            traverse_trie([pfx; ds_s[1]], pfxs, S, v, ms, p, a, D)
        end
    end

    # Construct dictionary of prefix and children and add unique IDs
    traverse_trie([(Int16(I[1]),Int16(I[2]))], pfxs, S, v, ms, p, I[3], D)

    return p, pfxs
end

function add_ext_state!(#
    p::Vector{Pair{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Vector{Tuple{Int16,Int16}},
            Vector{Float64},
            Vector{Tuple{Int16,Int16}},
            Vector{Float64},
            Vector{Float64},
            Vector{Int16},
            Vector{Int16},
            Vector{Int16}
        }
    }},
    pfx::Vector{Tuple{Int16,Int16}},
    ex_r::Vector{Int16}
)
    if !isempty(ex_r)
        v_ext = (Int16(0), Int16(0))
        push!(#
            p,
            [pfx; v_ext] => (#
                Vector{Tuple{Int16,Int16}}(),
                Vector{Float64}(),
                Vector{Tuple{Int16,Int16}}(),
                Vector{Float64}(),
                Vector{Float64}(),
                Vector{Int16}(),
                Vector{Int16}(),
                Vector{Int16}()
            )
        )
    end
end

function trie_standard(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    ms::Vector{String},
    a::String,
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}},
    pfx::Vector{Tuple{Int16,Int16}}
)
    # Initialize dictionary of prefixes with children and reaction indices
    p = Vector{Pair{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Vector{Tuple{Int16,Int16}}, # downstream states
            Vector{Float64},            # downstream fluxes
            Vector{Tuple{Int16,Int16}}, # upstream states
            Vector{Float64},            # upstream fluxes
            Vector{Float64},            # external fluxes
            Vector{Int16},              # downstream reactions
            Vector{Int16},              # upstream reactions
            Vector{Int16}               # external reactions
        }
    }}()
    function traverse_trie(#
        pfx::Vector{Tuple{Int16,Int16}},
        S::Matrix{Int16},
        v::Vector{<:Real},
        p::Vector{Pair{#
            Vector{Tuple{Int16,Int16}},
            Tuple{#
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Float64},
                Vector{Int16},
                Vector{Int16},
                Vector{Int16}
            }
        }},
        a::String,
        D::Dict{NTuple{4,Int64},Tuple{Int64,Int64}}
    )
        # Get downstream/upstream, external states, fluxes, and reactions
        ds_s = Vector{Tuple{Int16,Int16}}()
        ds_f = Vector{Float64}()
        ex_f = Vector{Float64}()
        ds_r = Vector{Int16}()
        ex_r = Vector{Int16}()
        connecting_r = findall_int16(<(0), @view S[pfx[end][1], :])
        for r in connecting_r
            if !isempty(ms[r])
                for c in 1:abs(S[pfx[end][1], r]) # substrate stoich copies
                    if D[(pfx[end][1], pfx[end][2], c, r)] ∉ ds_s
                        push!(ds_s, D[(pfx[end][1], pfx[end][2], c, r)])
                        push!(ds_f, v[r])
                        push!(ds_r, Int16(r))
                    else
                        idx = findfirst(#
                            ==(D[(pfx[end][1], pfx[end][2], c, r)]),
                            ds_s
                        )
                        ds_f[idx] += v[r]
                    end
                end
            else
                push!(ex_f, v[r])
                push!(ex_r, Int16(r))
            end
        end

        # Separate downstream/upstream reactions
        idx = findall(x -> x ∈ pfx, ds_s)
        us_s = ds_s[idx]
        us_f = ds_f[idx]
        deleteat!(ds_s, idx)
        deleteat!(ds_f, idx)
        us_r = ds_r[idx]
        deleteat!(ds_r, idx)

        @assert(#
            !(isempty(us_s) && isempty(ds_s) && isempty(ex_f)),
            join([#
                "Metabolite row index $(pfx[end][1]) does not participate ",
                "as a substrate in any reaction in S."
            ])
        )
        push!(p, pfx => (ds_s, ds_f, us_s, us_f, ex_f, ds_r, us_r, ex_r))
        add_ext_state!(p, pfx, ex_r)

        # Recursively enumerate CHMC and push downstream/upstream to dict
        for s in ds_s
            traverse_trie([pfx; s], S, v, p, a, D)
        end
    end

    # Construct dictionary of prefix and children and add unique IDs
    traverse_trie(pfx, S, v, p, a, D)

    return p
end

function trie_matrix_first_pass(#
    d::Dict{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Int64,
            Tuple{#
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Float64},
                Vector{Int16},
                Vector{Int16},
                Vector{Int16}
        }
        }
    },
    pfx::Vector{Tuple{Int64,Int64}}
)
    # Construct (sparse) CHMC transition matrix with COO format
    rt = Vector{Int64}()
    ct = Vector{Int64}()
    vt = Vector{Float64}()

    # Construct reaction index tuples (row i in T, col j in T, value T[i,j])
    R = Vector{NamedTuple{(:i,:j,:k),Tuple{Int64,Int64,Int16}}}()

    # Queue
    pfxs = Vector{Vector{Tuple{Int64,Int64}}}()

    function traverse_trie(#
        pfx::Vector{Tuple{Int64,Int64}},
        d::Dict{#
            Vector{Tuple{Int16,Int16}},
            Tuple{#
                Int64,
                Tuple{#
                    Vector{Tuple{Int16,Int16}},
                    Vector{Float64},
                    Vector{Tuple{Int16,Int16}},
                    Vector{Float64},
                    Vector{Float64},
                    Vector{Int16},
                    Vector{Int16},
                    Vector{Int16}
                }
            }
        },
        rt::Vector{Int64}, ct::Vector{Int64}, vt::Vector{Float64},
        R::Vector{NamedTuple{(:i,:j,:k),Tuple{Int64,Int64,Int16}}},
        pfxs::Vector{Vector{Tuple{Int64,Int64}}}
    )
        # Children, parent, external indices/reactions
        chd_s, chd_f, par_s, par_f, ext_f, chd_r, par_r, ext_r = last(d[pfx])
        i = first(d[pfx])

        # Fill matrix with external transitions
        for k in eachindex(ext_f)
            ii = first(d[[pfx; (0,0)]])
            push!(rt, i)
            push!(ct, ii)
            push!(vt, ext_f[k])

            push!(rt, ii)
            push!(ct, 1)
            push!(vt, ext_f[k])

            push!(R, (i = i,  j = ii, k = ext_r[k]))
            push!(R, (i = ii, j = 1,  k = ext_r[k]))
        end

        # Fill matrix with parent transitions
        for k in eachindex(par_s)
            j = first(d[pfx[1:findfirst(pfx .== Ref(par_s[k]))]])
            push!(rt, i)
            push!(ct, j)
            push!(vt, par_f[k])
            push!(R, (i = i, j = j, k = par_r[k]))
        end

        # Fill matrix with children transitions
        for k in eachindex(chd_s)
            j = first(d[[pfx; chd_s[k]]])
            push!(rt, i)
            push!(ct, j)
            push!(vt, chd_f[k])
            push!(R, (i = i, j = j, k = chd_r[k]))
            push!(pfxs, [pfx; chd_s[k]])
        end
        if !isempty(chd_s)
            p = (Int64(first(first(chd_s))), Int64(last(first(chd_s))))
            traverse_trie([pfx; p], d, rt, ct, vt, R, pfxs)
        end
    end
    traverse_trie(pfx, d, rt, ct, vt, R, pfxs)

    return rt, ct, vt, R, pfxs
end

function trie_matrix_standard(#
    d::Dict{#
        Vector{Tuple{Int16,Int16}},
        Tuple{#
            Int64,
            Tuple{#
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Tuple{Int16,Int16}},
                Vector{Float64},
                Vector{Float64},
                Vector{Int16},
                Vector{Int16},
                Vector{Int16}
            }
        }
    },
    pfx::Vector{Tuple{Int64,Int64}}
)
    # Construct (sparse) CHMC transition matrix with COO format
    rt = Vector{Int64}()
    ct = Vector{Int64}()
    vt = Vector{Float64}()

    # Construct reaction index tuples (row i in T, col j in T, value T[i,j])
    R = Vector{NamedTuple{(:i,:j,:k), Tuple{Int64,Int64,Int16}}}()

    function traverse_trie(#
        pfx::Vector{Tuple{Int64,Int64}},
        d::Dict{#
            Vector{Tuple{Int16,Int16}},
            Tuple{#
                Int64,
                Tuple{#
                    Vector{Tuple{Int16,Int16}},
                    Vector{Float64},
                    Vector{Tuple{Int16,Int16}},
                    Vector{Float64},
                    Vector{Float64},
                    Vector{Int16},
                    Vector{Int16},
                    Vector{Int16}
                }
            }
        },
        rt::Vector{Int64}, ct::Vector{Int64}, vt::Vector{Float64},
        R::Vector{NamedTuple{(:i,:j,:k),Tuple{Int64,Int64,Int16}}}
    )
        # Children, parent, external indices/reactions
        chd_s, chd_f, par_s, par_f, ext_f, chd_r, par_r, ext_r = last(d[pfx])
        i = first(d[pfx])

        # Fill matrix with external transitions
        for k in eachindex(ext_f)
            ii = first(d[[pfx; (0,0)]])
            push!(rt, i)
            push!(ct, ii)
            push!(vt, ext_f[k])
            push!(rt, ii)
            push!(ct, 1)
            push!(vt, ext_f[k])
            push!(R, (i = i,  j = ii, k = ext_r[k]))
            push!(R, (i = ii, j = 1,  k = ext_r[k]))
        end

        # Fill matrix with parent transitions
        for k in eachindex(par_s)
            j = first(d[pfx[1:findfirst(pfx .== Ref(par_s[k]))]])
            push!(rt, i)
            push!(ct, j)
            push!(vt, par_f[k])
            push!(R, (i = i, j = j, k = par_r[k]))
        end

        # Fill matrix with children transitions
        for k in eachindex(chd_s)
            j = first(d[[pfx; chd_s[k]]])
            push!(rt, i)
            push!(ct, j)
            push!(vt, chd_f[k])
            push!(R, (i = i, j = j, k = chd_r[k]))
            p = (Int64(first(chd_s[k])), Int64(last(chd_s[k])))
            traverse_trie([pfx; p], d, rt, ct, vt, R)
        end
    end
    traverse_trie(pfx, d, rt, ct, vt, R)

    return rt, ct, vt, R
end

# Main helper enumerate_efms
function construct_atomic_chmc_old(#
  S::Matrix{<:Real},
  ms::Vector{String},
  I::Tuple{Int64,Int64,String},
  D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
)
  # Error checking
  sanitize_stoichiometry(S)

  # Initialize dictionary of prefix and (children, reaction index)
  d = Dict{#
    Vector{Tuple{Int16,Int16}},
    Tuple{#
      Vector{Tuple{Int16,Int16}},
      Vector{Float64},
      Vector{Tuple{Int16,Int16}},
      Vector{Float64},
      Vector{Float64},
      Vector{Int16},
      Vector{Int16},
      Vector{Int16}
    }
  }()
  function traverse_trie(#
    prefix::Vector{Tuple{Int16,Int16}},
    S::Matrix{Int16},
    d::Dict{#
      Vector{Tuple{Int16,Int16}},
      Tuple{#
        Vector{Tuple{Int16,Int16}},
        Vector{Float64},
        Vector{Tuple{Int16,Int16}},
        Vector{Float64},
        Vector{Float64},
        Vector{Int16},
        Vector{Int16},
        Vector{Int16}
      }
    },
    a::String,
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
  )

    # Get downstream/upstream, external states and reactions
    ds_states = Vector{Tuple{Int16,Int16}}()
    ds_fluxes = Vector{Float64}()
    ex_fluxes = Vector{Float64}()
    ds_rxns = Vector{Int16}()
    ex_rxns = Vector{Int16}()
    connecting_rxns = findall_int16(<(0), @view S[prefix[end][1],:])
    for r in connecting_rxns
      if !isempty(ms[r])
        for c in 1:abs(S[prefix[end][1],r]) # substrate stoichiometry copies

          if D[(prefix[end][1], prefix[end][2], c, r)] ∉ ds_states
            push!(ds_states, D[(prefix[end][1], prefix[end][2], c, r)])
            push!(ds_fluxes, 1.0)
            push!(ds_rxns, Int16(r))
          else
            idx = findfirst(#
              ==(D[(prefix[end][1], prefix[end][2], c, r)]),
              ds_states
            )
            ds_fluxes[idx] += 1.0
          end
        end
      else
        push!(ex_fluxes, 1.0)
        push!(ex_rxns, Int16(r))
      end
    end

    # Separate downstream/upstream reactions
    idx = findall(x -> x ∈ prefix, ds_states)
    us_states = ds_states[idx]
    us_fluxes = ds_fluxes[idx]
    deleteat!(ds_states, idx)
    deleteat!(ds_fluxes, idx)
    us_rxns = ds_rxns[idx]
    deleteat!(ds_rxns, idx)

    @assert(#
      !(isempty(us_states) && isempty(ds_states) && isempty(ex_fluxes)),
      join([#
        "Metabolite row index $(prefix[end][1]) does not participate as a ",
        "substrate in any reaction in S."
      ])
    )
    d[prefix] = (ds_states, ds_fluxes, us_states, us_fluxes, ex_fluxes, ds_rxns, us_rxns, ex_rxns)

    # Recursively enumerate CHMC and push downstream/upstream to dict
    for s in ds_states
      #@info("State space: $(length(d))")
      traverse_trie([prefix; s], S, d, a, D)
    end
  end

  # Construct dictionary of prefix and children and add unique IDs
  traverse_trie([(Int16(I[1]),Int16(I[2]))], S, d, I[3], D)

  # Add external environment metabolites to the dictionary
  idx = .!isempty.(last.(collect(values(d))))
  ks  = deepcopy(collect(keys(d))[idx])
  for k in ks
    push!(k, (0, 0))
  end
  vs = Tuple{#
    Vector{Tuple{Int16,Int16}},
    Vector{Float64},
    Vector{Tuple{Int16,Int16}},
    Vector{Float64},
    Vector{Float64},
    Vector{Int16},
    Vector{Int16},
    Vector{Int16}
  }[]
  for i in 1:length(ks)
    push!(vs, ([], [], [], [], [], [], [], []))
  end
  d_ext = Dict(zip(ks, vs))
  d = merge(d, d_ext)

  # Ensure root node/prefix is variable 1
  v1 = ["X" * string(i) for i in 1:length(keys(d))]
  v2 = collect(values(d))
  root_idx = findfirst(==(1), length.(collect(keys(d))))
  switch_idx = findfirst(v1 .== "X1")
  v1[switch_idx] = v1[root_idx]
  v1[root_idx] = "X1"
  vs = [(v1[i], v2[i]) for i in 1:length(v1)]
  d = Dict(zip(keys(d), vs))

  # Construct CHMC transition probability matrix
  T = zeros(length(d), length(d)) # transition probabilities
  R = zeros(Int16, length(d), length(d)) # reaction indices

  # Traverse tree and populate new transition probability matrix with trie nodes
  function fill_chmc_matrix(#
    prefix::Vector{Tuple{Int16,Int16}},
    d::Dict{#
      Vector{Tuple{Int16,Int16}},
      Tuple{#
        String,
        Tuple{#
          Vector{Tuple{Int16,Int16}},
          Vector{Float64},
          Vector{Tuple{Int16,Int16}},
          Vector{Float64},
          Vector{Float64},
          Vector{Int16},
          Vector{Int16},
          Vector{Int16}
        }
      }
    },
    S::Matrix{Int16},
    T::Matrix{Float64},
    R::Matrix{Int16}
  )
    # Children, parent, external indices/reactions
    childs, child_fluxes, pars, par_fluxes, ext_fluxes, child_rxns, par_rxns, ext_rxns = last(d[prefix])
    i = parse(Int16, first(d[prefix])[2:end]) # current state

    # Fill matrix with external transitions
    for k in 1:length(ext_fluxes)
      ii = parse(Int16, first(d[[prefix; (0,0)]])[2:end])
      T[i,ii] = ext_fluxes[k] # to external environment state
      T[ii,1] = ext_fluxes[k] # External environment state back to root
      R[i,ii] = ext_rxns[k]
      R[ii,1] = ext_rxns[k]
    end

    # Fill matrix with parent transitions
    for k in 1:length(pars)
      idx = findfirst(prefix .== Ref(pars[k]))
      j = parse(Int16, first(d[prefix[1:idx]])[2:end])
      T[i,j] = par_fluxes[k]
      R[i,j] = par_rxns[k]
    end

    # Fill matrix with children transitions
    for k in 1:length(childs)
      j = parse(Int16, first(d[[prefix; childs[k]]])[2:end])
      T[i,j] = child_fluxes[k]
      R[i,j] = child_rxns[k]
      fill_chmc_matrix([prefix; childs[k]], d, S, T, R)
    end
  end
  fill_chmc_matrix([(Int16(I[1]), Int16(I[2]))], d, S, T, R)

  # Normalize flux adjacency matrix to transition probability matrix
  T = T ./ sum(T, dims=2)
  @assert(all(.!isnan.(T)), "")
  @assert(all(T .>= 0))

  # Dictionary relating single-valued states to (met row, atom formula index)
  ids = unique(vcat(collect(keys(d))...))
  invD = Dict(zip(1:length(ids), ids))

  # Dictionary with just keys/children values by unique IDs
  dd = Dict(zip(ids, 1:length(ids)))
  ks = Vector{Vector{Int16}}(undef, length(d)) # keys
  key_ids = collect(keys(d))
  for i in 1:length(key_ids)
    tmp = Vector{Int16}(undef, length(key_ids[i]))
    for j in 1:length(key_ids[i])
      tmp[j] = dd[key_ids[i][j]]
    end
    ks[i] = tmp
  end
  vs = Vector{Vector{Int16}}(undef, length(d)) # values
  val_ids = first.(last.(collect(values(d))))
  for i in 1:length(val_ids)
    tmp = Vector{Int16}(undef, length(val_ids[i]))
    for j in 1:length(val_ids[i])
      tmp[j] = dd[val_ids[i][j]]
    end
    vs[i] = tmp
  end
  # Combine keys, values, and additional unique prefix IDs
  ids = first.(collect(values(d)))
  vs = [(id=ids[i], children=vs[i]) for i in 1:length(ids)]
  d′ = Dict(zip(ks, vs))

  return T, d′, invD, R
end
function construct_atomic_chmc_sparse_old(#
  S::Matrix{<:Real},
  ms::Vector{String},
  I::Tuple{Int64,Int64,String},
  D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
)
  # Error checking
  sanitize_stoichiometry(S)

  # Initialize dictionary of prefix and (children, reaction index)
  d = Dict{#
    Vector{Tuple{Int16,Int16}},
    Tuple{#
      Vector{Tuple{Int16,Int16}},
      Vector{Float64},
      Vector{Tuple{Int16,Int16}},
      Vector{Float64},
      Vector{Float64},
      Vector{Int16},
      Vector{Int16},
      Vector{Int16}
    }
  }()
  function traverse_trie(#
    prefix::Vector{Tuple{Int16,Int16}},
    S::Matrix{Int16},
    d::Dict{#
      Vector{Tuple{Int16,Int16}},
      Tuple{#
        Vector{Tuple{Int16,Int16}},
        Vector{Float64},
        Vector{Tuple{Int16,Int16}},
        Vector{Float64},
        Vector{Float64},
        Vector{Int16},
        Vector{Int16},
        Vector{Int16}
      }
    },
    a::String,
    D::Dict{NTuple{4,Int64}, Tuple{Int64,Int64}}
  )

    # Get downstream/upstream, external states and reactions
    ds_states = Vector{Tuple{Int16,Int16}}()
    ds_fluxes = Vector{Float64}()
    ex_fluxes = Vector{Float64}()
    ds_rxns = Vector{Int16}()
    ex_rxns = Vector{Int16}()
    connecting_rxns = findall_int16(<(0), @view S[prefix[end][1],:])
    for r in connecting_rxns
      if !isempty(ms[r])
        for c in 1:abs(S[prefix[end][1],r]) # substrate stoichiometry copies

          if D[(prefix[end][1], prefix[end][2], c, r)] ∉ ds_states
            push!(ds_states, D[(prefix[end][1], prefix[end][2], c, r)])
            push!(ds_fluxes, 1.0)
            push!(ds_rxns, Int16(r))
          else
            idx = findfirst(#
              ==(D[(prefix[end][1], prefix[end][2], c, r)]),
              ds_states
            )
            ds_fluxes[idx] += 1.0
          end
        end
      else
        push!(ex_fluxes, 1.0)
        push!(ex_rxns, Int16(r))
      end
    end

    # Separate downstream/upstream reactions
    idx = findall(x -> x ∈ prefix, ds_states)
    us_states = ds_states[idx]
    us_fluxes = ds_fluxes[idx]
    deleteat!(ds_states, idx)
    deleteat!(ds_fluxes, idx)
    us_rxns = ds_rxns[idx]
    deleteat!(ds_rxns, idx)

    @assert(#
      !(isempty(us_states) && isempty(ds_states) && isempty(ex_fluxes)),
      join([#
        "Metabolite row index $(prefix[end][1]) does not participate as a ",
        "substrate in any reaction in S."
      ])
    )
    d[prefix] = (ds_states, ds_fluxes, us_states, us_fluxes, ex_fluxes, ds_rxns, us_rxns, ex_rxns)

    # Recursively enumerate CHMC and push downstream/upstream to dict
    for s in ds_states
      @info("State space: $(length(d))")
      traverse_trie([prefix; s], S, d, a, D)
    end
  end

  # Construct dictionary of prefix and children and add unique IDs
  traverse_trie([(Int16(I[1]),Int16(I[2]))], S, d, I[3], D)

  # Add external environment metabolites to the dictionary
  @info("Adding external environment metabolite to CHMC tree")
  idx = .!isempty.(last.(collect(values(d))))
  ks  = deepcopy(collect(keys(d))[idx])
  for k in ks
    push!(k, (0, 0))
  end
  vs = Tuple{#
    Vector{Tuple{Int16,Int16}},
    Vector{Float64},
    Vector{Tuple{Int16,Int16}},
    Vector{Float64},
    Vector{Float64},
    Vector{Int16},
    Vector{Int16},
    Vector{Int16}
  }[]
  for i in 1:length(ks)
    push!(vs, ([], [], [], [], [], [], [], []))
  end
  d_ext = Dict(zip(ks, vs))
  d = merge(d, d_ext)

  # Ensure root node/prefix is variable 1
  @info("Setting root node to index 1.")
  v1 = ["X" * string(i) for i in 1:length(keys(d))]
  v2 = collect(values(d))
  root_idx = findfirst(==(1), length.(collect(keys(d))))
  switch_idx = findfirst(v1 .== "X1")
  v1[switch_idx] = v1[root_idx]
  v1[root_idx] = "X1"
  vs = [(v1[i], v2[i]) for i in 1:length(v1)]
  d = Dict(zip(keys(d), vs))

  # Construct CHMC transition probability matrix
  @info("Constructing CHMC transition probability matrix.")
  T = spzeros(Float64, length(d), length(d)) # transition probabilities
  R = spzeros(Int16, length(d), length(d)) # reaction indices

  # Traverse tree and populate new transition probability matrix with trie nodes
  function fill_chmc_matrix(#
    prefix::Vector{Tuple{Int64,Int64}},
    d::Dict{#
      Vector{Tuple{Int16,Int16}},
      Tuple{#
        String,
        Tuple{#
          Vector{Tuple{Int16,Int16}},
          Vector{Float64},
          Vector{Tuple{Int16,Int16}},
          Vector{Float64},
          Vector{Float64},
          Vector{Int16},
          Vector{Int16},
          Vector{Int16}
        }
      }
    },
    S::Matrix{Int16},
    T::SparseMatrixCSC,
    R::SparseMatrixCSC,
  )
    # Children, parent, external indices/reactions
    childs, child_fluxes, pars, par_fluxes, ext_fluxes, child_rxns, par_rxns, ext_rxns = last(d[prefix])
    i = parse(Int64, first(d[prefix])[2:end]) # current state

    # Fill matrix with external transitions
    for k in 1:length(ext_fluxes)
      ii = parse(Int64, first(d[[prefix; (0,0)]])[2:end])
      T[i,ii] = ext_fluxes[k] # to external environment state
      T[ii,1] = ext_fluxes[k] # External environment state back to root
      R[i,ii] = ext_rxns[k]
      R[ii,1] = ext_rxns[k]
    end

    # Fill matrix with parent transitions
    for k in 1:length(pars)
      idx = findfirst(prefix .== Ref(pars[k]))
      j = parse(Int64, first(d[prefix[1:idx]])[2:end])
      T[i,j] = par_fluxes[k]
      R[i,j] = par_rxns[k]
    end

    # Fill matrix with children transitions
    for k in 1:length(childs)
      j = parse(Int64, first(d[[prefix; childs[k]]])[2:end])
      T[i,j] = child_fluxes[k]
      R[i,j] = child_rxns[k]
      fill_chmc_matrix([prefix; (Int64(childs[k][1]), Int64(childs[k][2]))], d, S, T, R)
    end
  end
  fill_chmc_matrix([(Int64(I[1]), Int64(I[2]))], d, S, T, R)

  # Normalize flux adjacency matrix to transition probability matrix
  @info("Normalizing transition probability matrix")
  #T = T ./ sum(T, dims=2)
  function normalize_rows(A::SparseMatrixCSC)
    A = permutedims(A, [2, 1])
    sums = sum(A, dims=1)
    I,J,V = findnz(A)
    for idx in 1:length(V)
      V[idx] /= sums[J[idx]]
    end
    return permutedims(sparse(I,J,V), [2, 1])
  end
  T = normalize_rows(T)
  #@assert(all(.!isnan.(T)), "")
  #@assert(all(T .>= 0))
  @info("Completed CHMC transition probability matrix.")

  # Dictionary relating single-valued states to (met row, atom formula index)
  @info("Constructing dictionary of atomic CHMC states.")
  ids = unique(vcat(collect(keys(d))...))
  invD = Dict(zip(1:length(ids), ids))

  # Dictionary with just keys/children values by unique IDs
  dd = Dict(zip(ids, 1:length(ids)))
  ks = Vector{Vector{Int16}}(undef, length(d)) # keys
  key_ids = collect(keys(d))
  for i in 1:length(key_ids)
    tmp = Vector{Int16}(undef, length(key_ids[i]))
    for j in 1:length(key_ids[i])
      tmp[j] = dd[key_ids[i][j]]
    end
    ks[i] = tmp
  end
  vs = Vector{Vector{Int16}}(undef, length(d)) # values
  val_ids = first.(last.(collect(values(d))))
  for i in 1:length(val_ids)
    tmp = Vector{Int16}(undef, length(val_ids[i]))
    for j in 1:length(val_ids[i])
      tmp[j] = dd[val_ids[i][j]]
    end
    vs[i] = tmp
  end
  # Combine keys, values, and additional unique prefix IDs
  ids = first.(collect(values(d)))
  vs = [(id=ids[i], children=vs[i]) for i in 1:length(ids)]
  d′ = Dict(zip(ks, vs))

  return T, d′, invD, R
end

# Print formatting of ErrorSummary struct
function string_errorsummary(s::CHMCAtomicErrorSummary)
    g(x) = join(x, ", ")
    h(x) = isempty(g(x)) ? "NONE." : g(x)
    arr = [#
        string(s.absolute_flux_error),
        h(s.reactions_duplicated),
        h(s.reactions_with_zero_flux),
        h(s.reactions_with_negative_flux),
        h(s.reactions_with_non_integer_stoichiometries),
        h(s.reactions_unimolecular_source_stoichiometry_one),
        h(s.reactions_unimolecular_source_stoichiometry_not_one),
        h(s.reactions_multimolecular_source_stoichiometry_one),
        h(s.reactions_multimolecular_source_stoichiometry_not_one),
        h(s.reactions_unimolecular_sink_stoichiometry_one),
        h(s.reactions_unimolecular_sink_stoichiometry_not_one),
        h(s.reactions_multimolecular_sink_stoichiometry_one),
        h(s.reactions_multimolecular_sink_stoichiometry_not_one),
        h(s.reactions_empty),
        h(s.unused_metabolites),
    ]
    return arr
end

# Helper for printing function
function pass_fail(res::CHMCAtomicErrorSummary)
    bool = repeat(["PASSED."], length(fieldnames(CHMCAtomicErrorSummary)))
    f = "FAILED."
    if res.absolute_flux_error > sqrt(eps())
        bool[1] = f
    end
    if !isempty(res.reactions_duplicated)
        bool[2] = f
    end
    if !isempty(res.reactions_with_zero_flux)
        bool[3] = f
    end
    if !isempty(res.reactions_with_negative_flux)
        bool[4] = f
    end
    if !isempty(res.reactions_with_non_integer_stoichiometries)
        bool[5] = f
    end
    if isempty(res.reactions_unimolecular_source_stoichiometry_one)
        bool[6] = f
    end
    if !isempty(res.reactions_unimolecular_source_stoichiometry_not_one)
        bool[7] = f
    end
    if !isempty(res.reactions_multimolecular_source_stoichiometry_one)
        bool[8] = f
    end
    if !isempty(res.reactions_multimolecular_source_stoichiometry_not_one)
        bool[9] = f
    end
    if isempty(res.reactions_unimolecular_sink_stoichiometry_one)
        bool[10] = f
    end
    if !isempty(res.reactions_unimolecular_sink_stoichiometry_not_one)
        bool[11] = f
    end
    if !isempty(res.reactions_multimolecular_sink_stoichiometry_one)
        bool[12] = f
    end
    if !isempty(res.reactions_multimolecular_sink_stoichiometry_not_one)
        bool[13] = f
    end
    if !isempty(res.reactions_empty)
        bool[14] = f
    end
    if !isempty(res.unused_metabolites)
        bool[15] = f
    end
    return bool
end

# Remove internal reactions and replace with source/sink fluxes
function convert_internal_rxn_to_sources_sinks(#
    S::Matrix{<:Real},
    v::Vector{<:Real},
    rxns::Vector{String},
    jj::Vector{Int64}
)

    S_concat = Any[]
    v_concat = Vector{Float64}()
    rxns_concat = Vector{String}()
    for j in jj
        Sj = S[:, j]
        neg_idx = findall(<(0), Sj)
        pos_idx = findall(>(0), Sj)
        neg_coeff = abs.(Sj[neg_idx])
        pos_coeff = abs.(Sj[pos_idx])

        # Repeat sink/source reaction names
        rink = Vector{String}()
        for i in eachindex(neg_idx)
            push!(rink, join([rxns[j], "_sink_met_", neg_idx[i]]))
        end
        rour = Vector{String}()
        for i in eachindex(pos_idx)
            push!(rour, join([rxns[j], "_source_met_", pos_idx[i]]))
        end

        # Construct sink reactions and assign fluxes
        Sink = zeros(Float64, size(S,1), length(neg_idx))
        vink = Vector{Float64}(undef, length(neg_idx))
        for i in 1:length(neg_idx)
            Sink[neg_idx[i],i] = -1
            vink[i] = neg_coeff[i] * v[j]
        end

        # Construct source reactions
        Sour = zeros(Float64, size(S, 1), length(pos_idx))
        vour = Vector{Float64}(undef, length(pos_idx))
        for i in eachindex(pos_idx)
            Sour[pos_idx[i], i] = +1
            vour[i] = pos_coeff[i] * v[j]
        end

        push!(S_concat, Sink)
        push!(S_concat, Sour)
        append!(v_concat, vink)
        append!(v_concat, vour)
        append!(rxns_concat, rink)
        append!(rxns_concat, rour)
    end

    if isempty(S_concat)
        S_concat = Matrix{Float64}(undef, size(S, 1), 0)
    else
        S_concat = reduce(hcat, S_concat)
    end
    return S_concat, v_concat, rxns_concat
end

# Compress linear pathways in the CHMC matrix
function compress_chmc_matrix(transition_matrix)
    # Identify states with a single incoming and outgoing transition
    n = size(transition_matrix, 1)
    middle_states = Int64[]
    root = 1
    for state in 1:n
        cond1 = length(findall(x -> x > 0, transition_matrix[state, :])) == 1
        cond2 = length(findall(x -> x > 0, transition_matrix[:, state])) == 1
        cond3 = findfirst(x -> x > 0, transition_matrix[state, :]) != root
        if cond1 && cond2 && cond3
            push!(middle_states, state)
        end
    end

    # Store middle states between each (from, to) pair
    middle_states_removed = Dict{Tuple{Int16, Int16}, Vector{Int16}}()

    # Update transition probabilities
    for middle in middle_states
        from = findfirst(x -> x > 0, transition_matrix[:, middle])
        to = findfirst(x -> x > 0, transition_matrix[middle, :])
        if from != to
            transition_matrix[from, to] = transition_matrix[from, middle]
            transition_matrix[from, middle] = 0
            transition_matrix[middle, to] = 0

            # Record middle states for each (from, to) pair
            if !haskey(middle_states_removed, (from, to))
                middle_states_removed[(from, middle)] = [to]
            else
                push!(middle_states_removed[(from, middle)], to)
            end
        end
    end

    # Assemble removed middle states into sequences
    function ismember(x::Pair{Tuple{Int16,Int16}, Vector{Int16}}, y::Vector{Int16})
        if x[1][1] == y[1] && x[1][2] == y[end]
            return true
        else
            return false
        end
    end
    assembly = Vector{Vector{Int16}}()
    while !isempty(middle_states_removed)
        x = pop!(middle_states_removed)
        i = findfirst(a -> ismember(x, a), assembly)
        if isnothing(i)
            push!(assembly, [x[1][1], x[1][2], x[2][1]])
        else
            push!(assembly[i], x[2][1])
        end
    end
    return transition_matrix, assembly
end
function compress(T::SparseMatrixCSC)
    return compress_chmc_matrix(T)
end
function compress(T::ExtendableSparseMatrix)
    T, U = compress_chmc_matrix(T)
    flush!(T)
    dropzeros!(T)
    return T, U
end


