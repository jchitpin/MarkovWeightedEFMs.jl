## Helpers

# Default panel image
function default_image()
    Drawing(900, 300, :svg)
    origin()
    background("white")
    sethue("black")
    fontsize(50)
    text("Click on a CHMC transition", halign=:center, valign=:center)
    finish()
    s = svgstring()
    return s
end

# Convert an SVG string to RGB matrix
function svg_to_rgb_matrix(svg::String)
    # Luxor
    svg_image = readsvg(svg)
    w = svg_image.width
    h = svg_image.height

    Drawing(w, h, :png)
    placeimage(svg_image)
    mat = image_as_matrix()
    finish()
    return mat
end

# Get atom to highlight
function rdkit_atom_to_highlight(#
    s::String,
    atom_of_interest::Tuple{Char, String}
)

    s = replace(s, "c" => "C")
    s = replace(s, "n" => "N")
    s = replace(s, "o" => "O")
    s = collect(s)

    filter!(isletter, s)
    return findall(s .== atom_of_interest[2])[atom_of_interest[1]] - 1
end

# Constructing reaction SMILES showing the highlighted atoms from LHS to RHS
function visualize_atomic_chmc_transition(#
    res::CHMCAtomicSummary,
    mets::Vector{String},
    S::Matrix{Int16},
    rxn::Int64,
    A::Matrix{Float64},
    node_labels::Vector{String},
    rss::String,
    idx::Int64;
    canvas_height::Int64=300,
    canvas_width::Int64=500
)

    # Reaction/transition selected in the atomic CHMC
    #rxn = parse(Int64, split(p.elabels[][idx], " ")[2][2:(end-1)])

    # Graph labels flanking selected reaction
    cis = findall(>(0), A)
    cis = getindex.(cis, [1 2])
    sorted_cis = cis[sortperm(cis[:,1]),:]
    a = sorted_cis[idx,:] # graph nodes R[a[1],a[2]]

    # CHMC states flanking selected reaction
    chmc_sub = parse(Int64, node_labels[a[1]])
    chmc_prod = parse(Int64, node_labels[a[2]])

    # MC states flanking selected reaction
    #chmc_ks = parse.(Int64, [i[2:end] for i in first.(collect(values(res.dchmc)))])
    chmc_ks = first.(collect(values(res.dchmc)))
    chmc_vs = collect(keys(res.dchmc))
    mc_sub = chmc_vs[findfirst(chmc_ks .== chmc_sub)][end]
    mc_prod = chmc_vs[findfirst(chmc_ks .== chmc_prod)][end]

    # Substrate/product indices flanking selected reaction
    sstate = res.dmc[mc_sub]
    pstate = res.dmc[mc_prod]

    # Substrate/product reaction string indices
    sub_idx = get_srsi(S[:,rxn], sstate[1], 1)
    prod_idx = get_prsi(S[:,rxn], pstate[1], 1)

    # Split reaction string
    subs, prods = split(rss, ">>")
    subss = split(subs, ".")
    if isempty(prods)
        prodss = Vector{String}()
    else
        prodss = split(prods, ".")
    end

  # Metabolite names ordered from LHS to RHS as they occur in S
    function met_names(Sj::Vector{Int16}, mets::Vector{String})
        idx = findall(<(0), Sj)
        sub_stoich = abs.(Sj[idx])
        subs = vcat(repeat.([[m] for m in mets[idx]], sub_stoich)...)
        idx = findall(>(0), Sj)
        prod_stoich = Sj[idx]
        prods = vcat(repeat.([[m] for m in mets[idx]], prod_stoich)...)
        return [subs; prods]
    end
    metnames = met_names(S[:,rxn], mets)

    # Initialize SVGs each metabolite in the reaction string
    svg_subs = repeat([""], length(subss))
    svg_prods = repeat([""], length(prodss))

    function atom_to_highlight(s::String, atom::Char, i::Int16)
        s = replace(s, "c" => "C")
        s = replace(s, "n" => "N")
        s = replace(s, "o" => "O")
        s = collect(s)
        filter!(isletter, s)
        filter!(x -> x != 'H', s)
        return findall(s .== atom)[i] - 1
    end

    # Construct substrate SVGs and highlight atom of interest
    for i in 1:length(svg_subs)
        args = Dict{String, Any}(#
            "height" => canvas_height,
            "width" => canvas_width,
            "highlightColour" => [137/256, 49/256, 239/256, 0.9]
        )
        if i == sub_idx
            hl = atom_to_highlight(string(subss[i]), only(string(res.i[3])), sstate[2])
            args["atoms"] = [hl]
        else
            if haskey(args, "atoms")
                delete!(args, "atoms")
            end
        end
        svg_subs[i] = get_svg(get_mol(subss[i]), args)
    end

    # Construct produce SVGs and highlight atom of interest
    for i in 1:length(svg_prods)
        args = Dict{String, Any}(#
            "height" => canvas_height,
            "width" => canvas_width,
            "highlightColour" => [137/256, 49/256, 239/256, 0.9]
        )
        if i == prod_idx
            hl = atom_to_highlight(string(prodss[i]), only(string(res.i[3])), pstate[2])
            args["atoms"] = [hl]
        else
            if haskey(args, "atoms")
                delete!(args, "atoms")
            end
        end
        svg_prods[i] = get_svg(get_mol(prodss[i]), args)
    end

    # Initialize combined plot dimensions
    arrow_width = 0
    if isempty(prodss)
        arrow_width = 500
    else
        arrow_width = 300
    end
    width = [#
        canvas_width * length(svg_subs),
        100 * (length(svg_subs) - 1),
        arrow_width,
        canvas_width * length(svg_prods),
        100 * (length(svg_prods) - 1)
    ]
    width = sum(width)
    height_padding = width/3 - canvas_height
    height = canvas_height + height_padding

    # Initialize plot and draw reaction string
    Drawing(width, height, :svg)
    background("white")
    w = 0
    for i in 1:length(svg_subs)
        placeimage(readsvg(svg_subs[i]), Point(w,height_padding/2))
        fontsize(20)
        sethue("black")
        text(metnames[i], Point(w+250,height_padding/2+20), halign=:center)
        w += 350
        if length(svg_subs) > 1 && i != length(svg_subs)
            fontsize(50)
            sethue("black")
            w += 150
            text("+", Point(w,height_padding/2 + 350/2-25), valign=:middle)
            w += 30
        end
    end
    fontsize(50)
    sethue("black")
    w += 150
    arrow(#
        Point(w,height_padding/2+350/2-25),
        Point(w+350,height_padding/2+350/2-25);
        linewidth = 5.0,
        arrowheadlength = 25,
        arrowheadangle = pi/8,
        decoration = 0.5
    )
    w += 350
    for i in 1:length(svg_prods)
        placeimage(readsvg(svg_prods[i]), Point(w,height_padding/2))
        fontsize(20)
        sethue("black")
        text(metnames[i+length(svg_subs)], Point(w+250,height_padding/2+20), halign=:center)
        w += 350
        if length(svg_prods) > 1 && i != length(svg_prods)
            fontsize(50)
            sethue("black")
            w += 150
            text("+", Point(w,height_padding/2+350/2-25), valign=:middle)
            w += 30
        end
    end

    finish()
    s = svgstring()
end

