#=
s = "[CH2:61]([CH:59]([CH:57]([CH:55]([CH:54]([C:2](=[O:1])[OH:3])[OH:53])[OH:56])[OH:58])[OH:60])[O:62][P:63](=[O:65])([OH:64])[OH:66].[OH:5][C@:6]1[C@:42]([OH:43])[C@@:9]([CH2:10][O:11][P:12]([OH:13])(=[O:14])[O:15][P:16]([OH:17])(=[O:18])[O:19][CH2:20][C@:21]2[O:22][C@:23]([N:24]3[CH:25]=[N:26][C:27]4=[C:28]([NH2:33])[N:29]=[CH:30][N:31]=[C:32]34)[C@@:34]([O:35][P:36]([OH:37])([OH:38])=[O:39])[C@@:40]2[OH:41])[O:8][C@@:7]1[N+:44]1=[CH:49][C:48]([C:50](=[O:51])[NH2:52])=[CH:47][CH:46]=[CH:45]1>>[O:1]=[C:2]=[O:3].[H+:4].[OH:5][C@:6]1[C@@:7]([N:44]2[CH:45]=[CH:46][CH2:47][C:48]([C:50](=[O:51])[NH2:52])=[CH:49]2)[O:8][C@@:9]([CH2:10][O:11][P:12]([OH:13])(=[O:14])[O:15][P:16]([OH:17])(=[O:18])[O:19][CH2:20][C@:21]2[O:22][C@:23]([N:24]3[CH:25]=[N:26][C:27]4=[C:28]([NH2:33])[N:29]=[CH:30][N:31]=[C:32]34)[C@@:34]([O:35][P:36]([OH:37])([OH:38])=[O:39])[C@@:40]2[OH:41])[C@:42]1[OH:43].[OH:53][CH2:54][C:55](=[O:56])[C:57]([OH:58])[C:59]([OH:60])[CH2:61][O:62][P:63]([OH:64])([OH:65])=[O:66]"
plot_mapped_reaction(s, "test.svg", view=false)
=#

include("helpers.jl")

"""
    function plot_mapped_reaction(#
        s::String,
        fname::String = "";
        view::Bool = false,
        canvas_width::Int64 = 3000,
        canvas_height::Int64 = 1000
)

Plot mapped reaction SMILES string `s` as an SVG and save to `fname` if specified.

`canvas_width` is the width of the SVG.

`canvas_height` is the height of the SVG.

`view=true` will plot the SVG assuming a plotting backend is specified. For
example, loading the `ElectronDisplay` package will plot the SVG to an
Electron window.
"""
function plot_mapped_reaction(#
    s::String,
    fname::String = "";
    view::Bool = false,
    canvas_width::Int64 = 3000,
    canvas_height::Int64 = 1000,
)

    # Reconstruct the original reaction strings
    #subs, _ = split(s, ">>")
    #max_subs = length(split(subs, "."))
    #rms = canonicalize_and_atom_map(s)
    #rss = join([#
        #join(first.(rms)[1:max_subs], "."),
        #">>",
        #join(first.(rms)[(max_subs+1):end], ".")
    #])

    # Plot reaction
    rxn = get_rxn(s)
    args = Dict{String, Any}(#
        "height" => canvas_height,
        "width" => canvas_width,
    )
    rxn_svg = get_rxn_svg(rxn, args)

    if isempty(fname)
        fname = :svg
    end
    Drawing(canvas_width, canvas_height, fname)
    origin()
    x = readsvg(rxn_svg)
    placeimage(x, centered = true)
    finish()
    if view == true
        preview()
    end
end

"""
    function plot_atomic_chmc(#
        res::CHMCAtomicSummary,
        S::Matrix{Int16},
        mets::Vector{String},
        rs::Vector{String};
        node_label_textsize::Real=15,
        edge_label_textsize::Real=12,
        arrow_shift::Real=0.85,
        x_pad::Real=0.75,
        y_pad::Real=0.75,
        tfactor::Real=0.15,
        tangents::Tuple{#
            Tuple{<:Real, <:Real},
            Tuple{<:Real, <:Real}
        } = ((1,0),(0,1)),
        show_all::Bool=false,
        width::Int64=620,
        height::Int64=310
)

Plot atomic cycle-history Markov chain.

`S` is the m by n stoichiometry matrix.

`mets` is the vector of metabolite names of length m.

`rs` is the vector of reaction SMILES strings of length n.

`node_label_textsize` is the text size of the node labels indexed from `T`.

`edge_label_textsize` is the text size of the edge labels taken from `T`.

`arrow_shift` is the percentage shift of the arrow head from src to dst.

`x_pad` is the left/right x coordinate padding of the plotting box.

`y_pad` is the up/down y coordinate padding of the plotting box.

`tfactor` scales the distance of the bezier control point relative to the
distance of the src and dst nodes.

`tangents` is the tangent of the src vertex and dst vertex.

`show_all=true` explicitly plots the upstream transition from all EFMs that
pass through the initial state/node `I`. By default, these arrows stemming
from the green nodes are omitted for visual clarity.

`width` and `height` specify the plotting window dimensions in pixel units.
"""
function plot_atomic_chmc(#
    res::CHMCAtomicSummary,
    S::Matrix{Int16},
    mets::Vector{String},
    rs::Vector{String};
    node_label_textsize::Real = 15,
    edge_label_textsize::Real = 12,
    arrow_shift::Real = 0.85,
    tfactor::Real = 0.15,
    x_pad::Real = 0.75,
    y_pad::Real = 0.75,
    tangents::Tuple{Tuple{<:Real, <:Real}, Tuple{<:Real, <:Real}} = ((1,0),(0,1)),
    show_all::Bool = false,
    width::Int64=620,
    height::Int64=310,
)

    # Initialize figure and set up 2-panel layout
    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98))
    main = fig[1, 1] = GridLayout()
    axtop = Axis(#
        main[1, 1],
        title = "Reaction mapping",
        width = width,
        height = height
    )
    axmain = Axis(#
        main[2, 1],
        title = "ACHMC rooted on $(res.i[3])$(res.i[2]) of $(mets[res.i[1]])",
        width = width,
        height = width,
    )
    colgap!(main, 0)
    rowgap!(main, 5)

    # New Markov chain transition probability matrix
    #A, E, t = MarkovWeightedEFMs.Plots.trie_adj_at_root(res.T, show_all, Int16(1))
    A, E, t = trie_adj_at_root(res.T, show_all, Int16(1))

    # Vector of node and edge colours
    #node_color = MarkovWeightedEFMs.Plots.color_nodes(t, res.T, 1)
    #edge_color = MarkovWeightedEFMs.Plots.color_edges(A, E)
    #node_labels = MarkovWeightedEFMs.Plots.label_nodes(t)
    #edge_labels = MarkovWeightedEFMs.Plots.label_edges(A, E)
    node_color = color_nodes(t, res.T, 1)
    edge_color = color_edges(A, E)
    node_labels = label_nodes(t)
    edge_labels = label_edges(A, E)
    for r in res.R
        if A[r.i, r.j] != 0
            A[r.i, r.j] = r.k
        end
    end
    #edge_attrs = MarkovWeightedEFMs.Plots.label_edge_fluxes(A)
    edge_attrs = label_edge_fluxes(A)
    edge_labels = edge_labels .* " (" .* string.(edge_attrs) .* ")"

    # Tangents specify curved edges for looped EFM transitions
    #tans = MarkovWeightedEFMs.Plots.tangent_edges(A, E, tangents)
    tans = tangent_edges(A, E, tangents)

    # Construct simple digraph
    G = SimpleDiGraph(A)

    # Construct atomic CHMC in main panel
    p = graphplot!(#
        axmain,
        G,
        nlabels = node_labels,
        elabels = edge_labels,
        labels_textsize = node_label_textsize,
        elabels_textsize = edge_label_textsize,
        arrow_show = true,
        arrow_shift = arrow_shift,
        node_color = node_color,
        tangents = tans,
        tfactor = repeat([tfactor], length(tans)),
        edge_color = edge_color
    )

    # Update main panel with custom tree layout
    #x1, x2, y1, y2 = MarkovWeightedEFMs.Plots.buchheim_layout(p, res.T, 1); # semicolon important!
    x1, x2, y1, y2 = buchheim_layout(p, res.T, 1); # semicolon important!

    # Padding boundary box
    x1 -= x_pad
    x2 += x_pad
    y1 += y_pad
    y2 -= y_pad

    # Remove figure clutter and correct aspect
    hidedecorations!(axmain)
    limits!(axmain, x1, x2, y2, y1)
    hidedecorations!(axtop)

    # Initialize default observable for axtop
    default_img = rotr90(svg_to_rgb_matrix(default_image()))
    o = Observable(default_img)

    # Set default image to axtop
    image!(axtop, o)
    w = size(default_img, 1)
    h = size(default_img, 2)
    limits!(axtop, 0, w, 0, h)

    # Click an edge to update axtop to that reaction
    function edge_click_action(idx, args...)
        # Reset any previously selected edge
        cidx = findall(==(:red), p.edge_color[])
        p.edge_color[][cidx] .= :gray

        # Show selected reaction
        if p.edge_color[][idx] != :red
            p.edge_color[][idx] = :red

            # Construct SVG of the reaction (needs refactoring with inputs)
            #@info(idx)
            rxn = parse(Int64, split(p.elabels[][idx], " ")[2][2:(end-1)])
            #@info(rxn)
            s = visualize_atomic_chmc_transition(#
                res, mets, S, rxn, A, node_labels, rs[rxn], idx
            )
            o[] = rotr90(svg_to_rgb_matrix(s))
            w = size(o.val, 1)
            h = size(o.val, 2)
            limits!(axtop, 0, w, 0, h)
            image!(axtop, o)
        else
            p.edge_color[][idx] = :gray
            o[] = default_img
            w = size(o.val, 1)
            h = size(o.val, 2)
            limits!(axtop, 0, w, 0, h)
            image!(axtop, o)
        end
        p.edge_color[] = p.edge_color[]
    end
    eclick = EdgeClickHandler(edge_click_action)
    register_interaction!(axmain, :eclick, eclick)

    # Don't forget to return the Makie plot
    fig
end

