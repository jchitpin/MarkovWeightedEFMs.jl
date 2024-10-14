# Structs
"""
Stores results of a standard CHMC.

`e` are the EFM arrays of metabolite index sequences.
`p` are the EFM probabilities normalized to unity.
`w` are the EFM weights.
"""
struct CHMCStandardSummary
    e::Vector{Vector{Int64}} # EFM metabolite index sequences
    p::Vector{Float64}
    w::Vector{Float64}
end

"""
Stores properties of the metabolic flux network inputs, including possible
problems that must be corrected to successfully construct the atomic CHMC.
The recorded problems are used to automatically corrected via
[`find_atomic_chmc_input_errors(S::Matrix{<:Real}, v::Vector{<:Real})`](@ref).
The properties may be printed via [`print(res::CHMCAtomicErrorSummary)`](@ref).
"""
struct CHMCAtomicErrorSummary
    absolute_flux_error::Float64
    reactions_duplicated::Vector{Int64}
    reactions_with_zero_flux::Vector{Int64}
    reactions_with_negative_flux::Vector{Int64}
    reactions_with_non_integer_stoichiometries::Vector{Int64}
    reactions_unimolecular_source_stoichiometry_one::Vector{Int64}
    reactions_unimolecular_source_stoichiometry_not_one::Vector{Int64}
    reactions_multimolecular_source_stoichiometry_one::Vector{Int64}
    reactions_multimolecular_source_stoichiometry_not_one::Vector{Int64}
    reactions_unimolecular_sink_stoichiometry_one::Vector{Int64}
    reactions_unimolecular_sink_stoichiometry_not_one::Vector{Int64}
    reactions_multimolecular_sink_stoichiometry_one::Vector{Int64}
    reactions_multimolecular_sink_stoichiometry_not_one::Vector{Int64}
    reactions_empty::Vector{Int64}
    unused_metabolites::Vector{Int64}
end

"""
Stores results of an atomic CHMC.

`i` is the initial starting state (source metabolite index, atom position, chemical symbol).
`e` are the atomic EFM sequences including the pairs of their simple cycle closures.
`p` are the atomic EFM probabilities normalized to unity.
`w` are the atomic EFM weights.
`dmc` is the dictionary of MC states.
`dchmc` is the dictionary of CHMC states.
`T` is the atomic CHMC transition probability matrix.
`R` are the entries in `T` that transform the atom in state `i` to state `j` by indices `k`.
"""
struct CHMCAtomicSummary
    i::Tuple{Int64,Int64,Symbol}
    e::Vector{#
        NamedTuple{#
            (:EFM,:Closures),
            Tuple{Vector{Int64},Vector{Tuple{Int64,Int64}}}
        }
    }
    p::Union{Nothing,Vector{Float64}}
    w::Union{Nothing,Vector{Float64}}
    dmc::Dict{Int64,Tuple{Int16,Int16}}
    dchmc::Dict{#
        Vector{Int16},
        NamedTuple{(:id,:children),Tuple{Int64,Vector{Int16}}}
    }
    T::SparseMatrixCSC{Float64,Int64}
    R::Vector{NamedTuple{(:i,:j,:k),Tuple{Int64,Int64,Int16}}}
end

"""
    findall_int16(f::Function, A)

Int16 version of base Julia `findall` returning matching indices as Int16.
"""
function findall_int16(f::Function, A)
    T = Int16
    gen = (first(p) for p in pairs(A) if f(last(p)))
    isconcretetype(T) ? collect(T, gen) : collect(gen)
end

"""
    check_open_closed(S::Matrix{<:Integer})

Returns "closed" or "open" if input stoichiometry matrix is closed or open loop.
"""
function check_open_closed(S::Matrix{<:Integer})
    j_src = findall(==(+1), vec(sum(S, dims = 1)))
    j_snk = findall(==(-1), vec(sum(S, dims = 1)))
    if isempty(j_src) && isempty(j_snk)
        return "closed"
    else
        j_sre = findall(==(+1), vec(sum(S, dims = 1)))
        j_snk = findall(==(-1), vec(sum(S, dims = 1)))
        bool = all([!isempty(j_src), !isempty(j_snk)])
        @assert(#
            bool,
            join([#
                "The network must be open-loop and contain at least one ",
                "source and sink reaction. An open network must contain at ",
                "least one source and sink reaction."
            ])
        )
        return "open"
    end
end

"""
    isrotation(efm::Vector{Int64}, scs::Vector{Vector{Int64}})

Identify all simple cycles in `scs` that are identical to or rotation of the
simple cycle `efm`. This function assumes the simple cycles in `scs` and `efm`
are broken such that the efm[1] != efm[end].
"""
function isrotation(efm::Vector{Int64}, scs::Vector{Vector{Int64}})
    idx = Vector{Int64}()
    ee = [efm; efm]
    for k in eachindex(scs)
        l = k
        if length(efm) == length(scs[k])
            i = findfirst(==(scs[k][1]), efm)
            if !isnothing(i)
                for j in 0:(length(scs[k]) - 1)
                    if ee[i+j] != scs[k][j + 1]
                        l = 0
                        break
                    end
                end
                if l != 0
                    push!(idx, k)
                end
            end
        end
    end

    return idx
end

function prefixes_mc_to_chmc(#
    d::Dict{Vector{Int16},NamedTuple{(:id,:children),Tuple{Int64,Vector{Int16}}}}
)
    pfx_chmc = Vector{Vector{Int64}}(undef, length(d))
    c = 0
    for i in eachindex(d)
        temp = Vector{Int64}()
        for j in eachindex(i)
            push!(temp, d[i[1:j]].id)
        end
        pfx_chmc[c+=1] = temp
    end
    return pfx_chmc
end

function parallel_prefixes_chmc_to_simple_cycles(#
    pfx_chmc::Vector{Vector{Int64}},
    T′::Union{Matrix{Float64}, SparseMatrixCSC};
    verbose::Bool = false
)
    if verbose == true
        p = Progress(length(pfx_chmc); dt = 1.0, desc = "")
    end

    sc_chmc = Vector{Vector{Vector{Int64}}}(undef, length(pfx_chmc))
    Threads.@threads for i in eachindex(pfx_chmc)
        upstream = findall(>(0), @view T′[pfx_chmc[i][end], :])
        tmp = Vector{Vector{Int64}}()
        for j in eachindex(upstream)
            idx = findfirst(pfx_chmc[i] .== upstream[j])
            if ~isnothing(idx)
                push!(tmp, [pfx_chmc[i][end]; pfx_chmc[i][idx:end]])
            end
        end
        sc_chmc[i] = tmp
        if verbose == true
            next!(p)
        end
    end
    if verbose == true
        finish!(p)
    end

    return reduce(vcat, sc_chmc)
end

function parallel_simple_cycles_chmc_to_mc(#
    sc_chmc::Vector{Vector{Int64}},
    d::Dict{Vector{Int16},NamedTuple{(:id,:children),Tuple{Int64,Vector{Int16}}}};
    verbose::Bool = false
)
    if verbose == true
        p = Progress(length(sc_chmc); dt = 1.0, desc = "")
    end

    e = Dict(d[k].id => k[end] for k in keys(d))
    sc_mc = Vector{Vector{Int64}}(undef, length(sc_chmc))
    Threads.@threads for i in eachindex(sc_chmc)
        temp = Vector{Int64}()
        for j in eachindex(sc_chmc[i])
            push!(temp, e[sc_chmc[i][j]])
        end
        sc_mc[i] = temp
        if verbose == true
            next!(p)
        end
    end
    if verbose == true
        finish!(p)
    end

    return sc_mc
end

function parallel_simple_cycles_to_efms(#
    sc_mc::Vector{Vector{Int64}},
    sc_chmc::Vector{Vector{Int64}};
    verbose::Bool = true
)
    if verbose == true
        p = Progress(Threads.nthreads(); dt = 1.0, desc = "")
    end

    # Set maximum number of chunks
    len = Threads.nthreads()
    if length(sc_mc) < Threads.nthreads()
        len = length(sc_mc)
    end

    x = Vector{Vector{Vector{Int64}}}(undef, len)
    y = Vector{Vector{Vector{Tuple{Int64,Int64}}}}(undef, len)
    sc_o2 = [i[2:end] for i in sc_mc]
    
    ids = collect(index_chunks(x; n=len))
    Threads.@threads for i in 1:len
        x[i], y[i] = group_simple_cycles(sc_o2[ids[i]], sc_mc[ids[i]], sc_chmc[ids[i]])
        if verbose == true
            next!(p)
        end
    end
    if verbose == true
        finish!(p)
    end
    x = reduce(vcat, x)
    y = reduce(vcat, y)

    res = NamedTuple{#
        (:EFM,:Closures),Tuple{Vector{Int64},Vector{Tuple{Int64,Int64}}}
    }[]
    x2 = [i[2:end] for i in x]
    while !isempty(x)
        ids = sort(isrotation(x2[1], x2))
        push!(res, (EFM = x[1], Closures = reduce(vcat, y[ids])))
        splice!(x2, ids)
        splice!(x, ids)
        splice!(y, ids)
    end

    return res
end

function group_simple_cycles(sc_o2, sc_mc, sc_chmc)
    g(x) = Tuple(x[1:2])
    x = Vector{Vector{Int64}}()
    y = Vector{Vector{Tuple{Int64,Int64}}}()
    while !isempty(sc_mc)
        ids = sort(isrotation(sc_o2[1], sc_o2))
        push!(x, sc_mc[1])
        push!(y, g.(sc_chmc[ids]))
        splice!(sc_mc, ids)
        splice!(sc_o2, ids)
        splice!(sc_chmc, ids)
    end

    return x, y
end

"""
    trie(#
        T::Union{#
            Matrix{<:Real},
            SparseMatrixCSC{Int16, Int64},
            SparseMatrixCSC{Float64, Int64}
        },
        I::Int64 = 1
    )

Constructs the CHMC dictionary from MC matrix and initial starting state.

`T` is the transition probability matrix.
`I` is the initial state to root the CHMC.
"""
function trie(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    I::Int16 = Int16(1)
)
    # First tumble down the tree for key/value pairs
    p, queue = trie_first_pass(T, I)

    # Traverse remaining prefixes in the queue
    p2 = Vector{Vector{Pair{Vector{Int16}, Vector{Int16}}}}(undef, length(queue))
    Threads.@threads for i in eachindex(queue)
        p2[i] = trie_standard(T, queue[i])
    end

    # Concatenate pairs and reshape to dictionary with CHMC state indices
    if !isempty(p2)
        p = [p; reduce(vcat, p2)]
    end

    return Dict(zip(first.(p), NamedTuple{(:id, :children)}.(tuple.(eachindex(p), last.(p)))))
end

function trie_first_pass(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    I::Int16 = Int16(1)
)
    p = Vector{Pair{Vector{Int16}, Vector{Int16}}}()
    pfxs = Vector{Vector{Int16}}()

    function traverse_trie(pfx, p, pfxs)
        ds = filter!(x -> x ∉ pfx, findall_int16(>(0), @view T[pfx[end], :]))
        push!(p, pfx => ds)
        if !isempty(ds)
            for d in ds[2:end]
                push!(pfxs, [pfx; d])
            end
            traverse_trie([pfx; ds[1]], p, pfxs)
        end
    end

    # Construct dictionary of prefixes
    traverse_trie(Int16[I], p, pfxs)

    return p, pfxs
end

function trie_standard(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    I::Vector{Int16}
)
    p = Vector{Pair{Vector{Int16}, Vector{Int16}}}()

    function traverse_trie(#
        pfx::Vector{Int16},
        T::Union{#
            Matrix{<:Real},
            SparseMatrixCSC{Int16, Int64},
            SparseMatrixCSC{Float64, Int64}
        },
        p::Vector{Pair{Vector{Int16}, Vector{Int16}}}
    )
        ds = filter!(x -> x ∉ pfx, findall_int16(>(0), @view T[pfx[end], :]))
        push!(p, pfx => ds)
        for d in ds
            traverse_trie([pfx; d], T, p)
        end
    end

    # Construct dictionary of prefixes
    traverse_trie(I, T, p)

    return p
end

"""
    trie_matrix(#
        T::Union{#
            Matrix{<:Real},
            SparseMatrixCSC{Int16, Int64},
            SparseMatrixCSC{Float64, Int64}
        },
        I::Int64 = 1,
        verbose::Bool=false
    )

Constructs the CHMC matrix from CHMC dictionary and returns both.

`T` is the transition probability matrix.
`I` is the initial state to root the CHMC.
"""
function trie_matrix(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    I::Int64=1;
    verbose::Bool=false
)
    # Construct trie to get number of nodes
    if verbose == true
        @info("       Constructing CHMC trie.")
    end
    d = trie(T, Int16(I))

    # First tumble down the tree to fill queue of prefixes
    if verbose == true
        @info("       Converting CHMC trie to CHMC probability matrix.")
    end
    rs, cs, vs, queue = trie_matrix_first_pass(T, d, I)

    # Traverse remaining prefixes in the queue
    rs2 = Vector{Vector{Int64}}(undef, length(queue))
    cs2 = Vector{Vector{Int64}}(undef, length(queue))
    vs2 = Vector{Vector{Float64}}(undef, length(queue))
    Threads.@threads for i in eachindex(queue)
        rs2[i], cs2[i], vs2[i] = trie_matrix_standard(T, d, queue[i])
    end

    # Initialize and return sparse or dense matrix
    if !isempty(rs2) && !isempty(cs2) && !isempty(vs2)
        rs = [rs; reduce(vcat, rs2)]
        cs = [cs; reduce(vcat, cs2)]
        vs = [vs; reduce(vcat, vs2)]
    end
    A = sparse(rs, cs, vs, length(d), length(d))
    if !issparse(T)
        return Matrix(A), d
    end

    return A, d
end

function trie_matrix_first_pass(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    d::Dict{#
        Vector{Int16},
        NamedTuple{(:id, :children), Tuple{Int64, Vector{Int16}}}
    },
    I::Int64 = 1
)
    # Construct (sparse) CHMC transition matrix with COO format
    rs = Vector{Int64}()
    cs = Vector{Int64}()
    vs = Vector{Float64}()

    # Queue
    pfxs = Vector{Vector{Int16}}()

    function traverse_trie(#
        pfx::Vector{Int16},
        d::Dict{#
            Vector{Int16},
            NamedTuple{(:id, :children), Tuple{Int64, Vector{Int16}}}
        },
        T::Union{Matrix{<:Real}, SparseMatrixCSC{Float64, Int64}},
        rs::Vector{Int64},
        cs::Vector{Int64},
        vs::Vector{Float64},
        pfxs::Vector{Vector{Int16}}
    )
        # Fill matrix with non-children (upstream) transitions
        parents = findall(>(0), @views T[pfx[end], :])
        filter!(x -> x∉ d[pfx].children, parents)
        for j in parents
            push!(rs, d[pfx].id)
            push!(cs, d[pfx[1:findfirst(x -> x == j, pfx)]].id)
            push!(vs, T[pfx[end], j])
        end

        # Fill matrix with children transitions
        for i in eachindex(d[pfx].children)
            push!(rs, d[pfx].id)
            push!(cs, d[[pfx; d[pfx].children[i]]].id)
            push!(vs, T[pfx[end], d[pfx].children[i]])
            if i > 1
                push!(pfxs, [pfx; d[pfx].children[i]])
            end
        end
        if !isempty(d[pfx].children)
            traverse_trie([pfx; d[pfx].children[1]], d, T, rs, cs, vs, pfxs)
        end
    end
    traverse_trie(Int16.([I]), d, T, rs, cs, vs, pfxs)

    return rs, cs, vs, pfxs
end

function trie_matrix_standard(#
    T::Union{#
        Matrix{<:Real},
        SparseMatrixCSC{Int16, Int64},
        SparseMatrixCSC{Float64, Int64}
    },
    d::Dict{#
        Vector{Int16},
        NamedTuple{(:id, :children), Tuple{Int64, Vector{Int16}}}
    },
    pfx::Vector{Int16}
)
    # Construct (sparse) CHMC transition matrix with COO format
    rs = Vector{Int64}()
    cs = Vector{Int64}()
    vs = Vector{Float64}()

    function traverse_trie(#
        pfx::Vector{Int16},
        d::Dict{#
            Vector{Int16},
            NamedTuple{(:id, :children), Tuple{Int64, Vector{Int16}}}
        },
        T::Union{Matrix{<:Real}, SparseMatrixCSC{Float64, Int64}},
        rs::Vector{Int64},
        cs::Vector{Int64},
        vs::Vector{Float64}
    )
        # Fill matrix with non-children (upstream) transitions
        parents = findall(>(0), @views T[pfx[end], :])
        filter!(x -> x∉ d[pfx].children, parents)
        for j in parents
            push!(rs, d[pfx].id)
            push!(cs, d[pfx[1:findfirst(x -> x == j, pfx)]].id)
            push!(vs, T[pfx[end], j])
        end

        # Fill matrix with children transitions
        for j in d[pfx].children
            push!(rs, d[pfx].id)
            push!(cs, d[[pfx; j]].id)
            push!(vs, T[pfx[end], j])
            traverse_trie([pfx; j], d, T, rs, cs, vs)
        end
    end
    traverse_trie(pfx, d, T, rs, cs, vs)

    return rs, cs, vs
end

"""
    enumerate_efms(#
        T′::Union{Matrix{Float64}, SparseMatrixCSC},
        d::Dict{#
            Vector{Int16},
            NamedTuple{(:id, :children), Tuple{String, Vector{Int16}}}
        };
        verbose::Bool = false
    )

Enumerates the EFMs from the CHMC.

`T′` is the CHMC transition probability matrix.
`d` is the CHMC dictionary.
"""
function enumerate_efms(#
    T′::Union{Matrix{Float64}, SparseMatrixCSC},
    d::Dict{#
        Vector{Int16},
        NamedTuple{(:id, :children), Tuple{Int64, Vector{Int16}}}
    };
    verbose::Bool = false
)
    # Prefixes with CHMC states
    if verbose == true
        @info("       Converting MC trie prefixes to CHMC prefixes.")
    end
    pfx_chmc = prefixes_mc_to_chmc(d)

    # For each CHMC prefix, check if the last state transitions back to a CHMC
    # state already contained in the prefix. (Transitioning "up" the prefix tree)
    # These simple cycles represent the EFMs
    if verbose == true
        @info("       Extracting simple cycles from all CHMC prefixes.")
    end
    sc_chmc = parallel_prefixes_chmc_to_simple_cycles(pfx_chmc, T′; verbose)

    # Convert simple cycle history states back to regular states
    if verbose == true
        @info("       Converting CHMC simple cycles to MC simple cycles.")
    end
    sc_mc = parallel_simple_cycles_chmc_to_mc(sc_chmc, d; verbose)

    # Aggregate CHMC simple cycles for each EFM
    if verbose == true
        @info("       Aggregating simple cycles across EFMs.")
    end
    return parallel_simple_cycles_to_efms(sc_mc, sc_chmc; verbose)
end

"""
    solve_efm_probabilities(#
        T′::Union{Matrix{<:Real}, SparseMatrixCSC},
        ϕ::Vector{#
            NamedTuple{#
                (:EFM, :Closures),
                Tuple{Vector{Int64},Vector{Tuple{Int64,Int64}}}
            }
        },
        solver::Symbol
    )

Solves the EFM probabilities from CHMC transition probability matrix `T′` and
array of EFMs and their simple cycles `ϕ`.

`solver` is the type used for eigenvector calculations. Default is `nothing`
and `LinearSolve` will pick the best solver. Consult `LinearSolve.jl` for 
specifying other solvers.
"""
function solve_efm_probabilities(#
    T′::Union{Matrix{<:Real}, SparseMatrixCSC},
    ϕ::Vector{#
        NamedTuple{#
            (:EFM, :Closures),
            Tuple{Vector{Int64},Vector{Tuple{Int64,Int64}}}
        }
    },
    solver = nothing
)
    π = linearsolve_eigenvector(T′, solver)
    @assert(all(.!isnan.(π)), "Eigenvector cannot contain NaNs.")
    @assert(all(.!isnan.(π)), "All eigenvector coefficients must be positive.")
    while true
        p = solve_probabilities(T′, ϕ, π)
        if all(.!isnan.(p)) && all(p .>= 0)
            return p
        end
    end
    #@assert(all(.!isnan.(p)), "Steady state CHMC probabilities cannot be NaN.")
end

function linearsolve_eigenvector(T::Union{Matrix{<:Real}, SparseMatrixCSC}, solver)
    A = (LinearAlgebra.I - T'[2:end, 2:end])
    b = Vector(T'[2:end, 1])
    lin_prob = LinearProblem(A, b)
    while true
        sol = solve(lin_prob, solver)
        if !any(isnan.(sol.u)) && all(sol.u .>= 0)
            π = [1; sol.u]
            return π / sum(π)
        end
    end
end

function solve_probabilities(#
    T′::Union{Matrix{<:Real}, SparseMatrixCSC},
    ϕ::Vector{#
        NamedTuple{#
            (:EFM, :Closures),
            Tuple{Vector{Int64},Vector{Tuple{Int64,Int64}}}
        }
    },
    π::Vector{Float64}
)
    p = Vector{Float64}(undef, length(ϕ))
    for i in eachindex(ϕ)
        for j in ϕ[i].Closures
            p[i] += π[j[1]] * T′[first(j), last(j)]
        end
    end

    return p
end

