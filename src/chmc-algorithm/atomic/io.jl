# Export precomputed atom tracing dictionary
"""
    export_atom_tracing_dictionary(#
        fname::String,
        D::Dict{NTuple{4, Int64}, Tuple{Int64, Int64}}
    )

Export atom tracing dictionary `D` to text file `fname`.

`fname` is the filename to write the atom tracing dictionary results.

`D` is the dictionary from [`precompute_atom_tracing_dictionary`](@ref).
"""
function export_atom_tracing_dictionary(#
    fname::String,
    D::Dict{NTuple{4, Int64}, Tuple{Int64, Int64}}
)
    CSV.write(fname, D, header=false)
end

# Import precomputed atom tracing dictionary
"""
    import_atom_tracing_dictionary(fname::String)

Import atomic tracing dictionary from text file `fname`.

`fname` is the filename containing the atom tracing dictionary results.
"""
function import_atom_tracing_dictionary(fname::String)
    lines = CSV.File(fname, header=false)
    d = Dict{NTuple{4, Int64}, Tuple{Int64, Int64}}()
    for l in lines
        d[eval(Meta.parse(l[1]))] = eval(Meta.parse(l[2]))
    end
    return d
end

