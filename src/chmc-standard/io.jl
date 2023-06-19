## Main functions
# Export CHMC results
function export_chmc(#
  fname::String,
  res::NamedTuple{#
    (:e, :p, :w),
    Tuple{#
      Vector{Vector{Int64}},
      Vector{Float64},
      Vector{Float64}
    }
  }
)
  open(fname, "w") do f
    # EFMs
    write(f, "## EFMs:\n")
    for e in res.e
      write(f, "$e\n")
    end
    write(f, "\n")

    # Atomic EFM probabilities
    write(f, "## EFM probabilities:\n")
    for p in res.p
      write(f, "$p\n")
    end
    write(f, "\n")

    # Atomic EFM weights
    write(f, "## EFM weights:\n")
    for w in res.w
      write(f, "$w\n")
    end
    write(f, "\n")
  end
end

# Import atomic CHMC results
function import_chmc(fname::String)
  # Load text file contents
  f = open(fname)
  lines = readlines(f)

  # Headers
  idx = findall(!isnothing, match.(r"## ", lines))

  # Extract data within each header and return CHMC results structure
  if length(idx) == 3 # standard CHMC
    efms =  eval.(Meta.parse.(lines[(idx[1]+1):(idx[2]-2)]))
    probs = parse.(Float64, lines[(idx[2]+1):(idx[3]-2)])
    weights = parse.(Float64, lines[(idx[3]+1):(length(lines)-1)])

    return (e=efms, p=probs, w=weights)
  elseif length(idx) == 6 # atomic CHMC
    init = eval(Meta.parse(lines[idx[1]+1]))
    aefms =  eval.(Meta.parse.(lines[(idx[2]+1):(idx[3]-2)]))
    probs = parse.(Float64, lines[(idx[3]+1):(idx[4]-2)])
    weights = parse.(Float64, lines[(idx[4]+1):(idx[5]-2)])
    d = split.(lines[(idx[5]+1):(idx[6]-2)], "\t")
    x = parse.(Int64, first.(d))
    y = eval.(Meta.parse.(last.(d)))
    dict = Dict(zip(x, tuple.(Int16.(first.(y)), Int16.(last.(y)))))
    s = parse(Int64, lines[idx[6]+1])

    return (e=aefms, p=probs, w=weights, d=dict, i=init, s=s)
  else
    @warn("Cannot import the specified text file.")
  end
end

