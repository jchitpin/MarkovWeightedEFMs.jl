# Export CHMC results
"""
    export_chmc(fname::String, res::CHMCStandardSummary)

Export CHMC results to text file `fname`.

`fname` is the filename to write the CHMC results.

`res` are the results from [`steady_state_efm_distribution`](@ref).
```
"""
function export_chmc(fname::String, res::CHMCStandardSummary)
  open(fname, "w") do f
    # EFMs
    write(f, "## EFMs:\n")
    for e in res.e
      write(f, "$e\n")
    end
    write(f, "\n")

    # EFM probabilities
    write(f, "## EFM probabilities:\n")
    for p in res.p
      write(f, "$p\n")
    end
    write(f, "\n")

    # EFM weights
    write(f, "## EFM weights:\n")
    for w in res.w
      write(f, "$w\n")
    end
    write(f, "\n")
  end
end

