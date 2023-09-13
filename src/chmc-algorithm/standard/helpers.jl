# Close network if open
function close_network(S::Matrix{<:Int64}, v::Vector{<:Real})
  # Error-checking for a unimolecular network with steady state fluxes
  sanitize_flux(v)
  sanitize_stoich_flux(S, v)

  # Indices of source/sink columns
  col_source = findall(==(1), vec(sum(S, dims=1)))
  col_sink = findall(==(-1), vec(sum(S, dims=1)))

  # Close off stoichiometry matrix
  rvec = zeros(Int64, size(S,2))
  rvec[col_source] .= -1
  rvec[col_sink] .= 1
  return vcat(rvec', S)
end

