"""
Close an open-loop, unimolecular network by introducing a pseudo-metabolite.
"""
function close_network(S::Matrix{<:Integer}, v::Vector{<:Real})
    # Error-checking for a unimolecular network with steady state fluxes
    sanitize_flux(v)
    sanitize_stoich_flux(S, v)

    # Indices of source/sink columns
    j_src = findall(==(+1), vec(sum(S, dims = 1)))
    j_snk = findall(==(-1), vec(sum(S, dims = 1)))

    # Close off stoichiometry matrix
    rvec = zeros(Int64, size(S, 2))
    rvec[j_src] .= -1
    rvec[j_snk] .= 1

    return vcat(rvec', S)
end

