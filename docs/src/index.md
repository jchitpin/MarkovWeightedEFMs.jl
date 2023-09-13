# MarkovWeightedEFMs.jl

MarkovWeightedEFMs.jl is a package for computing EFM probabilities and weights for steady state metabolic flux networks.

## Installation

To install this package, open a `julia` session and enter:

```julia
julia> ]
(@v1.6) pkg> add https://github.com/jchitpin/MarkovWeightedEFMs.jl.git
```

## Usage

Once installed, the package is loaded in a `julia` session by typing:

```julia
julia> using MarkovWeightedEFMs
```

Note: the one plotting function in this package relies on the GLMakie backend. Errors with `tree_plot()` are probably related to OpenGL. See <https://github.com/JuliaPlots/GLMakie.jl#troubleshooting-opengl> for troubleshooting tips.

## Citing MarkovWeightedEFMs.jl

Justin G Chitpin and Theodore J Perkins,
*A Markov constraint to uniquely identify elementary flux mode weights in unimolecular metabolic networks*, biorXiv preprint **biorXiv:2022.07.25.501464**, doi: https://doi.org/10.1101/2022.07.25.501464, 2022.

## License

This software is released under the MIT license.
