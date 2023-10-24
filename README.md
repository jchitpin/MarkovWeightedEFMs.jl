# MarkovWeightedEFMs.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jchitpin.github.io/MarkovWeightedEFMs.jl/dev/)

This package is used to decompose steady state metabolic fluxes onto elementary flux mode weights. Currently, this package only supports closed-loop networks of unimolecular reactions. Click the "docs" badge above to access package the documentation.

## Installation

To install this package, open a `julia` session and enter:

```julia
julia> ]
(@v1.6) pkg> add https://github.com/jchitpin/MarkovWeightedEFMs.jl.git
```

Alternatively, you can load the `Pkg` package and install by:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/jchitpin/MarkovWeightedEFMs.jl.git")
```


## Usage

Once installed, the package is loaded in a `julia` session by typing:

```julia
julia> using MarkovWeightedEFMs
```

Please read the docs for a guided tutorial and function descriptions.

## Citing MarkovWeightedEFMs.jl
Chitpin JG, Perkins TJ. *A Markov constraint to uniquely identify elementary flux mode weights in unimolecular metabolic networks*, J Theor Biol. 2023 Oct 5. doi: doi: 10.1016/j.jtbi.2023.111632. PMID: 37804942.

