# MarkovWeightedEFMs.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jchitpin.github.io/MarkovWeightedEFMs.jl/dev/)

This package is used to decompose steady state metabolic fluxes onto
elementary flux mode (EFM) weights and atomic elementary flux mode (AEFM)
weights. EFM flux decomposition works only for closed-loop networks of
unimolecular reactions. AEFM flux decomposition works on metabolic flux
networks with known metabolite structures.

Click the "docs" badge above to access package the documentation.

## Installation

To install this package, open a `julia` session and enter:

```julia
julia> ]
(@v1.10) pkg> add https://github.com/jchitpin/MarkovWeightedEFMs.jl.git
```

Alternatively, you can load the `Pkg` package and install by:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/jchitpin/MarkovWeightedEFMs.jl.git")
```

### Python dependencies for AEFM analysis

AEFM-specific analyses depend on the atom mapping program RXNMapper. This
package must be installed and built with PyCall.jl after installing
MarkovWeightedEFMs.jl. **Tested with Python version 3.10**

```console
$ pip install virtualenv
$ virtualenv --python="/usr/bin/python3.10" "virtualenv" # name of virtual environment
$ source virtualenv/bin/activate
(virtualenv) $ pip install rxnmapper
(virtualenv) $ pip install rdkit
(virtualenv) $ pip install requests
(virtualenv) $ pip install tdqm
(virtualenv) $ pip install bs4
(virtualenv) $ pip install CTSgetPy
(virtualenv) $ julia
```
```julia
julia> using Pkg, PyCall
julia> ENV["PYTHON"] = joinpath(ENV["VIRTUAL_ENV"], "bin", "python")
julia> Pkg.build("PyCall")
```

## Usage

Once installed, the package is loaded in a `julia` session by typing:

```julia
julia> using MarkovWeightedEFMs
```

Please read the docs for guided tutorials and function descriptions.

## Citing MarkovWeightedEFMs.jl

Justin G. Chitpin and Theodore J. Perkins,
*Atomic elementary flux modes explain the steady state flow of metabolites in flux networks*.
biorXiv preprint: https://doi.org/XX.XXXX/XXXX.XX.XX.XXXXXX

Justin G. Chitpin and Theodore J. Perkins,
*A Markov constraint to uniquely identify elementary flux mode weights in unimolecular metabolic networks*.
J Theor Biol. 2023 Nov 7;575:111632.
doi: 10.1016/j.jtbi.2023.111632. PMID: 37804942.

## License

This software is released under the MIT license.
