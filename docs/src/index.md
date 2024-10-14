# MarkovWeightedEFMs.jl

MarkovWeightedEFMs.jl is a package to decompose steady state metabolic
fluxes onto elementary flux mode (EFM) weights and atomic elementary flux
mode (AEFM) weights. EFM flux decomposition works only for closed-loop
networks of unimolecular reactions. AEFM flux decomposition works on
metabolic flux networks with known metabolite structures.

## Usage

See [Getting started](installation/installation.md) for package
installation and Python dependencies (RXNMapper) for AEFM enumeration and
flux decomposition. Once installed, the package is loaded in a `julia`
session by typing:

```julia
julia> using MarkovWeightedEFMs
```

See the tutorial sections for EFM or AEFM enumeration and flux decomposition.

## Citing MarkovWeightedEFMs.jl

Please cite the following papers if you use our method for (A)EFM enumeration
and flux decomposition.

Justin G. Chitpin and Theodore J. Perkins,
*Atomic elementary flux modes explain the steady state flow of metabolites in flux networks*.
biorXiv preprint doi: XX.XXXX/XXXX.XX.XX.XXXXXX

Justin G. Chitpin and Theodore J. Perkins,
*A Markov constraint to uniquely identify elementary flux mode weights in unimolecular metabolic networks*.
J Theor Biol. 2023 Nov 7;575:111632. doi: [10.1016/j.jtbi.2023.111632](https://doi.org/10.1016/j.jtbi.2023.111632)

## License

This software is released under the MIT license.

## Contact information

Justin G. Chitpin at jchit069@uottawa.ca for questions.
