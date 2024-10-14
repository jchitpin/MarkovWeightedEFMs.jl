# Getting started

To install this package, open a `julia` session and enter:

```julia
julia> ]
(@v1.10) pkg> add https://github.com/jchitpin/MarkovWeightedEFMs.jl.git
```

Alternatively, you can load `Pkg` and install by:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/jchitpin/MarkovWeightedEFMs.jl.git")
```

## Python dependencies for AEFM analysis

AEFM-specific analyses depend on the atom mapping program RXNMapper. This
package must be installed and built with PyCall.jl after installing
MarkovWeightedEFMs.jl. **Tested with python version 3.10**

This can be done by creating a python virtual environment, installing
RXNMapper, and setting the `PYTHON` environment variable to the python
executable in the virtual environment.

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

Note `PyCall.jl` will need to be rebuilt whenever you update your Julia version.
See [PyCall.jl documentation](https://github.com/JuliaPy/PyCall.jl) for
more options on setting up Python in Julia.

