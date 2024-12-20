import Pkg
Pkg.add("Documenter")
using Documenter
using MarkovWeightedEFMs

push!(LOAD_PATH, "../src/")

makedocs(#
    sitename = "MarkovWeightedEFMs.jl",
    authors = "Justin G. Chitpin",
    format = Documenter.HTML(prettyurls = true), # comment out for release
    pages = [#
        "Home" => Any[#
            "index.md"
        ],
        "Installation" => Any[#
            "installation/installation.md",
            "installation/python-dependencies.md"
        ],
        "Tutorials" => Any[#
            "tutorials/chmc-standard-metabolic-networks.md",
            "tutorials/chmc-standard-ion-channels.md",
            "tutorials/chmc-atomic-glucose.md",
            "tutorials/chmc-atomic-glucose-under-the-hood.md",
            "tutorials/boilerplate-for-bigg-gems.md",
        ],
        "Library" => Any[#
            "library/chmc-standard.md",
            "library/chmc-atomic.md",
            "library/chmc-plots-standard.md",
            "library/chmc-plots-atomic.md",
        ]
    ]
)

deploydocs(#
    repo="github.com/jchitpin/MarkovWeightedEFMs.jl.git"
)

