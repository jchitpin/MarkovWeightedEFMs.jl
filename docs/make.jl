import Pkg
Pkg.add("Documenter")
using Documenter
using MarkovWeightedEFMs

push!(LOAD_PATH, "../src/")

makedocs(#
  sitename = "MarkovWeightedEFMs.jl",
  authors = "Justin G. Chitpin",
  format = Documenter.HTML(prettyurls = false), # comment out for release
  pages = [#
    "Home" => Any[#
      "index.md"
    ],
    "Tutorials" => Any[#
      "tutorials/chmc-standard-metabolic-networks.md",
      "tutorials/chmc-standard-ion-channels.md",
      #"tutorials/chmc-atomic.md",
    ],
    "Library" => Any[#
      "library/chmc-standard.md",
      #"library/chmc-atomic.md",
      "library/chmc-plots.md",
    ]
  ]
)

deploydocs(#
    repo="github.com/jchitpin/MarkovWeightedEFMs.jl.git"
)

