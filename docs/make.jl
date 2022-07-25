import Pkg
Pkg.add("Documenter")
using Documenter
using MarkovWeightedEFMs

push!(LOAD_PATH, "../src/")

makedocs(#
  sitename = "MarkovWeightedEFMs.jl",
  #format = Documenter.HTML(prettyurls = false),
  pages = [#
    "Home" => Any[#
      "index.md"
    ],
    "Manual" => Any[#
      "tutorials/efm-estimation.md",
      "functions/cycle-history-markov-chain.md"
    ]
  ]
)

deploydocs(#
    repo="github.com/jchitpin/MarkovWeightedEFMs.jl.git"
)
