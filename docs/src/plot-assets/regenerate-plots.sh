#!/usr/bin/env bash

## Plot toy network as PNG
# Generate PNG
pdflatex -shell-escape toy-network-1.tex
pdflatex -shell-escape toy-network-1-chmc.tex

# Convert to PDF
convert toy-network-1.png -gravity East -chop 674x0 toy-network-1.png
convert toy-network-1.png -gravity West -chop 633x0 toy-network-1.png

rm toy-network-1.aux toy-network-1.log
rm toy-network-1-chmc.aux toy-network-1-chmc.log
mv toy-network-1.png ../assets/
mv toy-network-1-chmc.png ../assets/

julia toy-network-1-chmc.jl

convert toy-network-1-chmc-tree_plot.png -gravity East -chop 250x0 toy-network-1-chmc-tree_plot.png
convert toy-network-1-chmc-tree_plot.png -gravity West -chop 250x0 toy-network-1-chmc-tree_plot.png

mv toy-network-1-chmc-tree_plot.png ../assets/

