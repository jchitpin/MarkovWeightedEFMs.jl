#!/usr/bin/env bash

## Plot toy network as PNG
# Generate PNG
pdflatex -shell-escape toy-network-1.tex
pdflatex -shell-escape toy-network-1-prefix.tex

# Convert to PDF
#convert               \
   #-verbose           \
   #-density 500       \
   #-trim              \
    #toy-network-1.pdf \
   #-quality 100       \
   #-flatten           \
   #-sharpen 0x1.0     \
   #+repage      \
    #toy-network-1.png

convert toy-network-1.png -gravity East -chop 674x0 toy-network-1.png
convert toy-network-1.png -gravity West -chop 633x0 toy-network-1.png

rm toy-network-1.aux toy-network-1.log
rm toy-network-1-prefix.aux toy-network-1-prefix.log
mv toy-network-1.png ../assets/
mv toy-network-1-prefix.png ../assets/

julia toy-network-1-prefix.jl

convert prefix.png -gravity East -chop 250x0 prefix.png
convert prefix.png -gravity West -chop 250x0 prefix.png

mv prefix.png ../assets/

