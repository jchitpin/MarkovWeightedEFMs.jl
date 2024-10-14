#!/usr/bin/env bash

## Tutorial: CHMC (metabolic networks)
# Plot toy network and CHMC as PDF and PNG
pdflatex -shell-escape toy-network-1.tex
pdflatex -shell-escape toy-network-1-chmc.tex

# Plot toy network CHMC via MarkovWeightedEFMs
julia toy-network-1-chmc-analysis.jl # generates toy-network-1-chmc-makie.png

# Crop PNGs and set background solid white
magick toy-network-1.png -gravity East -chop 672x0 toy-network-1.png
magick toy-network-1.png -gravity West -chop 632x0 -background white -alpha remove -alpha off toy-network-1.png
magick toy-network-1-chmc.png -background white -alpha remove -alpha off toy-network-1-chmc.png
magick toy-network-1-chmc-makie.png -gravity East -chop 250x0 toy-network-1-chmc-makie.png
magick toy-network-1-chmc-makie.png -gravity West -chop 250x0 -background white -alpha remove -alpha off toy-network-1-chmc-makie.png

# Remove logs
rm toy-network-1.aux toy-network-1.log
rm toy-network-1-chmc.aux toy-network-1-chmc.log

# Move PNGs to assets
mv toy-network-1.png ../assets/
mv toy-network-1-chmc.png ../assets/
mv toy-network-1-chmc-makie.png ../assets/

## Tutorial: CHMC (ion channels)
# Plot Markov state model as PDF and PNG
pdflatex -shell-escape ion-channel-mc.tex

# Plot toy network CHMC via MarkovWeightedEFMs
julia ion-channel-chmc-analysis.jl # generates ion-channel-chmc-makie.png

# Set background solid white
magick ion-channel-mc.png -background white -alpha remove -alpha off ion-channel-mc.png

# Remove logs
rm ion-channel-mc.aux ion-channel-mc.log

# Move PNGs to assets
mv ion-channel-mc.png ../assets/
mv ion-channel-chmc-makie.png ../assets/

## Tutorial ACHMC (one-carbon)
# Plot toy network as PDF and PNG
pdflatex -shell-escape toy-network-2-achmc.tex

# Set background solid white
magick toy-network-2-achmc.png -background white -alpha remove -alpha off toy-network-2-achmc.png

# Remove logs
rm toy-network-2-achmc.aux toy-network-2-achmc.log

# Move PNGs to assets
mv toy-network-2-achmc.png ../assets/

## Tutorial: ACHMC (glucose)
# Plot network as PDF and PNG
pdflatex -shell-escape toy-network-1-achmc.tex

# Set background solid white
magick toy-network-1-achmc.png -background white -alpha remove -alpha off toy-network-1-achmc.png

# Remove logs
rm toy-network-1-achmc.aux toy-network-1-achmc.log

# Move PNGs to assets
mv toy-network-1-achmc.png ../assets/

## Logo
pdflatex -shell-escape logo.tex

# Remove logs
rm logo.aux logo.log

# Convert to SVG
inkscape --export-type=svg logo.pdf

# Move logo to assets
mv logo.svg ../assets/

