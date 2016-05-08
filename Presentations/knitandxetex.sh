#!/bin/bash

if [ -n $1 ]; then
  Rscript -e "library(methods);library(knitr); knit('$1.Rnw')"; xelatexmk.sh $1.tex
fi
