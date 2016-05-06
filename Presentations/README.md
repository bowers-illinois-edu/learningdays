
To make the presentation, I do the following:

```
Rscript -e "library(methods);library(knitr); knit('presentation.Rnw',output='presentation.tex')"
latexmk -pdflatex='xelatex --shell-escape' -pdf presentation.tex
```

