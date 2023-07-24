# Compiling the JOSS Paper

To build the paper with [pandoc](https://pandoc.org/installing.html) (3.0.1):

```shell
pandoc --citeproc --bibliography=refs.bib -s paper.md -o gx_minimax.pdf --template latex.template
```

with a latex templated sourced from a Github [source]( https://github.com/Wandmalfarbe/pandoc-latex-template/blob/master/eisvogel.tex)
rather than the [JOSS template](https://github.com/sigsep/open-unmix-paper-joss/blob/master/latex.template).
