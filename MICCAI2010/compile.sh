# Compiles paper properly, deletes temporary files and open the generated pdf

rm paper.aux paper.bbl paper.blg paper.log

pdflatex paper.tex
bibtex paper.aux
pdflatex paper.tex
pdflatex paper.tex

rm paper.aux paper.bbl paper.blg paper.log

open paper.pdf