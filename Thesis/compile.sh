# Compiles paper properly, deletes temporary files and open the generated pdf (tested on Mac OS)

rm -f *.{blg,bbl,toc,out,log,aux,pdf,maf,brf,mtc,mtc0,mtc1,mtc2,mtc3,nlo,lot,lof,nls,ilg}

pdflatex Thesis
bibtex Thesis
makeindex Thesis.nlo -s nomencl.ist -o Thesis.nls
pdflatex Thesis
pdflatex Thesis

rm -f *.{blg,bbl,toc,out,log,aux,maf,brf,mtc,mtc0,mtc1,mtc2,mtc3,nlo,lot,lof,nls,ilg}

open Thesis.pdf