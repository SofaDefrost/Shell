# variables:
#   draft       1 to make draft version
#   comments    1 to show comments in the output
#   itype       dvi/pdf - type of intime build
#
# Rtargets is timestamp file used to prevent R based images from rebuilding all
# the time
#
# If you get 'No rule to make target foo.pdf needed for all' error check the bib file name

FILES:=$(shell ls *.tex | sed 's/\.tex$$//')
BIBFILE:=bibdata


TMPDIR=tmp
LATEXARGS=-halt-on-error -interaction=nonstopmode -file-line-error -output-directory="$(TMPDIR)"

LATEX=/usr/bin/latex $(LATEXARGS)
PDFLATEX=/usr/bin/pdflatex $(LATEXARGS)
DVIPS=/usr/bin/dvips
DVIPDF=/usr/bin/dvipdf
BIBTEX=/usr/bin/bibtex

SHELL=/bin/bash
RM=rm -fv
CP=cp -v

# Build type for intime target
ifeq ($(itype),)
itype=dvi
endif

# Display for view-* targets
ifneq ($(disp),)
DISP=DISPLAY=$(disp)
endif

# Temporary directory
ifeq ($(TMPDIR),)
TMPDIR="tmp"
endif

# 'comments' variable
ifneq ($(comments),1)
comments=0
endif

TEXCMDS=\newcount\DoComments\DoComments=$(comments)


TEXFILES=$(patsubst %,%.tex, $(FILES)) 
PDFFILES=$(patsubst %,%.pdf, $(FILES)) 
PSFILES=$(patsubst %,%.ps, $(FILES)) 

comma:= ,
empty:=
space:= $(empty) $(empty)
ifeq ($(space),$(findstring $(space),$(FILES)))
    FILEPAT={$(subst $(space),$(comma),$(FILES))}
else
    FILEPAT=$(FILES)
endif


all: $(PDFFILES) #$(PSFILES)

##  Real targets  ##############################################################

# TODO: dvips removes .dvi file if I do: make clean && makes foo.ps
%.ps: %.dvi
	HOME=`pwd` $(DVIPS) -o $@ $<

# NOTE:  we need to remove old .toc file or hyperref breaks down terribly
# 			-- 2010-09-15 is it still valid?
%.pdf: %.tex $(BIBFILE).bib
	export HOME=`pwd`; \
	export TEXMFLOCAL=`pwd`/texmf:`kpsexpand '$$TEXMFLOCAL'`; \
	mkdir -p "$(TMPDIR)" && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' ; \
	BIBINPUTS=.:$(TMPDIR) $(BIBTEX) $(TMPDIR)/$* && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' && \
	$(CP) $(TMPDIR)/$@ .

%.dvi: %.tex $(BIBFILE).bib
	export HOME=`pwd`; \
	export TEXMFLOCAL=`pwd`/texmf:`kpsexpand '$$TEXMFLOCAL'`; \
	mkdir -p "$(TMPDIR)" && \
	$(LATEX) '$(TEXCMDS)\input $*' ; \
	BIBINPUTS=.:$(TMPDIR) $(BIBTEX) $(TMPDIR)/$* && \
	$(LATEX) '$(TEXCMDS)\input $*' && \
	$(LATEX) '$(TEXCMDS)\input $*' && \
	$(CP) $(TMPDIR)/$@ .

# Default Arguments

INKSCAPE_EPS=
INKSCAPE_PDF=--export-background=\#ffffff --export-dpi=800

## Dia
%.eps: %.dia
	dia -t eps $<

%.png: %.dia 
	dia -t png $<

## Images via ImageMagick
%.eps: %.png
	convert $< $@

%.eps: %.jpg
	convert $< $@

%.eps: %.pdf
	convert $< $@

%.pdf: %.png
	convert $< $@

## Inkscape
#TODO: inkscape does not return error values!
%.eps: %.svg
	inkscape $(INKSCAPE_EPS) --export-eps=$@ $<; [ -f "$@" ]

%.pdf: %.svg
	inkscape $(INKSCAPE_PDF) --export-pdf=$@ $<; [ -f "$@" ]

%.pdf_tex: %.svg
	inkscape $(INKSCAPE_PDF) --export-latex --export-pdf=$*.pdf $<; [ -f "$@" ]


##  Other targets  #############################################################

view: view-dvi

# TODO: add 'file' argument
view-dvi:
	@if [ -z "$(file)" ] ; then \
		echo "Don't know which file, use 'make file=<name> $@' to tell"; \
		false; \
	fi
	HOME=`pwd` $(DISP) exec xdvi -q -watchfile 5 $(file).dvi

# TODO: add 'file' argument
view-ps:
	@if [ -z "$(file)" ] ; then \
		echo "Don't know which file, use 'make file=<name> $@' to tell"; \
		false; \
	fi
	#HOME=`pwd` $(DISP) gv -center -watch -widgetless $(FILE).ps
	$(DISP) exec gv -center -watch -widgetless $(file).ps

rebuild: clean all

#$(TMPDIR)/$(FILE).{aux,log,ps,dvi,pdf,toc,bbl,blg} 
# TODO: fix this, doesn't work if TMPDIR = '.'
clean:
	$(RM) -r $(FILEPAT).{ps,dvi,pdf} missfont.log .texlive 
	[ "." != "$(TMPDIR)" ] && $(RM) -r "$(TMPDIR)" 

# following method is little bit consuming (lots of stat/read invocations) but
# it's most responsive so far
intime:
	@FPAT="$(file)"; [ -z "$$FPAT" ] && FPAT="$(FILEPAT)"; \
	echo "Building $$FPAT.$(itype)  (modify with: make itype=<type> intime)"; echo; \
	export WATCH="$${FPAT}.tex $(BIBFILE).bib"; \
	declare -i last cur; \
	last=0; \
	cur=0; \
	while true ; do \
		while [ $$last -eq $$cur ] ; do \
			cur=`stat -c'%Z' $$(eval echo $$WATCH) | cut -d. -f1 | sort -r | head -n1`; \
			read -t 1 x && last=cur-1; \
		done ; \
		make $${FPAT}.$(itype); \
		last=$$cur; \
	done

.PHONY: all clean rebuild view view-dvi view-ps intime
