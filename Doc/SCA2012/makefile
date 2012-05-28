# variables:
#   comments    1 to show comments in the output
#   itype       dvi/pdf - type of in-time build
#
# Rtargets is timestamp file used to prevent R based images from rebuilding all
# the time

FILE:=sca2012
SUBFILES:=
BIBFILE:=sca2012


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

# Build type for in-time target
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


TEXFILES=$(patsubst %,%.tex, $(FILE))# $(patsubst %,%.tex, $(SUBFILES)) 
PDFFILES=$(patsubst %,%.pdf, $(FILE)) 
PSFILES=$(patsubst %,%.ps, $(FILE)) 

comma:= ,
empty:=
space:= $(empty) $(empty)
FILEPAT={$(subst $(space),$(comma),$(FILE))}


all: $(PDFFILES) #$(PSFILES)

##  Real targets  ##############################################################

# TODO: dvips removes .dvi file if I do: make clean && make foo.ps
%.ps: %.dvi
	HOME=`pwd` $(DVIPS) -o $@ $<

%.pdf: %.tex $(TEXFILES) $(BIBFILE).bib
	export HOME=`pwd`; \
	export TEXMFLOCAL=`pwd`/texmf:`kpsexpand '$$TEXMFLOCAL'`; \
	mkdir -p "$(TMPDIR)" && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' ; \
	BIBINPUTS=.:$(TMPDIR) $(BIBTEX) $(TMPDIR)/$* && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' && \
	$(PDFLATEX) '\def\GraphicsType{pdftex}$(TEXCMDS)\input $*' && \
	$(CP) $(TMPDIR)/$@ .

%.dvi: %.tex $(TEXFILES) $(BIBFILE).bib
	export HOME=`pwd`; \
	export TEXMFLOCAL=`pwd`/texmf:`kpsexpand '$$TEXMFLOCAL'`; \
	mkdir -p "$(TMPDIR)" && \
	$(LATEX) '$(TEXCMDS)\input $*' ; \
	BIBINPUTS=.:$(TMPDIR) $(BIBTEX) $(TMPDIR)/$* && \
	$(LATEX) '$(TEXCMDS)\input $*' && \
	$(LATEX) '$(TEXCMDS)\input $*' && \
	$(CP) $(TMPDIR)/$@ .

# TODO: pdf doesn't work? 
%.eps %.png: %.asy
	asy -f png $<
	asy -f eps $<

%.eps %.png: %.dia
	dia -t eps $<
	dia -t png $<

%.eps: %.png
	convert $< $@

##  Other targets  #############################################################

view: view-dvi

view-dvi:
	@HOME=`pwd` $(DISP) exec xdvi -q -watchfile 5 $(FILE).dvi

view-ps:
	#HOME=`pwd` $(DISP) gv -center -watch -widgetless $(FILE).ps
	@$(DISP) exec gv -center -watch -widgetless $(FILE).ps

rebuild: clean all

#$(TMPDIR)/$(FILE).{aux,log,ps,dvi,pdf,toc,bbl,blg} 
# TODO: fix this, doesn't work if TMPDIR = '.'
clean:
	$(RM) -r $(FILEPAT).{ps,dvi,pdf} missfont.log .texlive 
	[ "." != "$(TMPDIR)" ] && $(RM) -r "$(TMPDIR)" 

# following method is little bit consuming (lots of stat/read invocations) but
# it's most responsive so far
in-time:
	@echo "Building $(FILE).$(itype)  (modify with: make itype=<type> in-time)"; echo; \
	export WATCH="$(FILE).tex $(BIBFILE).bib"; \
	declare -i last cur; \
	last=0; \
	cur=0; \
	while true ; do \
		while [ $$last -eq $$cur ] ; do \
			cur=`stat -c'%Z' $$(eval echo $$WATCH) | cut -d. -f1 | sort -r | head -n1`; \
			read -t 1 x && last=cur-1; \
		done ; \
		make $(FILE).$(itype); \
		last=$$cur; \
	done

.PHONY: all clean rebuild view view-dvi view-ps in-time
