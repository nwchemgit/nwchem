.SUFFIXES:    .tex .dvi .ps .log .aux

PSFLAGS	      = -D 300

TEX	      = user.tex prog.tex

DVI           = user.dvi prog.ps

PSFILE	      = user.ps prog.ps


#
#
#


all:		$(PSFILE)



postscript:	$(PSFILE)



print:		postscript
		lp $(PSFILE) 


dvi:		$(DVI)



clean:		
		@rm -f $(DVI) $(PSFILE) *.toc *.aux *.log *.bbl *.blg



.tex.dvi:
		latex $*
		latex $*
		latex $*



.dvi.ps:
		dvips $(PSFLAGS) -o $@ $<


####
user.dvi:	fullpage.sty dft.tex
prog.dvi:	fullpage.sty
