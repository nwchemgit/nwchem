.SUFFIXES:    .tex .dvi .ps .log .aux

PSFLAGS	      = -D 300

TEX	      = user.tex

DVI           = user.dvi

PSFILE	      = user.ps


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
		latex \&psfont $*
		dvips $(PSFLAGS) -o $@ $<


####
user.dvi:	fullpage.sty dft.tex



