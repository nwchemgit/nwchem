.SUFFIXES:    .tex .dvi .ps .log .aux

PSFLAGS	      = -D 300

TEX	      = user.tex prog.tex md.tex
 
DVI           = user.dvi prog.dvi md.dvi

PSFILE	      = user.ps prog.ps md.ps


#
#
#


all:		$(PSFILE)



postscript:	$(PSFILE)



print:		postscript
		lp $(PSFILE) 


dvi:		$(DVI)



clean:		
		@rm -f $(DVI)  *.toc *.aux *.log *.bbl *.blg
realclean:	clean
		@rm -f $(PSFILE)


.tex.dvi:
		latex $*
		latex $*
		latex $*



.dvi.ps:
		dvips $(PSFLAGS) -o $@ $<


####
user.dvi:	fullpage.sty dft.tex
prog.dvi:	fullpage.sty integral_api.tex


