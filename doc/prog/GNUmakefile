
 SOURCE = prog.tex

all postscript prog.ps:	$(SOURCE)
	latex prog.tex
	latex prog.tex
	latex prog.tex
	dvips -D 300 -N0 -o prog.ps prog.dvi

quick:		$(SOURCE)
	latex prog.tex
	dvips -D 300 -N0 -o prog.ps prog.dvi

html:	
	../htmlize prog

view:	prog.ps
	ghostview -bg white -fg black -bd black prog.ps

realclean clean:
	/bin/rm -f *dvi *ps *toc *aux *log *~
	/bin/rm -rf prog
