
 SOURCE = prog.tex

all postscript prog.ps:	$(SOURCE)
	latex prog.tex
	latex prog.tex
	latex prog.tex
	dvips -D 300 -o prog.ps prog.dvi

html:	
	../htmlize prog

realclean clean:
	/bin/rm -f *dvi *ps *toc *aux *log *~
	/bin/rm -rf prog
