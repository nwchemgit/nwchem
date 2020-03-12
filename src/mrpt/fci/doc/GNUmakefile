
SOURCE   = matrixelements.tex \
           mspingraph.eps  orbitalgraph.eps spingraph.eps \
           orbitalarcs.eps spinarcs.eps

matrixelements.pdf:	$(SOURCE)
	pdflatex matrixelements.tex
	pdflatex matrixelements.tex
	pdflatex matrixelements.tex

clean:
	/bin/rm -f *toc *aux *log *dvi *pdf *~ core *bbl *blg 
