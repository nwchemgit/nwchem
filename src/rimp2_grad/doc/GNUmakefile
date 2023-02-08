# $Id$

rimp2-grad.dvi dvi: rimp2-grad.tex 
	latex rimp2-grad
	latex rimp2-grad

rimp2-grad.ps ps postscript: rimp2-grad.dvi
	dvips -orimp2-grad.ps rimp2-grad

changebar.sty: changebar.ins changebar.dtx
	latex changebar.ins

subeqnarray.sty: subeqnarray.ins subeqnarray.dtx
	latex subeqnarray.ins

clean:
	/bin/rm -f *.dvi *.ps *.toc *.aux *.log *.cb *~

realclean: clean
	/bin/rm -f changebar.sty subeqnarray.sty subeqnarray.tex
