# $Id: GNUmakefile,v 1.8 1997-02-24 19:42:05 d3e129 Exp $

 SOURCE = ack.tex execution.tex mcscf.tex sample.tex symex.tex \
          basis.tex fourindex.tex scf.tex  top-level.tex \
          dft.tex  geometry.tex mp2.tex  scfgrad.tex user.tex \
          disclaimer.tex intro.tex plnwv.tex selci.tex vib.tex \
          dplot.tex knownbas.tex property.tex stepper.tex nwargos.tex \
          ecp.tex functionality.tex ccsd.tex qmmm.tex

 SOURCE_NOT_USED = md.tex

all postscript user.ps:	$(SOURCE)
	latex user.tex
	latex user.tex
	latex user.tex
	dvips -D 300 -N0 -o user.ps user.dvi
odd:	$(SOURCE)
	latex user.tex
	latex user.tex
	latex user.tex
	dvips -D 300 -N0 -A -o userodd.ps user.dvi 
	dvips -D 300 -N0 -B -o usereven.ps user.dvi 

quick:		$(SOURCE)
	latex user.tex
	dvips -D 300 -N0 -o user.ps user.dvi
html:	
	../htmlize user

view:	user.ps
	ghostview -bg white -fg black -bd black user.ps

realclean clean:
	/bin/rm -f *.dvi *.ps *.toc *.aux *.log *~
	/bin/rm -rf user

