# $Id: GNUmakefile,v 1.14 2000-02-04 21:44:57 edo Exp $

 SOURCE = ack.tex execution.tex mcscf.tex sample.tex symex.tex \
          basis.tex fourindex.tex scf.tex  top-level.tex \
          dft.tex gapss.tex geometry.tex mp2.tex  scfgrad.tex user.tex \
          disclaimer.tex intro.tex plnwv.tex selci.tex vib.tex \
          dplot.tex knownbas.tex property.tex stepper.tex nwargos.tex \
          ecp.tex functionality.tex ccsd.tex qmmm.tex

 SOURCE_NOT_USED = md.tex

 DVIPS_OPTS = -D300 -N0 -i -S200

all postscript user.ps:	$(SOURCE)
	rm -f user.[0-9][0-9][0-9] user.[0-9][0-9][0-9].ps
	latex user.tex
	latex user.tex
	latex user.tex
	dvips user ${DVIPS_OPTS} -o
	./userpsfix
odd:	$(SOURCE)
	rm -f user.[0-9][0-9][0-9] user.[0-9][0-9][0-9].ps usereven.ps userodd.ps
	latex user.tex
	latex user.tex
	latex user.tex
	dvips user -D300 -N0 -A -o userodd.ps 
	dvips user -D300 -N0 -B -o usereven.ps 

quick:		$(SOURCE)
	rm -f user.[0-9][0-9][0-9] user.[0-9][0-9][0-9].ps
	latex user.tex
	dvips user ${DVIPS_OPTS} -o
	./userpsfix
html:	
	@../htmlize user
	@../update_www user

html_only:
	@../htmlize user

view:	user.ps
	ghostview -bg white -fg black -bd black user.001.ps

pdf:	$(SOURCE) .IGNORE
	rm -f user.out user.aux
	pdflatex -interaction=nonstopmode user
	pdflatex -interaction=nonstopmode user
	pdflatex -interaction=nonstopmode user

.IGNORE:

clean:
	/bin/rm -f *.dvi *.ps *.toc *.aux *.log *~
	/bin/rm -rf user
realclean:
	$(MAKE) clean
