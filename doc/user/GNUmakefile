
 SOURCE = ack.tex execution.tex mcscf.tex sample.tex symex.tex \
          basis.tex fourindex.tex md.tex  scf.tex  top-level.tex \
          dft.tex  geometry.tex mp2.tex  scfgrad.tex user.tex \
          disclaimer.tex intro.tex plnwv.tex selci.tex vib.tex \
          dplot.tex knownbas.tex property.tex stepper.tex nwargos.tex

all postscript user.ps:	$(SOURCE)
	latex user.tex
	latex user.tex
	latex user.tex
	dvips -D 300 -o user.ps user.dvi

html:	
	../htmlize user

realclean clean:
	/bin/rm -f *dvi *ps *toc *aux *log *~
	/bin/rm -rf user

