
#######################################################################
# --original SOURCE list
# SOURCE = prog.tex            disclaimer.tex         cvsintro.tex \
#          codingsty.tex       makefile_libs.tex      newmodules.tex \
#          common_errors.tex   rtdb.tex               ga.tex \
#          ffinp.tex           print_control.tex      geom.tex \
#          basis.tex           symmetry.tex           context.tex \
#          util.tex            util_ga.tex            integral_api.tex \
#          pstat.tex           integral_io.tex        int_api_spec.tex
# -- my SOURCE list (JMC -- 8/14/98)
 SOURCE = prog.tex            about.tex              generic.tex  \
          basis.tex           codingsty.tex          cvsintro.tex \
          design.tex          disclaimer.tex         ffinp.tex    \
          ga.tex              geom.tex               glossary.tex \
          instaport.tex       int_api_spec.tex       la.tex       \
          ma.tex              makefile_libs.tex      messages.tex \
          modules.tex         newdoc.tex             newmodules.tex \
          nwarch.tex          print_control.tex      prog_man_intro.tex \
          pstat.tex           rtdb.tex               utilities.tex \
          symmetry.tex        testing.tex            util.tex     \
          ChemIO.tex          ecce.tex

 SOURCE_NOT_USED = prog.old.tex

 TARG = prog

all:	$(SOURCE)
	pdflatex $(TARG).tex
	pdflatex $(TARG).tex
	pdflatex $(TARG).tex
postscript $(TARG).ps:	$(SOURCE)
	latex $(TARG).tex
	latex $(TARG).tex
	latex $(TARG).tex
	dvips -D300 -N0 -i -S200 -o $(TARG).ps $(TARG).dvi
	./progpsfix
odd:	$(SOURCE)
	latex $(TARG).tex
	latex $(TARG).tex
	latex $(TARG).tex
	dvips -D300 -N0 -i -S200 -A -o $(TARG).odd.ps $(TARG).dvi 
	dvips -D300 -N0 -i -S200 -B -o $(TARG).even.ps $(TARG).dvi 

quick:		$(SOURCE)
	latex $(TARG).tex
	dvips -D300 -N0 -i -S200 -o $(TARG).ps $(TARG).dvi

html:	
	@../htmlize prog
	@../update_www prog

view:	$(TARG).ps
	ghostview -bg white -fg black -bd black $(TARG).ps

pdf:    $(SOURCE) .IGNORE
	pdflatex -interaction=nonstopmode progpdf
	pdflatex -interaction=nonstopmode progpdf
	pdflatex -interaction=nonstopmode progpdf

.IGNORE:

clean:
	/bin/rm -f *.dvi *.ps *.toc *.aux *.log *~
	/bin/rm -rf prog
realclean:
	$(MAKE) clean 
	rm -f progpdf.pdf progpdf.out
