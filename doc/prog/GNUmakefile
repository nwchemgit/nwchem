
#######################################################################
 SOURCE = prog.tex            disclaimer.tex         cvsintro.tex \
          codingsty.tex       makefile_libs.tex      newmodules.tex \
          common_errors.tex   rtdb.tex               ga.tex \
          ffinp.tex           print_control.tex      geom.tex \
          basis.tex           symmetry.tex           context.tex \
          util.tex            util_ga.tex            integral_api.tex \
          pstat.tex           integral_io.tex        int_api_spec.tex

 SOURCE_NOT_USED = prog.old.tex

 TARG = prog

all postscript $(TARG).ps:	$(SOURCE)
	latex $(TARG).tex
	latex $(TARG).tex
	latex $(TARG).tex
	dvips -D 300 -N0 -o $(TARG).ps $(TARG).dvi

odd:	$(SOURCE)
	latex $(TARG).tex
	latex $(TARG).tex
	latex $(TARG).tex
	dvips -D 300 -N0 -A -o $(TARG).odd.ps $(TARG).dvi 
	dvips -D 300 -N0 -B -o $(TARG).even.ps $(TARG).dvi 

quick:		$(SOURCE)
	latex $(TARG).tex
	dvips -D 300 -N0 -o $(TARG).ps $(TARG).dvi

html:	
	@../htmlize prog
	@../update_www prog

view:	$(TARG).ps
	ghostview -bg white -fg black -bd black $(TARG).ps

clean:
	/bin/rm -f *.dvi *.toc *.aux *.log 
realclean:	clean
	/bin/rm -f *~ *.ps
	/bin/rm -rf prog
