user.ps:	user.tex
	latex user
	dvips user -o 
view:	user.ps
	ghostview user.ps
hardcopy:	user.ps
	lpr -Pm1137 user.ps
clean:
	/bin/rm -f *.dvi *.aux *.log *.ps 
