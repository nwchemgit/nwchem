all:	user.ps author.ps
user.ps:	user.tex	fullpage.sty
	latex user
	latex user
	latex user
	dvips user -o 
uview:	user.ps
	ghostview user.ps
author.ps:	author.tex	fullpage.sty
	latex author
	dvips author -o 
aview:	author.ps
	ghostview author.ps
hardcopy:	user.ps	author.ps
	lpr -Pm1137 user.ps
	lpr -Pm1137 author.ps
both_view:	uview	aview
clean:
	/bin/rm -f *.dvi *.aux *.log *.ps 
