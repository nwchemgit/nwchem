# attempt to determine the type of a file from its name

# $Id: filetype.awk,v 1.2 1995-02-02 18:09:23 d3g681 Exp $

/\.fh$/			{print "Fortran-header"; next;}

/\.F$/ || /\.f$/	{print "Fortran"; next;}

/\.h$/			{print "C-header"; next;}

/\.c$/			{print "C"; next;}

/[Mm]akefile$/		{print "Makefile"; next;}

			{print "Unknown"; next;}
