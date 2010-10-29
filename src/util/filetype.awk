# attempt to determine the type of a file from its name

# $Id$

/\.fh$/			{print "Fortran-header"; next;}

/\.F$/ || /\.f$/	{print "Fortran"; next;}

/\.h$/			{print "C-header"; next;}

/\.c$/			{print "C"; next;}

/[Mm]akefile$/		{print "Makefile"; next;}

			{print "Unknown"; next;}
