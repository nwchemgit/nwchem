
#
# Nasty little awk script to insert a CVS id comment line
# before the second executable statement
#

# $Id: fortranid.awk,v 1.3 1995-02-02 23:26:45 d3g681 Exp $

BEGIN {
   FIRST = 0;
   DONE  = 0;
}

DONE == 1  {print; next;}

/^[ \t][ \t][ \t][ \t][ \t][ \t]/ {
		if (FIRST) {
			printf("C$Id: fortranid.awk,v 1.3 1995-02-02 23:26:45 d3g681 Exp $\n");
			DONE = 1;
		} else {
			FIRST = 1;
		}
		print;
		next;
	}

		{print;}
