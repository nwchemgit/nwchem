
#
# Nasty little awk script to insert a CVS id comment line
# before the second executable statement
#

BEGIN {
   FIRST = 0;
   DONE  = 0;
}

DONE == 1  {print; next;}

/^[ \t][ \t][ \t][ \t][ \t][ \t]/ {
		if (FIRST) {
			printf("C$id$\n");
			DONE = 1;
		} else {
			FIRST = 1;
		}
		print;
		next;
	}

		{print;}
