#
# This script is used by the makefile in building util_version.f.
# It expects on standard input lines containing
#
# a CVS ID    -> produce a write statement with the ID
# module name -> produce a write statement for the module name
# other       -> produce nothing
# 
# $Id: ids.awk,v 1.3 1995-02-02 20:21:59 d3g681 Exp $

BEGIN {
		underline = "---------------------------------------------------------";
		printf("      subroutine util_version\n");
		printf("      write(6,*)\n");
		printf("      write(6,*) ' Software version information'\n");
		printf("      write(6,*) ' ----------------------------'\n");
}

/^module/			{
				 len = length($0);
				 printf("      write(6,*)\n");
				 printf("      write(6,*) ' %s'\n",$0);
				 printf("      write(6,*) ' %s'\n",substr(underline,1,len));
				}

/\$\I\d: [^\n]*Exp \$/		{
				 i = index($0, "$Id: ") + 5;
				 j = index($0, "Exp $") - 1;
				 n = j - i + 1;
				 if (n > (72 - 20)) n = 72 - 20;
				 printf("      write(6,*) ' %s'\n",substr($0,i,n));
				}

END {
		printf("      write(6,*)\n");
		printf("      call util_flush(6)\n");
		printf("      end\n");
}

	


