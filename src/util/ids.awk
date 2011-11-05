#
# This script is used by the makefile in building util_version.f.
# It expects on standard input lines containing
#
# a CVS ID    -> produce a write statement with the ID ("Exp" at the end)
# a SVN ID    -> produce a write statement with the ID (no "Exp" at the end)
# module name -> produce a write statement for the module name
# prev_str    -> skip to eliminate duplicates
# other       -> produce nothing
# 
# $Id$

BEGIN {
		underline = "---------------------------------------------------------";
		printf("      subroutine util_version\n");
                printf("#include \"global.fh\"\n");
                printf("      if (ga_nodeid().eq.0) then\n");
		printf("      write(6,*)\n");
		printf("      write(6,*) ' Software version information'\n");
		printf("      write(6,*) ' ----------------------------'\n");
		prev_str = underline;
}

/^module/			{
				 len = length($0);
				 printf("      write(6,*)\n");
				 printf("      write(6,*) ' %s'\n",$0);
				 printf("      write(6,*) ' %s'\n",substr(underline,1,len));
				 prev_str = underline;
				}

0 != match($0,prev_str)		{
				  $0="";
				}

/\$\I\d: [^\n]*Exp \$/		{
				 i = index($0, "$Id: ") + 5;
				 j = index($0, "Exp $") - 1;
				 n = j - i + 1;
				 if (n > (72 - 20)) n = 72 - 20;
				 printf("      write(6,*) ' %s'\n",substr($0,i,n));
				 prev_str = substr($0,i,n);
				 $0=""
				}

/\$\I\d: [^\n]* \$/		{
				 i = index($0, "$Id: ") + 5;
				 n = index(substr($0,i), " $") - 1;
				 if (n > (72 - 20)) n = 72 - 20;
				 printf("      write(6,*) ' %s'\n",substr($0,i,n));
				 prev_str = substr($0,i,n);
				}

END {
		printf("      write(6,*)\n");
		printf("      call util_flush(6)\n");
                printf("      endif\n");
                printf("      call ga_sync()\n");
		printf("      end\n");
}

	


