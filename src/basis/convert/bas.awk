BEGIN {element = "nosuchatom"; what="basis"; firstone=1;}

/^#/    {next;}  #Discard full line comments.  Too much work to treat correctly.  Just refer to the EMSL web page for full info.

/^END/  {next;}

/^[ \t]$/  {next;}

/^Effective Core Potentials/  {what="ecp"; next;}
/^--------/  {next;}

/^BASIS/  {what="basis"; type = "";
           if (($NF ~ /SPHERICAL/) || ($NF ~ /CARTESIAN/)) {
              type = $NF;
              $NF = "";
           } 
           $1=""; gsub("\" ","\"",$0); gsub("\"","",$0); gsub(" ","_",$0); basis=$0; element="nosuchatom";next;}

/[ \t]*[A-Z][a-z]*[ \t]+[SPDFGHInuL]+/  {
       if ((length($1) != length(element)) || (!match($1,element)) ) {
          element=$1;
          if(!firstone) printf("end\n"); firstone=0;
          printf("%s \"%s%s\"",what,element,basis);
          if (! (what ~ /ecp/)) printf(" %s",type);
          printf("\n");}
        }

    #   /^[ \t]*[0-9.]+[ \t]+[-+0-9.]+[ \t]*/ {if ($2 != 0.0) print; next;}
/^[ \t]*[0-9.]+[ \t]+[-+0-9.]+[ \t]*/ {print; next;}

        {printf("  ");print;}


END {printf("end\n");}
