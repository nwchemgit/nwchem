/*---------------------------------------------------------*\
$Id$
\*---------------------------------------------------------*/

/*------------------------------------------------------------*\
  s e e t e x

  code to parse tex directives from specific comment lines in
  NWChem fortran source code, fortran include files, and latex source with
  seetex block begin and end statements.  

  Target syntax for FORTRAN Latex comments
  *:tex-         means print this line (stripping *:tex-)
  *:tex-\begin   means print this line and all that follow until the 
                 termination delimeter is found.  strip all *:tex
  *:tex-\end     terminate tex block of information. print this line
                 stripping all *:tex strings
  these are for allowing tex files to be read as well
  *:tex-begin    means print lines after this one until *:tex-end
                 is found
  *:tex-end      terminates *:tex-end
  *:tex-cvsid    inserts a latex specific cvs id tag in the output file.
  c:tex- is equivalent to *:tex in all cases
  !:tex- is equivalent to *:tex but it applies to only single line 
         Latex strings.

  Assumptions:
  1) files have the appropriate suffix ".F", ".f", ".fh", and ".th"
     ".th" -> seetex-ed tex include file with blocked begin and ends:
     *:tex-begin
     %normal latex source
     *:tex-end
  2) output file has the appropriate ".tex" suffix
  3) output file will be appended to if it exists and created if not
  4) performance of seetex is not critical
  5) files are parsed in command line order.

  Future Enhancements:
  1) add c code?? (hard)

  Usage: seetex file1.F [file2.f file3.fh ....] file.tex
  Flags: -h or -help

  Designed and Written by:
     Ricky A. Kendall
     High Performance Computational Chemistry Group
     Environmental Molecular Sciences Laboratory
     Pacific Northwest National Laboratory 
     Richland, WA 99352-0999

  Dates:
    Initial Version:    April 1996.
    Usage update   :    September 1996.
\*------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "seetex.h"

int main (int argc, char *argv[])
{
  int i;
  char *file_tex;
  char *file_source;
  FILE *texid;
  FILE *srcid;
    

  if (DEBUG_MODE >= 10)
    {
      for (i=0;i<argc;i++)
	(void) printf(" seetex argument: %d  is <%s>\n",i,argv[i]);
      (void) printf("\n");
    }

  if (argc == 1) 
      {
	  (void) printf(" zero arguments: print usage/syntax only \n\n");
	  (void) print_syntax();
	  (void) usage(0);
      }
  for (i=0;i<argc;i++) 
      {
	  if (! strncmp(argv[i],"-h",(size_t) 2)) 
	      {
		  (void) printf(" [-h|-help] flag found: print usage/syntax only\n\n");
		  (void) print_syntax();
		  (void) usage(0);
	      }
      }
  file_tex = argv[(argc-1)];
  if (DEBUG_MODE >= 10) (void) printf(" seetex output file <%s>\n\n",file_tex);
  
  if (is_tex_file(file_tex))  /* is file *.tex ? */
    {
/* check existance and/or permissions */
      if (access(file_tex,F_OK))
	{if (DEBUG_MODE >= 10) (void) printf("Creating %s\n",file_tex);}
      else
	{if (DEBUG_MODE >= 10) (void) printf("Appending to %s\n",file_tex);}
      
      texid = fopen(file_tex,"a");
    }
  else
    (void) seetex_error("tex file not specified in proper order",-1);

  for(i=1;i<(argc-1);i++)
    { 
      file_source = argv[i];
      if (is_fortran_file(file_source) && (!access(file_source,R_OK)))
	{
	  if (DEBUG_MODE >= 10) (void) printf("fortran or seetex include file %s is readable\n",file_source);
	  srcid = fopen(file_source,"r");
	}
      else
	(void) seetex_error
	  (" fortran source or seetex include file does not exist or \n              does not have the proper '.F', '.f', '.fh' suffix",-2);
      if (DEBUG_MODE >= 10)
	{
	  if (is_fortran_file(file_source))
	    (void) printf("%s is a fortran or seetex include file\n",file_source);
	  else
	    (void) printf("%s is NOT a fortran or seetex include file\n",file_source);
	}
      (void) printf(" tex output from <%s> written to <<%s>>\n",file_source,file_tex);
      if (! seetex_process(srcid,texid))
	(void) seetex_error(" error processing fortran source or seetex include file\n",1);
      (void) fflush(texid);
      (void) fclose(srcid);
    }

  (void) fflush(texid);
  (void) fclose(texid);
/* return no error condition */
  return FALSE;
}

int get_line(FILE *file, char *buf, int size)
{
  char *ptmp;
  int length;
  int i;

  ptmp = fgets(buf, (size-1), file);

  if (ptmp == (char *) NULL) 
    { 
      if (DEBUG_MODE >=10)(void) printf("get_line EOF\n");
      return FALSE;}

  if (DEBUG_MODE >= 10)
    (void) printf("before strip:\nbuf  is <%s>\nptmp is <%s>\n",buf,ptmp);

/* replace '\n' with ' ' in string */  

  length = strlen(ptmp);
  
  for (i=0;i<length;i++)
    {
      if (*ptmp == '\n') *ptmp = ' ';
      ptmp++;
    }
  if (DEBUG_MODE >= 10)
    (void) printf("after strip:\nbuf  is <%s>\nptmp is <%s>\n",buf,ptmp);

  return TRUE;
}

int seetex_end(char *line)
{
  if (! strncmp(line,"*:tex-\\end",(size_t) 10)) 
    return TRUE;
  if (! strncmp(line,"c:tex-\\end",(size_t) 10)) 
    return TRUE;
  if (! strncmp(line,"*:tex-end",(size_t) 9)) 
    return TRUE;
  if (! strncmp(line,"c:tex-end",(size_t) 9)) 
    return TRUE;
  return FALSE;
}
int seetex_begin(char *line)
{
  if (! strncmp(line,"*:tex-\\begin",(size_t) 12)) 
    return TRUE;
  if (! strncmp(line,"c:tex-\\begin",(size_t) 12)) 
    return TRUE;
  if (! strncmp(line,"*:tex-begin",(size_t) 11)) 
    return TRUE;
  if (! strncmp(line,"c:tex-begin",(size_t) 11)) 
    return TRUE;
  return FALSE;
}
void put_line_strip(char *line,FILE *id)
{
  char *ptmp;

  ptmp = line;
  if (!strncmp(ptmp,"c:tex-begin",(size_t) 11)) return;
  if (!strncmp(ptmp,"c:tex-end",(size_t) 9)) return;
  if (!strncmp(ptmp,"*:tex-begin",(size_t) 11)) return;
  if (!strncmp(ptmp,"*:tex-end",(size_t) 9)) return;
  if (!strncmp(ptmp,"*:tex-cvsid",(size_t) 11)) 
    { (void) fprintf(id,"%%\n%% $");
      (void) fprintf(id,"Id$\n%%\n"); return;}
  if (!strncmp(ptmp,"c:tex-cvsid",(size_t) 11)) 
    { (void) fprintf(id,"%%\n%% $");
      (void) fprintf(id,"Id$\n%%\n"); return;}
  if (!strncmp(ptmp,"*:tex-",(size_t) 6)) ptmp += 6;
  if (!strncmp(ptmp,"c:tex-",(size_t) 6)) ptmp += 6;
  (void) fprintf(id,"%s\n",ptmp);
  (void) fflush(id);
}

void put_line_strip_check(char *line,FILE *id)
{
  char *ptmp;
  int length,i;
  int print_it=0;

  ptmp = line;

  if (!strncmp(ptmp,"c:tex-begin",(size_t) 11)) return;
  if (!strncmp(ptmp,"c:tex-end",(size_t) 9)) return;
  if (!strncmp(ptmp,"*:tex-begin",(size_t) 11)) return;
  if (!strncmp(ptmp,"*:tex-end",(size_t) 9)) return;
  if (!strncmp(ptmp,"*:tex-cvsid",(size_t) 11)) 
    { (void) fprintf(id,"%%\n%% $");
      (void) fprintf(id,"Id$\n%%\n"); return;}
  if (!strncmp(ptmp,"c:tex-cvsid",(size_t) 11)) 
    { (void) fprintf(id,"%%\n%% $");
      (void) fprintf(id,"Id$\n%%\n"); return;}
  if (!strncmp(ptmp,"*:tex-",(size_t) 6)) 
      {ptmp += 6;print_it++;}
  if (!strncmp(ptmp,"c:tex-",(size_t) 6))
      {ptmp += 6;print_it++;}
  if (print_it) 
      {
	  (void) fprintf(id,"%s\n",ptmp);
	  (void) fflush(id);
	  return;
      }
  length = strlen(ptmp);
  for (i=0;i<length;i++)
    {
      if (!strncmp(ptmp,"!:tex-",(size_t) 6))
	{ptmp += 6;print_it++;
	 break;}
      ptmp++;
    }
  if (print_it) (void) fprintf(id,"%s\n",ptmp);
  (void) fflush(id);
}

int is_tex_file(char *filename)
{
  char *ptr;
  int length;
  int return_value;

  length = strlen(filename);
  if (length < 5) 
    ptr = filename;
  else
    ptr = filename + length - 4;

  if (DEBUG_MODE >= 10)
    (void) printf(".... is_tex_file .... \nfilename: <%s>\nptr     : <%s>\n",
		  filename,ptr);
  if (length > 4 && !strcmp(ptr,".tex"))
    return_value = TRUE;
  else
    return_value = FALSE;  

  if (DEBUG_MODE >= 10) 
    (void) printf("is_tex_file: return_value is %d\n\n",
		  return_value);

  return return_value;
}

void seetex_error(char *errmsg, int error_code)
{
  (void)printf("\nseetex:error:\n");
  (void)printf("seetex:error:%s \n",errmsg);
  (void)printf("seetex:error:\n\n");
  if (error_code < 0)
    (void) usage(error_code);
  else
    (void) exit((int) error_code);
}

void print_syntax()
{
    (void) printf("s e e t e x\n");
    (void) printf("\n");
    (void) printf("code to parse tex directives from specific comment lines in\n");
    (void) printf("NWChem fortran source code, fortran include files, and latex source with\n");
    (void) printf("seetex block begin and end statements.  \n");
    (void) printf("\n");
    (void) printf("Target syntax for FORTRAN Latex comments\n");
    (void) printf("*:tex-         means print this line (stripping *:tex-)\n");
    (void) printf("*:tex-\\begin   means print this line and all that follow until the \n");
    (void) printf("termination delimeter is found.  strip all *:tex\n");
    (void) printf("*:tex-\\end     terminate tex block of information. print this line\n");
    (void) printf("                stripping all *:tex strings\n");
    (void) printf("these are for allowing tex files to be read as well\n");
    (void) printf("*:tex-begin    means print lines after this one until *:tex-end\n");
    (void) printf("               is found\n");
    (void) printf("*:tex-end      terminates *:tex-end\n");
    (void) printf("c:tex- is equivalent to *:tex in all cases\n");
    (void) printf("!:tex- is equivalent to *:tex but it applies to only single line \n");
    (void) printf("       Latex strings.\n");
    (void) printf("\n");
    (void) printf("Assumptions:\n");
    (void) printf("1) files have the appropriate suffix \".F\", \".f\", \".fh\", and \".th\"\n");
    (void) printf("  \".th\" -> seetex-ed tex include file with blocked begin and ends:\n");
    (void) printf("   *:tex-begin\n");
    (void) printf("   %%normal latex source\n");
    (void) printf("   *:tex-end\n");
    (void) printf("2) output file has the appropriate \".tex\" suffix\n");
    (void) printf("3) output file will be appended to if it exists and created if not\n");
    (void) printf("4) performance of seetex is not critical\n");
    (void) printf("5) files are parsed in command line order.\n");
    (void) printf("\n");
    (void) printf("Future Enhancements:\n");
    (void) printf("1) add c code?? (hard)\n");
    (void) printf("\n");
    (void) printf("Usage: seetex file1.F [file2.f file3.fh ....] file.tex\n");
    (void) printf("Flags: -h or -help\n");
    (void) printf("\n");
    (void) printf("Designed and Written by:\n");
    (void) printf("   Ricky A. Kendall\n");
    (void) printf("   High Performance Computational Chemistry Group\n");
    (void) printf("   Environmental Molecular Sciences Laboratory\n");
    (void) printf("   Pacific Northwest National Laboratory \n");
    (void) printf("   Richland, WA 99352-0999\n");
    (void) printf("\n");
    (void) printf("Dates:\n");
    (void) printf("   Initial Version:    April 1996.\n");
    (void) printf("   Usage update   :    September 1996.\n");
    (void) printf("\n");
}

void usage(int error_code)
{
  (void) printf("Usage: seetex file1.F [file2.f file3.fh ...] file.tex\n\n");
  (void) printf("Flags: -h or -help prints the usage and syntax info\n");
  (void) printf("\nNote: tex file must be specified last with the .tex extension\n");
  (void) exit((int) error_code);
}
int is_fortran_file(char *filename)
{
  char *ptr;
  int length;
  int stat_F;
  int stat_f;
  int stat_fh;
  int stat_th;
  int return_value;

  length = strlen(filename);
  ptr = filename + length;
  while (TRUE)
    if (!strncmp(ptr,".",(size_t) 1))
      break;
    else
      ptr--;

  if (DEBUG_MODE >= 10) 
    (void) printf(" is_fortran_file: ptr is <%s> \n",ptr);

  stat_F  = ! strcmp(ptr,".F");
  stat_f  = ! strcmp(ptr,".f");
  stat_fh = ! strcmp(ptr,".fh");
  stat_th = ! strcmp(ptr,".th");

  if (stat_F || stat_f || stat_fh || stat_th) 
    return_value = TRUE;
  else
    return_value = FALSE;
    
  if (DEBUG_MODE >= 10) 
    (void) printf("is_fortran_file: return_value is %d\n\n",
		  return_value);

  return return_value;
}

int seetex_process(FILE *srcid, FILE *texid)
{
  char line[MAX_LINE_LEN];
  char *pline = (char *)NULL;
  int eof_found = 0;
  int in_seetex_block = 0;
  int line_count = 0;
  
  pline = &line[0];
  
  while (! eof_found)
    {
      if (DEBUG_MODE >= 10)
	(void) printf("seetex_process:<%d>\n",line_count);	
      eof_found = ! get_line(srcid,line,MAX_LINE_LEN);
      if (DEBUG_MODE >= 10)
	(void) printf("seetex_process:<%s>\n",pline);
      if (eof_found)
	{ if (DEBUG_MODE >= 10)
	    (void) printf("EOF found\n");
	  break;}
      
      line_count++;
      if (DEBUG_MODE >= 10)
	(void) printf("seetex_process: read line %d \n",line_count);
      if (seetex_end(pline)) 
	if (in_seetex_block > 0)
	  in_seetex_block--;
	else
	  {
	    (void) printf("Found *:tex-{\\}end when not in a tex block\n");
	    return FALSE;
	  }
      
      if (seetex_begin(pline)) 
/*	if (in_seetex_block == 0)*/
	  in_seetex_block++;
/*	else
	  {
	    (void) printf("Found *:tex-{\\}begin while in a tex block\n");
	    return FALSE;
	  }*/

      if (in_seetex_block)
	(void) put_line_strip(pline,texid);
      else
	(void) put_line_strip_check(pline,texid);
    }
  if (in_seetex_block == 0)
    return TRUE;
  else
    {
      (void) printf("EOF found but still in tex block\n");
      return FALSE;
    }
}


