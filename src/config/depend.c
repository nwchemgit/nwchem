/*
 $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

extern char *strdup(const char *);
char *makefile;
char backup[] = "makefile.bak";

void copy_truncate_makefile(const char *backup)
{
    char buf[8192];
    FILE *in;
    int i, j, ninbuf;
    char line[] = 
	"# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY";
    int len = strlen(line);

    if (!(in = fopen(backup,"r"))) {
	fprintf(stderr,"failed to open %s for reading\n",backup);
	perror(" ");
	exit(1);
    }
    
    ninbuf = 0;
    while ((i = getc(in)) != EOF) {
	if (i == '\n') {
	    if (ninbuf && (strncmp(buf,line,len) == 0)) {
		ninbuf = 0;
		break;		/* Truncate here */
	    }
	    else {
		for (j=0; j<ninbuf; j++)
		    putchar(buf[j]);
		putchar('\n');
		ninbuf = 0;
	    }
	}
	else {
	    if (ninbuf >= sizeof(buf)) {
		fprintf(stderr,"dependencies: makefile line too long\n");
		exit(1);
	    }
	    buf[ninbuf++] = i;
	}
    }
    
    if (ninbuf) { /* Last line missed a carriage return */
	for (j=0; j<ninbuf; j++)
	    putchar(buf[j]); 
	putchar('\n');
    }
    printf("%s\n",line);
    fclose(in);
}
		    
void copy_file(const char *input, const char *output)
{
    char buf[8192];
    int nread;
    FILE *in, *out;

    if (!(in = fopen(input,"r"))) {
	fprintf(stderr,"failed to open %s for reading\n",input);
	perror(" ");
	exit(1);
    }
    if (!(out = fopen(output,"w"))) {
	fprintf(stderr,"failed to open %s for writing\n",output);
	perror(" ");
	exit(1);
    }

    while ((nread = fread(buf, 1, sizeof(buf), in)))
	if (nread != fwrite(buf, 1, nread, out)) {
	    fprintf(stderr,"failed writing %s\n",output);
	    perror(" ");
	    exit(1);
	}
    fclose(in);
    fclose(out);
}

void restore_and_error_exit()
{
    fclose(stdout);		/* was going to makefile */
    copy_file(backup, makefile);
    exit(1);
}

int compar(const void *a, const void *b)
{
    return strcmp(*((char **) a), *((char **) b));
}

void skip_to_eol(FILE *file)
{
    int i;

    while ((i = getc(file)) != EOF)
	if (i == '\n')
	    break;
}

void skip_white_space(FILE *file)
{
    int i;

    while ((i = getc(file)) != EOF)
	if ((i != ' ') && (i != '\t')) {
	    ungetc(i,file);
	    break;
	}
}

char *include_directive(FILE *file)
{
#define MAXBUF 256
    char tmp[MAXBUF];
    int n = 0;
    int i;

    skip_white_space(file);	/* Allow for #   include */
    
    if ((i = getc(file)) == EOF) return 0;

    if (i != 'i')		/* Not an include directive */
	return 0;

    while ((i = getc(file)) != EOF)  {/* Find first quote */
	if (i == '"' || i =='<')
	    break;
	if (i == '\n') return 0; /* Don't wrap lines in bad code */
    }

    if (i == EOF) return 0;

    while ((i = getc(file)) != EOF)  {/* Find first quote */
	if (i == '\n')
	    return 0;		/* Mangled code */
	else if (i == '"') 
	    break;
	if (n >= (MAXBUF-1)) {
	    fprintf(stderr,"dependencies: include file name too long\n");
	    exit(1);
	}
	tmp[n++] = i;
    }
    tmp[n] = 0;
    return strdup(tmp);
}

int main(int argc, const char *argv[])
{
#define MAXINCDIR 1024
#define MAXFILE 16384
    const char *incdirlist[MAXINCDIR];
    char *filelist[MAXFILE];
    int nincdir, nfiles;
    DIR *dirp;
    struct dirent *ent;

    /* This code used to append dependencies to the makefiles.
       Now it generates dependencies on STDOUT and does not
       mess with the makefiles at all.  However, the old stuff
       is left in here just in case.   To recover the old version
       recompile with -DAPPENDMAKEFILE */

#ifdef ADDPENDMAKEFILE
    /* Figure out the name of the makefile */

    if (access("makefile",R_OK|W_OK) == 0)
	makefile = "makefile";
    else if (access("Makefile",R_OK|W_OK) == 0)
	makefile = "Makefile";
    else {
	fprintf(stderr,"dependencies: makefile/Makefile not found\n");
	exit(1);
    }

    /* Make a backup copy */

    copy_file(makefile, backup);

    /* Reopen the makefile as standard output, copy the backup
       over it truncating before the magic line.  Write the magic
       line and then do the dependency analysis */

    if (!freopen(makefile, "w", stdout)) {
	perror("dependencies: failed reopening stdout");
	exit(1);
    }

    copy_truncate_makefile(backup);
#endif

    /* Extract list of include paths */

    nincdir = 1; incdirlist[0] = ".";

    while (argc > 1) {
	argc--;
	if (strncmp("-I",argv[argc],2) == 0) {
	    if (nincdir >= MAXINCDIR) {
		fprintf(stderr,"dependencies: too many include directories\n");
		restore_and_error_exit();
	    }
	    incdirlist[nincdir++] = argv[argc]+2;
	}
	else {
	    fprintf(stderr,"usage: dependencies [-Idir] [-Idir] ...\n");
	    restore_and_error_exit();
	}
    }

    /* Generate list of .F or .c files in the current directory 
       ... don't acutually need to store the list but it is convenient */

    if (!(dirp = opendir("."))) {
	perror("dependencies: failed to open current directory");
	restore_and_error_exit();
    }

    nfiles = 0;
    while ((ent = readdir(dirp))) {
	char *name = ent->d_name;
	int len = strlen(name);
	if ( (strncmp(".c",name+len-2,2) == 0) ||
             (strncmp(".F",name+len-2,2) == 0)) {
	    if (nfiles >= MAXFILE) {
		fprintf(stderr,"dependencies: too many files\n");
		restore_and_error_exit();
	    }
	    filelist[nfiles++] = strdup(name);
	}
    }

    (void) closedir(dirp);

    /* Loop thru files */

    while (nfiles--) {
	char *name = filelist[nfiles];
	char *objname = strdup(name);
	FILE *file;
	int i;
#define MAXINCFILE 8192
	char *incfiles[MAXINCFILE];
	int nincfile;
	int prev;

	objname[strlen(objname)-1] = 'o';

	/* Form sorted list of unique include files */

	/* printf("%d\t%s\n",nfiles,name); */
	if (!(file = fopen(name,"r"))) {
	    fprintf(stderr,"dependencies: failed to open %s : ",name);
	    perror(" ");
	    restore_and_error_exit();
	}

	nincfile = 0;
	while ((i = getc(file)) != EOF) {
	    char *tmp;
/************************************************************************\
Original code:
	    if (i == '#') 
		if ((tmp = include_directive(file))) {
		    if (nincfile >= MAXINCFILE) {
			fprintf(stderr,"dependencies: too many includes\n");
			restore_and_error_exit();
		    }
		    incfiles[nincfile++] = tmp;
		}
	    else
	       skip_to_eol(file);
\***********************************************************************/

/* makeing sure that # is the first character takes care of fortran comments
   but it does not take care of c style comments or includes that are dependent
   upon cpp flags.  The latter is the difficult one.
*/
	    if (i == '#') 
		{
		    if ((tmp = include_directive(file))) {
			if (nincfile >= MAXINCFILE) {
			    fprintf(stderr,"dependencies: too many includes\n");
			    restore_and_error_exit();
			}
			incfiles[nincfile++] = tmp;
		    }
		}
/*          else                  *\
  removed this "else" because 
  include_directive stops after 
  the second '"' and not the end 
  of the line.  Therefore all 
  lines parsed must go to EOL.
\*                                */
               skip_to_eol(file);
	}
	fclose(file);

	if (nincfile == 0)	/* If no included files skip this file */
	    continue;

	printf("$(LIBRARY_PATH)(%s):\t",objname);
	
	qsort(incfiles, nincfile, sizeof(char *), compar);
	
	prev = 0;		/* Remove and free duplicates */
	for (i=1; i<nincfile; i++) {
	    if (strcmp(incfiles[i],incfiles[prev]) == 0) {
		free(incfiles[i]);
		incfiles[i] = 0; /* For paranoia */
	    }
	    else
		incfiles[++prev] = incfiles[i];
	}
	nincfile = prev + 1;
	
	/* Loop thru the include files figuring out where each is */
	
    while (nincfile--) {
        char *incname = incfiles[nincfile];
        char path[256];

        for (i=0; i<nincdir; i++) {
            (void) sprintf(path, "%s/%s", incdirlist[i], incname);
            if (access(path, R_OK) == 0) {
                break;
            }
            path[0] = 0;
        }
        if (!path[0]) {
            (void) sprintf(path, "$(INCDIR)/%s", incname);
        }
        /* Skip mpif.h header since it is an external header. */
        if (0 != strcmp(incname, "mpif.h")) {
            printf("%s ",path);
        }
        free(incname);
    }
	printf("\n");
	free(objname);
    }
    return 0;
}
