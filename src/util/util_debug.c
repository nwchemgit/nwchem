#if defined(CRAY) || defined(WIN32)

#include "typesf2c.h"

#ifdef CRAY
#define util_debug_ UTIL_DEBUG
#include <fortran.h>
#endif
#ifdef WIN32
#define util_debug_ UTIL_DEBUG
#endif

void FATR util_debug_(Integer *rtdb)
{
  ga_error("Don't know how to debug on this machine", 0);
}

#else

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "macdecls.h"
#include "global.h"
#include "rtdb.h"

#include "typesf2c.h"

void FATR util_debug_(Integer *rtdb)
{
  int i;
  pid_t child;
  char *argv[20];
  char display[256], path[256], title[256], pid[256], xterm[256];

  sprintf(title, "Debugging NWChem process %d", ga_nodeid_());
  sprintf(pid, "%d", getpid());

  if (!rtdb_get(*rtdb, "dbg:path", MT_CHAR, sizeof(path), path)) 
    strcpy(path, "nwchem");
  if (!rtdb_get(*rtdb, "dbg:xterm", MT_CHAR, sizeof(xterm), xterm)) 
    xterm[0] = 0;
  if (!rtdb_get(*rtdb, "dbg:display", MT_CHAR, sizeof(display), display)) {
    char *disp = getenv("DISPLAY");
    if (!disp)
      disp = "unix:0.0";
    strcpy(display, disp);
  }

  argv[1] = "-T";
  argv[2] = title;
  argv[3] = "-display";
  argv[4] = display;
  argv[5] = "-e";
#ifdef SOLARIS
  argv[6] = "dbx";
  argv[7] = path;
  argv[8] = pid;
  argv[9] = 0;
  if (!xterm[0])
    strcpy(xterm, "/usr/openwin/bin/xterm");
#elif defined(AIX) || defined(IBM) || defined(IBMSP) || defined(LAPI)
  argv[6] = "dbx";
  argv[7] = "-f";
  argv[8] = "-a";
  argv[9] = pid;
  argv[10] = 0;
  if (!xterm[0])
    strcpy(xterm, "/usr/bin/X11/xterm");
#elif defined(LINUX)
  argv[6] = "gdb";
  argv[7] = path;
  argv[8] = pid;
  argv[9] = 0;
  if (!xterm[0])
    strcpy(xterm, "/usr/X11R6/bin/xterm");
#elif defined(HPUX)
  argv[6] = "gdb";
  argv[7] = path;
  argv[8] = pid;
  argv[9] = 0;
  if (!xterm[0])
    strcpy(xterm, "/usr/bin/X11/xterm");
#elif defined(SGI_N32) || defined(SGI) || defined(SGIFP)
  argv[6] = "dbx";
  argv[7] = "-p";
  argv[8] = pid;
  argv[9] = 0;
  if (!xterm[0])
    strcpy(xterm, "/usr/bin/X11/xterm");
#else
  ga_error("Don't know how to debug on this machine", 0);
#endif

  argv[0] = xterm;
  if (ga_nodeid_() == 0) {
    int i;
    printf("\n Starting xterms with debugger using command\n\n    ");
    for (i=0; argv[i]; i++)
      printf("'%s' ", argv[i]);
    printf("\n\n");
    fflush(stdout);
  }
  ga_sync_();

  child = fork();

  if (child < 0) {
    ga_error("util_debug: fork failed?", 0);
  }
  else if (child > 0) {
    sleep(5);			/* Release cpu while debugger starts*/
  }
  else {
    execv(xterm, argv);
    perror("");
    ga_error("util_debug: execv of xterm with debugger failed", 0);
  }
}
#endif

