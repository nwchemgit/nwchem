/*$Id: ieeetrap.c,v 1.2 1995-02-02 17:51:26 d3g681 Exp $*/
#include <floatingpoint.h>
#include <stdio.h>
#include <signal.h>

static void catchit()
{
  printf("!!  Floating point interrupt caught  !!\n");
  fflush(stdout);
/*  (void) signal(SIGIOT, SIG_DFL);*/
  abort();
}

void ieeetrap_()
{
  (void) ieee_handler("set","common", SIGFPE_ABORT);

/*
  (void) ieee_handler("set","inexact", SIGFPE_IGNORE);
  (void) ieee_handler("set","underflow", SIGFPE_IGNORE);
  (void) ieee_handler("set","invalid", SIGFPE_IGNORE);
*/

/*
 (void) ieee_handler("set","inexact", catchit);
 (void) ieee_handler("set","underflow", catchit);
 (void) ieee_handler("set","invalid", catchit);
*/

}
