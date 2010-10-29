/*$Id$*/
#include <floatingpoint.h>
#include <stdio.h>
#include <signal.h>
#include <ucontext.h>

static int nunderflow;
static int ninexact;
static int ninvalid; 

static void catchunderflow(int sig, siginfo_t *sip, ucontext_t *uap)
{
    nunderflow++;
    printf(" !! #underflows = %d\n", nunderflow);
}

static void catchinexact(int sig, siginfo_t *sip, ucontext_t *uap)
{
    ninexact++;
    printf(" !! #inexacts = %d\n", ninexact);
}

static void catchinvalid(int sig, siginfo_t *sip, ucontext_t *uap)
{
    ninvalid++;
    printf(" !! #invalids = %d\n", ninvalid);
}

void ieeetrap_()
{
    (void) ieee_handler("set","common", SIGFPE_ABORT); 

/*
  (void) ieee_handler("set","inexact", SIGFPE_IGNORE);
  (void) ieee_handler("set","underflow", SIGFPE_IGNORE);
  (void) ieee_handler("set","invalid", SIGFPE_IGNORE);
*/

/* (void) ieee_handler("set","inexact", catchinexact);
 (void) ieee_handler("set","invalid", catchinvalid); */

 (void) ieee_handler("set","common", SIGFPE_ABORT);

    /* (void) ieee_handler("set","underflow", catchunderflow); */

}

