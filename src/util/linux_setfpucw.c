/*
 $Id: linux_setfpucw.c,v 1.8 2003-09-11 00:10:03 edo Exp $
 */
#include <stdio.h>
#ifdef __CYGWIN__
#include <mingw/fenv.h>
#else
#include <fenv.h>
#endif
#define FTZ 0x8000
#define FP_QUIET

void linux_trapfpe_(void) { 
#ifdef LINUXIA64
//#include <linux/prctl.h>
#include </usr/src/linux/include/linux/prctl.h>
/* grabbed from ftp://linux.hpl.hp.com/pub/linux-ia64/prctl-1.3.tar.gz  */
int retval;
/* this causes a Sigbus for each unaligned access */
    retval = prctl(PR_SET_UNALIGN,  PR_UNALIGN_SIGBUS);
	if (retval == -1) fprintf(stderr, "Failed to sigbus align \n");
#ifdef PR_SET_FPEMU
#ifdef FP_QUIET
/* this causes no print for each fpswa */
    retval = prctl(PR_SET_FPEMU,  PR_FPEMU_NOPRINT);
#else
/* this causes a SigFPE for each fpswa */
    retval = prctl(PR_SET_FPEMU,  PR_FPEMU_SIGFPE);
#endif
        if (retval == -1) fprintf(stderr, "Failed to sigfpe fpswa \n");
#endif

#else
            feenableexcept(FE_OVERFLOW | FE_DIVBYZERO);
#endif
 }

