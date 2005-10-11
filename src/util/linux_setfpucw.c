/*
 $Id: linux_setfpucw.c,v 1.12 2005-10-11 23:53:46 edo Exp $
 */
#include <stdio.h>
#ifdef __CYGWIN__
#include <mingw/fenv.h>
#elif  __FreeBSD__
#include <ieeefp.h>
#else
#define __USE_GNU
#include <fenv.h>
#endif
/*#define FPSWAMOD /* this modifies fpswa behavior on ia64 */

void linux_trapfpe_(void) { 
int retval;
#if defined(LINUXIA64) && defined(FPSWAMOD)
#define FP_QUIET /* this prevents fp assist fault messages in syslog */
#include <linux/prctl.h>
//#include </usr/src/linux/include/linux/prctl.h>
/* grabbed from ftp://linux.hpl.hp.com/pub/linux-ia64/prctl-1.3.tar.gz  */
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
#endif

#elif  __FreeBSD__
    retval=fpsetmask(FP_X_DZ|FP_X_INV|FP_X_OFL);
#else
            retval = feenableexcept(FE_OVERFLOW | FE_DIVBYZERO|FE_INVALID);
#endif
        if (retval == -1) fprintf(stderr, "Failed to enable sigfpe  \n");
 }

