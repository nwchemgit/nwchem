/*
 $Id: linux_setfpucw.c,v 1.5 2003-04-23 01:37:07 edo Exp $
 */

#ifdef CYGNUS

/* Make this a dummy routine under CYGNUS since fpu_control.h doesn't exist */
void linux_trapfpe_(void) { }
#elif defined(LINUXIA64)
#include <fenv.h>
#include <stdio.h>
#include <linux/prctl.h>
//#include </usr/src/linux-2.4.19-hp2_pnnl6_Lv15irqstkM/include/linux/prctl.h>
/* grabbed from ftp://linux.hpl.hp.com/pub/linux-ia64/prctl-1.3.tar.gz  */
void linux_trapfpe_(void) { 
int retval;
/* this causes a Sigbus for each unaligned access */
    retval = prctl(PR_SET_UNALIGN,  PR_UNALIGN_SIGBUS);
	if (retval == -1) fprintf(stderr, "Failed to sigbus align \n");
#ifdef PR_SET_FPEMU
/* this causes a SigFPE for each fpswa */
    retval = prctl(PR_SET_FPEMU,  PR_FPEMU_SIGFPE);
	if (retval == -1) fprintf(stderr, "Failed to sigfpe fpswa \n");
#endif
}

#else

/* Regular LINUX case*/
#include <fpu_control.h>
#include <math.h>
#include <stdio.h>
/* static void __attribute__ ((constructor)) trapfpe (void)*/
  void  linux_trapfpe_ (void)
 {
   /* fpu_control_t cw = 
    (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(cw);*/
 fpu_control_t cw = 
#ifdef IFCLINUX
 _FPU_DEFAULT & ~(_FPU_MASK_ZM |_FPU_MASK_OM | _FPU_RC_ZERO); _FPU_SETCW(cw);
#else
 _FPU_DEFAULT & ~(_FPU_MASK_ZM |_FPU_MASK_OM | _FPU_RC_ZERO); _FPU_SETCW(cw);
 // _FPU_DEFAULT & ~(_FPU_MASK_IM |_FPU_MASK_ZM |_FPU_MASK_OM | _FPU_RC_ZERO); _FPU_SETCW(cw);
#endif
 }

#endif
