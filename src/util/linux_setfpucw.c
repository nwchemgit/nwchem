/*
 $Id: linux_setfpucw.c,v 1.4 2003-04-22 01:55:50 edo Exp $
 */

#ifdef CYGNUS

/* Make this a dummy routine under CYGNUS since fpu_control.h doesn't exist */
void linux_trapfpe_(void) { }
#elif defined(LINUXIA64)
#include <fenv.h>
#include <stdio.h>
#include <linux/prctl.h>
# define FE_NONIEEE_ENV ((__const fenv_t *) 0xc009a04d0270037fUL)               
/* http://devresource.hp.com/devresource/Docs/TechTips/ia64linuxTips.html#ia64lptip2 */
void linux_trapfpe_(void) { 
int retval;
//    fesetenv (FE_NONIEEE_ENV);
  //fedisableexcept (FE_ALL_EXCEPT); 
//__asm__ __volatile__ (";; mov.m ar.fpsr=%0;;" :: "r"(st));
/* this causes a Sigbus for each unaligned access */
    retval = prctl(PR_SET_UNALIGN,  PR_UNALIGN_SIGBUS);
	if (retval == -1) fprintf(stderr, "Failed to sigbus align \n");
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
