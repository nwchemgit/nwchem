/*
 $Id: linux_setfpucw.c,v 1.3 2001-06-28 00:16:13 edo Exp $
 */

#ifdef CYGNUS

/* Make this a dummy routine under CYGNUS since fpu_control.h doesn't exist */
void linux_trapfpe_(void) { }
#elif defined(LINUXIA64)
#include <fenv.h>
# define FE_NONIEEE_ENV ((__const fenv_t *) 0xc009a04d0270037fUL)               
/* http://devresource.hp.com/devresource/Docs/TechTips/ia64linuxTips.html#ia64lptip2 */
void linux_trapfpe_(void) { 
    fesetenv (FE_NONIEEE_ENV);
  //fedisableexcept (FE_ALL_EXCEPT); 
//__asm__ __volatile__ (";; mov.m ar.fpsr=%0;;" :: "r"(st));

    //asm volatile ("mov ar.fpsr=%0" :: "r"(0x9804c0270037f));

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
