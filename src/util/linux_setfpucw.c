/*
 $Id: linux_setfpucw.c,v 1.2 2001-03-12 15:43:22 bjohnson Exp $
 */

#ifdef CYGNUS

/* Make this a dummy routine under CYGNUS since fpu_control.h doesn't exist */
void linux_trapfpe_(void) { }

#else

/* Regular LINUX case*/
#include <fpu_control.h>
/* static void __attribute__ ((constructor)) trapfpe (void)*/
  void  linux_trapfpe_ (void)
 {
 fpu_control_t cw = 
 _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM |_FPU_RC_ZERO); _FPU_SETCW(cw);
 }

#endif
