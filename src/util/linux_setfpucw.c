/*
 $Id: linux_setfpucw.c,v 1.1 2001-01-18 20:40:10 edo Exp $
 */
#include <fpu_control.h>
/* static void __attribute__ ((constructor)) trapfpe (void)*/
  void  linux_trapfpe_ (void)
 {
 fpu_control_t cw = 
 _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM |_FPU_RC_ZERO); _FPU_SETCW(cw);
 }
  
