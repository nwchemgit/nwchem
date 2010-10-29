/*
 $Id$
 */
#if (defined(__ppc__) || defined(__ppc64__)) 
/* from http://developer.apple.com/documentation/Performance/Conceptual/Mac_OSX_Numerics/Mac_OSX_Numerics.pdf */
#define fegetenvd(x) asm volatile("mffs %0" : "=f" (x)); 
#define fesetenvd(x) asm volatile("mtfsf 255,%0" : : "f" (x)); 
enum { 
  FE_ENABLE_INEXACT = 0x00000008, 
  FE_ENABLE_DIVBYZERO = 0x00000010, 
  FE_ENABLE_UNDERFLOW = 0x00000020, 
  FE_ENABLE_OVERFLOW = 0x00000040, 
  FE_ENABLE_INVALID = 0x00000080, 
  FE_ENABLE_ALL_EXCEPT = 0x000000F8 };
typedef union { 
  struct { 
    unsigned long hi; 
        unsigned long lo; 
  } i; 
  double d; 
} hexdouble; 

void macx_trapfpe_(void) { 
  hexdouble t; 
  fegetenvd(t.d); 
  /* Enable hardware trapping for all exceptions */ 
  //  t.i.lo |= FE_ENABLE_ALL_EXCEPT; 
  t.i.lo = FE_ENABLE_OVERFLOW | FE_ENABLE_DIVBYZERO|FE_ENABLE_INVALID;
  fesetenvd(t.d); 
 }
#elif (defined (__i386__) || defined( __x86_64__ ))
#include <fenv.h>
#include <stdio.h>

#define _FPU_MASK_IM  0x01
#define _FPU_MASK_ZM  0x04
#define _FPU_MASK_OM  0x08
#define _FPU_RESERVED 0xF0C0 

#ifndef _FE_DIVBYZERO
#define _FE_DIVBYZERO FE_DIVBYZERO
#endif
#ifndef _FE_OVERFLOW
#define _FE_OVERFLOW  FE_OVERFLOW
#endif
#ifndef _FE_INVALID
#define _FE_INVALID FE_INVALID
#endif
typedef unsigned int fpu_control_t __attribute__ ((__mode__ (__HI__)));
extern fpu_control_t __fpu_control;

#define _FPU_GETCW(cw) asm volatile ("fnstcw %0" : "=m" (*&cw))
#define _FPU_SETCW(cw) asm volatile ("fldcw %0" : : "m" (*&cw))
#define _FPU_GETMXCSR(cw_sse) asm volatile ("stmxcsr %0" : "=m" (cw_sse))
#define _FPU_SETMXCSR(cw_sse) asm volatile ("ldmxcsr %0" : : "m" (cw_sse))

int macx_trapfpe_()
{
  fpu_control_t mode, mode_sse;

  _FPU_GETCW (mode) ;
  mode &= (_FPU_RESERVED | _FE_DIVBYZERO |  _FE_OVERFLOW |  _FE_INVALID) ;
  _FPU_SETCW (mode) ;

  _FPU_GETMXCSR (mode_sse) ;
  mode_sse &= (0xFFFF0000 | (_FPU_RESERVED | _FE_DIVBYZERO | _FE_OVERFLOW | _FE_INVALID) <<7);
  _FPU_SETMXCSR (mode_sse) ;

  return 1 ;
}
#else
#error arch not ready
#endif
