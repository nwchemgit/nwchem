/*
 $Id: macx_trapfpe.c,v 1.1 2005-01-12 00:06:06 edo Exp $
 */
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

