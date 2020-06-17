#include <stdio.h>
#ifdef __CYGWIN__
#include <fenv.h>
#elif  __FreeBSD__
#include <ieeefp.h>
#else
#define __USE_GNU
#include <fenv.h>
#endif

void linux_trapfpe_(void)
{
    int retval;
#if  __FreeBSD__
    retval=fpsetmask(FP_X_DZ|FP_X_INV|FP_X_OFL);
#else
    retval = feenableexcept(FE_OVERFLOW | FE_DIVBYZERO|FE_INVALID);
#endif
    if (retval == -1) fprintf(stderr, "Failed to enable sigfpe  \n");
 }

