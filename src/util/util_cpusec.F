      double precision function util_cpusec()
      implicit none
c $Id: util_cpusec.F,v 1.3 1995-11-03 20:44:13 d3g681 Exp $
*
* routine to give cpu seconds since start of execution
*  delta times are handled by all applications
*
#if defined(IBM) || defined(SP1)
#define __UTIL_CPUSEC__
      double precision ibm_cputime
      external ibm_cputime
      util_cpusec = ibm_cputime()
#endif
#if defined(KSR)
#define __UTIL_CPUSEC__
*--KSR-start
      real ksr_cputime
      external ksr_cputime
c
      util_cpusec = ksr_cputime()
*--KSR-end
#endif
#if defined(SUN) || defined(SGI) || defined(SGITFP)|| defined(__PARAGON__)
#define __UTIL_CPUSEC__
      real tarray(2)
      real etime
      external etime
c
      util_cpusec = dble(etime(tarray))
c      
#endif
#if defined(DELTA)
#define __UTIL_CPUSEC__
      double precision dclock
      external dclock
c
      util_cpusec = dclock()
#endif
#if defined(LINUX)
#define __UTIL_CPUSEC__
      double precision linux_cputime
      external linux_cputime
c
      util_cpusec = linux_cputime()
#endif
*
#if !defined(__UTIL_CPUSEC__)
      util_cpusec = 0.0d00
#endif
      end
