      double precision function util_cpusec()
      implicit none
c $Id: util_cpusec.F,v 1.12 1999-11-13 03:24:14 bjohnson Exp $
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
#if defined(SUN) || defined(__PARAGON__) || defined(DECOSF) || defined(SOLARIS) || defined(HPUX)|| defined(FUJITSU_VPP)
#define __UTIL_CPUSEC__
      real*4 tarray(2)
      real*4 etime
      external etime
c
      util_cpusec = etime(tarray)
c      
#endif

#if defined(SGI)
#define __UTIL_CPUSEC__
#if defined(SGITFP)
      real*4 tarray(2)        ! nb: real == real*8 on TFP
      real*4 etime
      external etime
#else
      real tarray(2)
      real etime
      external etime
#endif
      util_cpusec = etime(tarray)
#endif  

#if defined(DELTA)
#define __UTIL_CPUSEC__
      double precision dclock
      external dclock
      double precision first
      logical ofirst
      data ofirst /.true./
c
      if (ofirst) then
	first = dclock()
	ofirst = .false.
      endif
c
      util_cpusec = dclock() - first
#endif
#if defined(LINUX)
#if defined(CYGNUS)
#else
#define __UTIL_CPUSEC__
      double precision linux_cputime
      external linux_cputime
c
      util_cpusec = linux_cputime()
#endif
#endif
#if defined(CRAY)
#define __UTIL_CPUSEC__
      double precision secondr
      external secondr
      double precision first
      logical ofirst
      save first, ofirst
      data ofirst /.true./
c
      if (ofirst) then
	first = secondr()
	ofirst = .false.
      endif
c
      util_cpusec = secondr() - first ! Actually wall time
#endif
*
#if !defined(__UTIL_CPUSEC__)
      util_cpusec = 0.0d00
#endif
      end


