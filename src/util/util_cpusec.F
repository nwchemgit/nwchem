      double precision function util_cpusec()
      implicit none
c $Id$
*
* routine to give cpu seconds since start of execution
*  delta times are handled by all applications
*
#if defined(BGL) || defined(BGP) || defined(BGQ)
      EXTERNAL MPI_WTIME
      REAL*8 MPI_WTIME
      util_cpusec = MPI_WTIME()
#elif defined(IBM) || defined(SP1) ||defined(HPUX)
      double precision ibm_cputime
      external ibm_cputime
      util_cpusec = ibm_cputime()
#elif defined(SUN) || defined(SOLARIS) || defined(LINUXALPHA)
      real*4 tarray(2)
      real*4 etime
#ifndef GCC4
      external etime
#endif
      util_cpusec = etime(tarray)
#elif (defined(LINUX) || defined(CYGNUS)) && !defined(LINUXIA64) && !defined(CATAMOUNT)
      double precision linux_cputime
      external linux_cputime
      util_cpusec = linux_cputime()
#elif defined(WIN32) &&!defined(__MINGW32__)
      double precision win32_cputime
      external win32_cputime
      util_cpusec = win32_cputime()
#elif defined(CRAY)
      double precision secondr
      external secondr
      double precision first
      logical ofirst
      save first, ofirst
      data ofirst /.true./
      if (ofirst) then
        first = secondr()
        ofirst = .false.
      endif
      util_cpusec = secondr() - first ! Actually wall time
#elif defined(FUJITSU_SOLARIS) || defined(CATAMOUNT)
#ifdef FUJITSU_SOLARIS
      intrinsic cpu_time
#endif
      call cpu_time (util_cpusec)
#else
#include "tcgmsg.fh"
      util_cpusec = TCGTIME()
#endif
      end
