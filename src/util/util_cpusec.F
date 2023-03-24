      double precision function util_cpusec()
      implicit none
#ifndef NO_WALLTIME
! we implement cpu time as wall time because cpu time
! is often unhelpful, especially when threads are used.
! there were issues with the implementation of cpu time
! as well, that motivated deleting the code.
      double precision util_wallsec
      external util_wallsec
      util_cpusec=util_wallsec()
#else
*
* routine to give cpu seconds since start of execution
*  delta times are handled by all applications
*
#if 1
! this is a Fortran 95 standard intrinsic.
      call cpu_time (util_cpusec)
#elif (defined(LINUX) || defined(CYGNUS))
      double precision linux_cputime
      external linux_cputime
      util_cpusec = linux_cputime()
#elif defined(WIN32) && !defined(__MINGW32__)
      double precision win32_cputime
      external win32_cputime
      util_cpusec = win32_cputime()
#else
      EXTERNAL MPI_WTIME
      REAL*8 MPI_WTIME
      util_cpusec = MPI_WTIME()
#endif
#endif
      end
