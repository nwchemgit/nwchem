      double precision function util_cpusec()
      implicit none
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
      end
