      double precision function util_wallsec()
      implicit none
#ifdef NO_MPIF
#include "tcgmsg.fh"
#endif
c $Id$
*
* routine to give wall clock seconds since start of execution
*  delta times are handled by all applications
*
*
      logical first_call
      double precision first_time,last_time,cur_time,diff
#ifndef NO_MPIF
      double precision MPI_Wtime
      external MPI_Wtime
#endif
c
      save first_time
      save first_call
      data first_call /.true./
c
#ifdef NO_MPIF
      util_wallsec = TCGTIME()
#else
      util_wallsec = MPI_Wtime()
#endif
      if (first_call) then
#ifdef NO_MPIF
         first_time = TCGTIME()
#else
         first_time = MPI_Wtime()
#endif
         first_call = .false.
      endif

#ifdef NO_MPIF
      cur_time = TCGTIME()
#else
      cur_time = MPI_Wtime()
#endif
      util_wallsec = cur_time - first_time

*
      end


