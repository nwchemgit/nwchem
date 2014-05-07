      double precision function util_wallsec()
      implicit none
c $Id$
*
* routine to give wall clock seconds since start of execution
*  delta times are handled by all applications
*
*
      logical first_call
      double precision first_time,last_time,cur_time,diff
      double precision MPI_Wtime
      external MPI_Wtime
c
      save first_time
      save first_call
      data first_call /.true./
c
      util_wallsec = MPI_Wtime()
      if (first_call) then
         first_time = MPI_Wtime()
         first_call = .false.
      endif

      cur_time = MPI_Wtime()
      util_wallsec = cur_time - first_time

*
      end


