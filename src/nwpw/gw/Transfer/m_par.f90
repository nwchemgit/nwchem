!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_par.f90,v $:
! $Revision: 1.1 $
! $Author: prange $
! $Date: 2012/10/16 15:23:25 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module par

      implicit none
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm

      end module
