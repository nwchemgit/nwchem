      subroutine util_setmaxthreads()
c
c     sets the number of openmp threads as
c     number of hw cores / no. of process per node
c
      implicit none
#include "global.fh"
#include "stdio.fh"
      integer ncpus,util_getncpus,max_thr,ppn
      character*255 env_out
      integer err8
c
c     skip if OMP_NUM_THREADS is set
c
      call get_environment_variable(
     A     "OMP_NUM_THREADS", value=env_out, status=err8)
      if(err8.ne.0) then
         ncpus=util_getncpus()
         if(ncpus.ne.0) then
            call util_getppn(ppn)
            max_thr=ncpus/ppn
            if(ga_nodeid().eq.0)write(luout,97)
     N           ncpus,max_thr
            call omp_set_num_threads(max_thr)
         endif
      endif
 97   format(2x,'NWChem w/ OpenMP: # of HW cores   = ',i4,
     F     /2x,'NWChem w/ OpenMP: nprocs/node     = ',i4)
      return
      end
