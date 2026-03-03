      subroutine util_setmaxthreads()
      implicit none
#include "global.fh"
#include "stdio.fh"
      integer ncpus,util_getncpus,max_thr,ppn
      external util_getncpus
c
      ncpus=util_getncpus()
      if(ncpus.ne.0) then
         call util_getppn(ppn)
         max_thr=ncpus/ppn
         if(ga_nodeid().eq.0)write(luout,*) ' ### ncpus max_thr',
     N        ncpus,max_thr
         call omp_set_num_threads(max_thr)
      endif
      return
      end
