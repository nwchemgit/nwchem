c
c     $Id$
c
      double precision function perfm_flop()
      implicit none
#include "perfm.fh"
#ifdef USE_PAPI
      integer ii
      integer*4 arg1
      if(lfailed) then
	perfm_flop=0d0
        return
      endif
#ifdef PAPI_FPINS
      call PAPIf_flips( real_time, proc_time, flpins, mflops, retval )
      perfm_flop=flpins*conv2m
!      if (retval .NE. PAPI_OK) call 
!     $ perfm_error(retval,
!     $ ' perfm_start: could not init papi flpins ')
#else
c
c     addtional papi events
c
#endif
#else
      perfm_flop=0d0
#endif
      return
      end
