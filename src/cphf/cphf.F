      logical function cphf(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c     
      integer rtdb              ! [input] database handle
      logical status
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(ga_nodeid().eq.0)then
         call util_print_centered(6,'NWCHEM CPHF Module',40,.true.)
         write(6,*)
         write(6,*) 
      endif
      call cphf_solve(rtdb)
      cphf=.true.
      end


