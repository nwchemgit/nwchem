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
      if(.not.context_push('cphf'))
     $     call errquit('cphf: context_push failed',0)
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(ga_nodeid().eq.0)
     $     call util_print_centered(6, 'NWCHEM CPHF Module', 40, .true.)
      write(6,*)
      write(6,*)
c
      call cphf_poliz(rtdb)
c     
      if(.not.context_pop('cphf'))
     $     call errquit('context_pop failed',0)
      end


