      logical function cphf(rtdb)
*
* $Id: cphf.F,v 1.7 2001-11-09 23:09:19 windus Exp $
*
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb              ! [input] database handle
      logical status
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(util_print('information',print_default) .and.
     $     ga_nodeid().eq.0) then
         call util_print_centered(6,'NWChem CPHF Module',40,.true.)
         write(6,*)
         write(6,*) 
      endif
      call cphf_solve(rtdb)
c
      cphf=.true.
c
      end
      logical function cphf2(rtdb)
*
* This is temporary until the cphf is reunited under one umbrella
*
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb              ! [input] database handle
      logical status
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(util_print('information',print_default) .and.
     $     ga_nodeid().eq.0) then
         call util_print_centered(6,'NWChem CPHF Module',40,.true.)
         write(6,*)
         write(6,*) 
      endif
      call cphf_solve2(rtdb)
c
      cphf2=.true.
c
      end
