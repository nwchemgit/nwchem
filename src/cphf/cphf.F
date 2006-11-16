      logical function cphf(rtdb)
*
* $Id: cphf.F,v 1.8 2006-11-16 18:22:50 jochen Exp $
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


c ... jochen: now I have added even a third version of cphf to
c     deal with frequency dependent properties

      logical function cphf3(rtdb, omega, lifetime, gamwidth)
*
* This is temporary until the cphf is reunited under one umbrella
*
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c     
      integer rtdb              ! [input] database handle
      double precision omega    ! [input] frequency
      logical lifetime          ! [input] dampinfg switch
      double precision gamwidth ! [input] damping parameter

      logical status
c     
      status=rtdb_parallel(.true.) ! Broadcast reads to all processes
c     
      if(util_print('information',print_default) .and.
     $     ga_nodeid().eq.0) then
       
        if (ga_nodeid().eq.0) write (luout,*)
         call util_print_centered(6,'NWChem Dynamic CPHF Module',
     &     40,.true.)
         write(6,*)
         write(6,*) 
      endif
      call cphf_solve3(rtdb, omega, lifetime, gamwidth)
c
      cphf3=.true.
c
      end
