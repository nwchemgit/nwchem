      logical function cphf(rtdb)
*
* $Id$
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
c
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
c
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
c
c       if (ga_nodeid().eq.0) then
c         write(6,*) "omega: ",omega
c         write(6,*) "lifetime: ",lifetime
c         write(6,*) "gamwidth: ",gamwidth
c       end if
c
      call cphf_solve3(rtdb, omega, lifetime, gamwidth)
c
      cphf3=.true.
c
      end

      logical function cphf4(rtdb)
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
         call util_print_centered(6,'NWChem CPHF4 Module',40,.true.)
         write(6,*)
         write(6,*) 
      endif
c
      call cphf_solve4(rtdb)
c
      cphf4=.true.
c
      end
