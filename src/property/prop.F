c***********************************************************************
      logical function property(rtdb)
c***********************************************************************
c
c>>>  Inquire rtdb for properties to be evaluated.
c
*
* $Id$
*
      implicit none
      integer rtdb              ! [input]
      integer nbofile,aimfile,moldenfile,polfromsos
      logical status
      logical ostart,ocontinue,orestart
      integer restr
      logical do_fockbuild
      logical  hnd_property
      external hnd_property
#include "context.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
c     Push context down to prop.
c
      call ga_sync()

      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
*      if(ga_nodeid().eq.0)
*     $     call util_print_centered(LUout,'Properties',40,.true.)
c
c     start ecce property output module
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'prop')
      call ecce_print_module_entry('task property')
c
c     Find out if this is a restart
c
      restr = 0
      call util_get_rtdb_state(rtdb,ostart,ocontinue,orestart)
      if (orestart) restr = 1
      if (restr.gt.0.and. ga_nodeid().eq.0)
     &    write(luout,*) 'Restarting calculation'
c
      call solver_setup(rtdb,restr)  ! restart parameters for solver (ga_lkain())
      if (.not.rtdb_get(rtdb,'prop:do_fockbuild',MT_LOG,1,do_fockbuild))
     D     do_fockbuild=.true.
c
c     hondo property calculation routines 
c
      if(do_fockbuild)  then
         status = hnd_property(rtdb)
      endif
c
c     nbo input file
c
      if (rtdb_get(rtdb,'prop:nbofile',MT_INT,1,nbofile)) then
        if(nbofile.eq.1) call wnbofile(rtdb)
      endif 
      if (rtdb_get(rtdb,'prop:aimfile',MT_INT,1,aimfile)) then
        if(aimfile.eq.1) call waimfile(rtdb)
      endif
      if (rtdb_get(rtdb,'prop:moldenfile',MT_INT,1,moldenfile)) then
        if(moldenfile.eq.1) call prop_moldenfile(rtdb)
      endif
      if (rtdb_get(rtdb,'prop:polfromsos',MT_INT,1,polfromsos)) then
        if(polfromsos.eq.1) call prop_polfromsos(rtdb)
      endif
c
c     finish ecce property output module
c
      call ecce_print_module_exit('task property','ok')
      call util_print_pop
      property = .true.
      return
      end
