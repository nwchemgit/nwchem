C> \ingroup stpr_api
C> @{
C>
      integer function stpr_walk(rtdb)
c     $Id$
      implicit none
#include "errquit.fh"
c     
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "nwc_const.fh"
********************************************************************#include "geomP.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "pstat.fh"
#include "util.fh"
#include "stdio.fh"
c
c     
      integer rtdb              ! [input] run-time-data-base handle
c     
c     parameters read from rtdb
c     prefix:energy (later) rohf and uhf
c     
      integer geom              ! geometry handle
      integer needed            ! core needed
      integer nat               ! number of atoms 
      character*255 title       ! title buffer
      double precision energy   ! energy
c     
      integer h_grad, k_grad    ! MA handle/index of gradient
      integer h_coord, k_coord  ! MA handle/index of coords
      integer h_core, k_core    ! MA handle/index of scratch array for stepper
      integer h_chg, k_chg      ! MA handle/index of charge array
      integer int_val           ! value to propigate status (stpr_walk value)
      character*40 new_geom_name ! new geometry name for next step
      character*16 stpr_tags(nw_max_atom) ! local tags array
      integer i_timewalk
      logical o_timewalk, okeepgeom
      logical status
      logical lstpr_walk
      double precision convggm, convgg, convge
      integer mitob1
#include "cstprfiles.fh"
      integer first_call, step_number, astep_number
      data first_call /0/
c     
      call util_print_push
      call util_print_rtdb_load(rtdb,'stepper')
c     
      call stpr_filenames()
c     
      call ecce_print_module_entry ('stepper')
c
      if (first_call.eq.0) then
         first_call = 1
      endif
c     
#if defined(KSR)||defined(SUN)||defined(__PARAGON__)
      o_timewalk = util_print('time',print_default)
#else 
      o_timewalk = .false.
#endif
c     
      if (o_timewalk) then
         if (.not.(pstat_allocate
     &        ('stepper time',Pstat_Qstat,0,0,i_timewalk)))
     &        call errquit('stpr_walk: pstat allocate failed',911,
     &       MEM_ERR)
         call pstat_on(i_timewalk)
      endif
c     
      lstpr_walk = .false.
c....................................set broadcast off
      status = rtdb_parallel(.false.) 
c     
*     rak:      status = rtdb_print(rtdb,.true.)
      if (ga_nodeid().eq.0)  then
         if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $        title = ' '
         if (util_print('information',print_low)) then
            call util_print_centered(6, 'NWChem STEPPER Module', 
     $           40, .true.)
            write(6,*)
            write(6,*)
            if (title .ne. ' ') then
               call util_print_centered(6, title, 40, .false.)
               write(6,*)
               write(6,*)
            endif
         endif
c     
         if (.not. rtdb_get(rtdb, 'task:energy', MT_DBL, 1, energy))
     &        call errquit
     &        ('stpr_walk: failed to read energy from rtdb', 911,
     &       RTDB_ERR)
c     
         if (.not. geom_create(geom, 'geometry'))
     &        call errquit('stpr_walk: geom_create?', 911, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &        call errquit('stpr_walk: no geometry ', 911, RTDB_ERR)
c     
*     rak:        write(6,*)' after load in stpr_walk'
*     rak:        if (.not.geom_print(geom)) stop 'tmp error'
c     
         if (.not. geom_ncent(geom,nat))
     &        call errquit('stpr_walk: geom_ncent?',911, GEOM_ERR)
c     
         call stpr_cneed(needed,nat) ! determine core needed
         needed = (needed + 1023)*1024/1024 ! nearest 1024 block
c     
c.....space for core array
         if (.not.MA_Push_Get(MT_DBL,needed,'stpr core',h_core,k_core))
     &        call errquit('stpr_walk: allocation for core failed?',911,
     &       MA_ERR)
c.....space for coords
         if (.not.
     &        MA_Push_Get(MT_DBL,(3*nat),'stpr coords',h_coord,k_coord))
     &        call errquit
     &        ('stpr_walk: allocation for coords failed?',911, MA_ERR)
c.....space for gradient
         if (.not.
     &        MA_Push_Get(MT_DBL,(3*nat),'stpr gradient',h_grad,k_grad))
     &        call errquit
     &        ('stpr_walk: allocation for gradient failed?',911, MA_ERR)
c.....space for charge
         if (.not.
     &        MA_Push_Get(MT_DBL,nat,'stpr geom charge',h_chg,k_chg))
     &        call errquit
     &        ('stpr_walk: allocation for geom charge failed?',911,
     &       MA_ERR)
c     
c.....read gradients
         if (.not. rtdb_get(rtdb, 'task:gradient', MT_DBL, 
     &        (3*nat),dbl_mb(k_grad)))
     &        call errquit
     &        ('stpr_walk: reading gradients failed',911, RTDB_ERR)
c     
c.....get coordinates
         if (.not. geom_cart_get(geom,nat,
     &        stpr_tags,
     &        dbl_mb(k_coord),
     &        dbl_mb(k_chg)))
     &        call errquit
     &        ('stpr_walk: geom_cart_get failed ?',911, GEOM_ERR)
c     
c...  take a step
c     
         convggm= 8.0d-04
         convgg = 1.0d-02
         convge = 1.0d-04
         call stpr_stepper(
     &        dbl_mb(k_core),
     &        needed,
     &        dbl_mb(k_grad),
     &        dbl_mb(k_coord),
     &        energy,nat,convge,convgg,convggm,lstpr_walk,
     &        rtdb,step_number)
c.....set coordinates
         if (.not. geom_cart_set(geom,nat,
     &        stpr_tags,
     &        dbl_mb(k_coord),
     &        dbl_mb(k_chg)))
     &        call errquit
     &        ('stpr_walk: geom_cart_set failed ?',911, GEOM_ERR)
C        output to ECCE
         astep_number = abs(step_number)
         call ecce_print1 ('iteration', mt_int, astep_number, 1)
         call ecce_print2 ('geometry', mt_dbl, dbl_mb(k_coord), 3,3,nat)
c
c     optionally delete the old geometry (done by default)
c     
         if (.not. rtdb_get(rtdb,'stepper:keepgeom', mt_log, 1, 
     $        okeepgeom)) okeepgeom = .false.
         if (step_number.gt.1 .and. (.not. okeepgeom)) then
            call stpr_new_geom_name(new_geom_name, step_number-1)
            if (.not. geom_rtdb_delete(rtdb,new_geom_name))
     $           call errquit('stpr_walk: rtdb corrupt? geom delete',
     $           step_number-1, RTDB_ERR)
         endif
c
c     copy the new geometry to the database with a new name
c
         call stpr_new_geom_name(new_geom_name, step_number)
c     
         call sym_geom_project(geom, 1d-4)
c     
         if (.not.geom_rtdb_store(rtdb,geom,new_geom_name))
     &        call errquit
     &        ('stpr_walk: geom_rtdb_store (of copy) failed',911,
     &       RTDB_ERR)
         if (util_print('new coordinates',print_low)) then
           if (lstpr_walk) then
             write(luout,11001)
           else
             write(luout,11002)
           endif
           if (.not.geom_print(geom)) call errquit(
     &         'stpr_walk: geom_print failed',911, GEOM_ERR)
         endif
         if (util_print('distances',print_default)) then
           if (.not.geom_print_distances(geom)) call errquit(
     &         'stpr_walk: geom_print_distances failed',911, GEOM_ERR)
         endif
         if (util_print('angles',print_default)) then
           if (.not.geom_print_angles(geom)) call errquit(
     &       'stpr_walk: geom_print_angles failed',911, GEOM_ERR)
         endif
         if (util_print('dihedrals',print_default)) then
           if (.not.geom_print_dihedrals(geom)) call errquit(
     &         'stpr_walk: geom_print_angles failed',911, GEOM_ERR)
         endif
c     
C**** #define ECCE
#if defined(ECCE)
         write(luout,*)'<stepper updated geometry> start '
         if(.not. geom_print(geom)) call errquit
     &        ('stpr_walk:geom print failed ',911, GEOM_ERR)
         write(luout,*)'<stepper updated geometry> end'
#endif
c     
         if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &        call errquit
     &        ('stpr_walk: rtdb_cput failed ?',911, RTDB_ERR)
c     
*     if(.not.rtdb_print(rtdb,.true.)) stop 'tmp error'
c     
*     rak---- old copy over geometry - start
*     rak:c
*     rak:c..... update rtdb with new coordinates        
*     rak:c
*     rak:*        if (.not. geom_print(geom))
*     rak:*     &      call errquit('stpr_walk: geom_print ?',911)
*     rak:        if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
*     rak:     &      call errquit
*     rak:     &      ('stpr_walk: geom_rtdb_store failed ?',911)
*     rak:c
*     rak---- old copy over geometry - end
c...  clean up
         if (.not.geom_destroy(geom))
     &        call errquit('stpr_walk: geom_destroy failed?',911,
     &       GEOM_ERR)
*     if (.not.geom_destroy(new_geom))
*     &      call errquit('stpr_walk: geom_destroy failed?',911)
         if (.not.MA_Pop_Stack(h_chg))
     &        call errquit('stpr_walk: ma_pop_stack h_chg failed',911,
     &       MA_ERR)
         if (.not.MA_Pop_Stack(h_grad))
     &        call errquit('stpr_walk: ma_pop_stack h_grad failed',911,
     &       MA_ERR)
         if (.not.MA_Pop_Stack(h_coord))
     &        call errquit('stpr_walk: ma_pop_stack h_coord failed',911,
     &       MA_ERR)
         if (.not.MA_Pop_Stack(h_core))
     &        call errquit('stpr_walk: ma_pop_stack h_core failed',911,
     &       MA_ERR)
c     
         int_val = 0 
         if (lstpr_walk) int_val = 1
         if(step_number.le.0) then
            write(6,*)' maximum number of geometry optimization ',
     &           'iterations exceeded:  restartable shutdown will occur'
            int_val = -1
         endif
      endif
c
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(Msg_Stpr_Walk+MSGINT, int_val, mitob1, 0) ! Propagate walker status
c     
      lstpr_walk = .false.
      if (int_val.gt.0) lstpr_walk = .true.
      stpr_walk = int_val
c     
c.... restore parallel access to rtdb
      status = rtdb_parallel(.true.) 
c     
      if (o_timewalk) then
         call pstat_off(i_timewalk)
         call pstat_print_handles(i_timewalk,1)
         write(6,*)' '
         o_timewalk = .false.
         if (.not.pstat_free(i_timewalk))
     &        call errquit('stpr_walk: pstat_free failed?',911, MEM_ERR)
      endif
c     
      
      if (stpr_walk.lt.0) then
        call ecce_print_module_exit
     &      ('stepper','failed : maxiter reached')
      else
        call ecce_print_module_exit ('stepper','ok')
      endif
      call util_print_pop
      call ga_sync()
c     
11001 format(/,/,/,'  ',16('-'),' Converged geometry ',16('-'))
11002 format(/,/,/,'  ',19('-'),' New geometry ',19('-'))
      end
      subroutine stpr_walk_reset
      implicit none
c
      integer iwctr, natom, icalc
      COMMON / CFACE / IWCTR,NATOM,ICALC
c
      iwctr = 0
      natom = 0
      icalc = 0
c
      end
C>
C> @}
C> \ingroup stpr_priv
C> @{
C>
      subroutine stpr_new_geom_name(new_geom_name, step_number)
      implicit none
c
      character*(*) new_geom_name
      integer step_number
c
c     Construct a name for the geometry at this step
c
      new_geom_name = ' '
      if (abs(step_number) .lt. 10) then
         write(new_geom_name,'(''g-'',i1,''-step'')') 
     $        abs(step_number)
      else if (abs(step_number) .lt. 100) then
         write(new_geom_name,'(''g-'',i2,''-step'')') 
     $        abs(step_number)
      else if (abs(step_number) .lt. 1000) then
         write(new_geom_name,'(''g-'',i3,''-step'')') 
     $        abs(step_number)
      else
         write(new_geom_name,'(''g-'',i4,''-step'')') 
     $        abs(step_number)
      endif
c
      end
      logical function stpr_walk_task(rtdb)
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb
*
      double precision gradw_time
      integer grad_time
      logical grad_status
      logical me_write
      integer value_stpr_walk
      integer stpr_walk
      external stpr_walk
      logical task_gradient
      external task_gradient
*
      stpr_walk_task = .false.
      call stpr_walk_reset()
*
      me_write = ga_nodeid().eq.0
00010 continue
      gradw_time = (-1.0d00) * util_wallsec()
      grad_status = task_gradient(rtdb)
      gradw_time = (gradw_time + util_wallsec())*1.2d00 ! add 20 percent
      grad_time = int(gradw_time)
      if (.not. grad_status) then
        if (me_write) then
          write(luout,*)' stpr_walk_task: energy+gradient failed',
     &              ' ... optimization incomplete'
          call util_flush(luout)
        endif
      endif
      value_stpr_walk = stpr_walk(rtdb)
c
c     Disable printing of movecs to ecce after the first step
c     ... reenable at the end
c
      call movecs_ecce_print_off()
c
      if (value_stpr_walk .eq. 0) then
        if (util_test_time_remaining(rtdb,grad_time)) then
          goto 00010
        else
          stpr_walk_task = .false.
          if (me_write) then
            write(luout,*)
     &          ' remaining time is not long enough for next step'
            call util_flush(luout)
          endif
        endif
      else if (value_stpr_walk.eq.1) then
        stpr_walk_task = .true.
      endif
c
      call movecs_ecce_print_on()
c
      end
C>
C> @}
