      logical function stpr_walk(rtdb)
c     $Id: stpr_walk.F,v 1.21 1996-07-26 16:39:34 d3g681 Exp $
      implicit none
c     
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "pstat.fh"
#include "util.fh"
#include "stdio.fh"
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
      character*16 stpr_tags(max_cent) ! local tags array
      integer i_timewalk
      logical o_timewalk
      logical status
      double precision convgg, convge
#include "cstprfiles.fh"
      integer first_call, step_number
      data first_call /0/
c     
      call util_print_push
      call util_print_rtdb_load(rtdb,'stepper')
c     
      call stpr_filenames()
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
     &        call errquit('stpr_walk: pstat allocate failed',911)
         call pstat_on(i_timewalk)
      endif
c     
      stpr_walk = .false.
c....................................set broadcast off
      status = rtdb_parallel(.false.) 
c     
*     rak:      status = rtdb_print(rtdb,.true.)
      if (ga_nodeid().eq.0)  then
         if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $        title = ' '
         if (util_print('information',print_low)) then
            call util_print_centered(6, 'NWCHEM STEPPER Module', 
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
     &        ('stpr_walk: failed to read energy from rtdb', 911)
c     
         if (.not. geom_create(geom, 'geometry'))
     &        call errquit('stpr_walk: geom_create?', 911)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &        call errquit('stpr_walk: no geometry ', 911)
c     
*     rak:        write(6,*)' after load in stpr_walk'
*     rak:        if (.not.geom_print(geom)) stop 'tmp error'
c     
         if (.not. geom_ncent(geom,nat))
     &        call errquit('stpr_walk: geom_ncent?',911)
c     
         call stpr_cneed(needed,nat) ! determine core needed
         needed = (needed + 1023)*1024/1024 ! nearest 1024 block
c     
c.....space for core array
         if (.not.MA_Push_Get(MT_DBL,needed,'stpr core',h_core,k_core))
     &        call errquit('stpr_walk: allocation for core failed?',911)
c.....space for coords
         if (.not.
     &        MA_Push_Get(MT_DBL,(3*nat),'stpr coords',h_coord,k_coord))
     &        call errquit
     &        ('stpr_walk: allocation for coords failed?',911)
c.....space for gradient
         if (.not.
     &        MA_Push_Get(MT_DBL,(3*nat),'stpr gradient',h_grad,k_grad))
     &        call errquit
     &        ('stpr_walk: allocation for gradient failed?',911)
c.....space for charge
         if (.not.
     &        MA_Push_Get(MT_DBL,nat,'stpr geom charge',h_chg,k_chg))
     &        call errquit
     &        ('stpr_walk: allocation for geom charge failed?',911)
c     
c.....read gradients
         if (.not. rtdb_get(rtdb, 'task:gradient', MT_DBL, 
     &        (3*nat),dbl_mb(k_grad)))
     &        call errquit
     &        ('stpr_walk: reading gradients failed',911)
c     
c.....get coordinates
         if (.not. geom_cart_get(geom,nat,
     &        stpr_tags,
     &        dbl_mb(k_coord),
     &        dbl_mb(k_chg)))
     &        call errquit
     &        ('stpr_walk: geom_cart_get failed ?',911)
c     
c...  take a step
c     
         convgg = 1.0d-04
         convge = 1.0d-07
         call stpr_stepper(
     &        dbl_mb(k_core),
     &        needed,
     &        dbl_mb(k_grad),
     &        dbl_mb(k_coord),
     &        energy,nat,convge,convgg,stpr_walk,
     &        rtdb,step_number)
c.....set coordinates
         if (.not. geom_cart_set(geom,nat,
     &        stpr_tags,
     &        dbl_mb(k_coord),
     &        dbl_mb(k_chg)))
     &        call errquit
     &        ('stpr_walk: geom_cart_set failed ?',911)
c     
c.....copy geom on rtdb with new name
c     
*     tmpname = ' '
*     write(tmpname,'(i4)')abs(step_number)
*     new_geom_name='g-'//tmpname//'-step'
         new_geom_name = ' '
         if (abs(step_number) .lt. 10) then
            write(new_geom_name,'(''g-'',i1,''-step'')') 
     $           abs(step_number)
         else if (abs(step_number) .lt. 100) then
            write(new_geom_name,'(''g-'',i2,''-step'')') 
     $           abs(step_number)
         else if (abs(step_number) .lt. 1000) then
            write(new_geom_name,'(''g-'',i3,''-step'')') 
     $           abs(step_number)
         else
            write(new_geom_name,'(''g-'',i4,''-step'')') 
     $           abs(step_number)
         endif
c     
         call sym_geom_project(geom, 1d-4)
c     

         if (.not.geom_rtdb_store(rtdb,geom,new_geom_name))
     &        call errquit
     &        ('stpr_walk: geom_rtdb_store (of copy) failed',911)
c     
C**** #define ECCE
#if defined(ECCE)
         write(luout,*)'<stepper updated geometry> start '
         if(.not. geom_print(geom)) call errquit
     &        ('stpr_walk:geom print failed ',911)
         write(luout,*)'<stepper updated geometry> end'
#endif
c     
         if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &        call errquit
     &        ('stpr_walk: rtdb_cput failed ?',911)
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
     &        call errquit('stpr_walk: geom_destroy failed?',911)
*     if (.not.geom_destroy(new_geom))
*     &      call errquit('stpr_walk: geom_destroy failed?',911)
         if (.not.MA_Pop_Stack(h_chg))
     &        call errquit('stpr_walk: ma_pop_stack h_chg failed',911)
         if (.not.MA_Pop_Stack(h_grad))
     &        call errquit('stpr_walk: ma_pop_stack h_grad failed',911)
         if (.not.MA_Pop_Stack(h_coord))
     &        call errquit('stpr_walk: ma_pop_stack h_coord failed',911)
         if (.not.MA_Pop_Stack(h_core))
     &        call errquit('stpr_walk: ma_pop_stack h_core failed',911)
c     
         int_val = 0 
         if (stpr_walk) int_val = 1
         if(step_number.le.0) then
            write(6,*)' maximum number of geometry optimization ',
     &           'iterations exceeded:  graceful shutdown will occur'
            int_val = -1
         endif
      endif
c     
      call ga_brdcst(Msg_Stpr_Walk+MSGINT, int_val, mitob(1), 0) ! Propagate walker status
c     
      stpr_walk = .false.
      if (int_val.ne.0) stpr_walk = .true.
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
     &        call errquit('stpr_walk: pstat_free failed?',911)
      endif
c     
      call util_print_pop
      call ga_sync()
c     
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

