      logical function stpr_gen_hess_at(rtdb,delta)
      implicit none
*     
*     $Id: stpr_gh_at.F,v 1.10 1997-06-12 22:52:25 d3e129 Exp $
*     routine to generate hessian at current geometry using default basis set
*     ra kendall 7/95 PNL
*     
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
c::   functions
      logical stpr_check_genat_restart
      external stpr_check_genat_restart
c::   passed
      integer rtdb
      double precision delta
c::   local
      integer whoiam, master
      integer geom
      integer nat, rank_hess
      integer h_gradm, k_gradm  ! handle/index for central gradient
      integer h_gradp, k_gradp  ! handle/index for delta gradient
      integer h_hess, k_hess    ! handle/index for hessian matrix
      integer h_ddipole, k_ddipole
      integer iatom, ixyz
      integer iatom_start, ixyz_start
      integer h_active, i_active, iat_test, nat_active
      logical all_atom, skip_atom
      logical ignore_status
      logical restart
      logical xyz_pass
      logical no_x, no_y, no_z
      logical dipole_okay
      double precision xyz(3), chg
      double precision pdipole(3), mdipole(3)
      character*255 tmpname
      character*16 tag_name
      character*40 new_geom_name
      logical task_gradient, task_save_state, task_restore_state,
     $     task_delete_state, task_energy
      external task_gradient, task_save_state, task_restore_state,
     $     task_delete_state, task_energy
c     
      dipole_okay = .true.
      stpr_gen_hess_at = .false.
      whoiam = ga_nodeid()
      master = 0
      ignore_status = rtdb_parallel(.false.)
c
      call stpr_filenames       ! Initialize common with all stepper filename
c
      if (whoiam.eq.master) then
c     
c     create/load reference geometry
c     
         if (.not.geom_create(geom,'geometry')) call errquit
     $        ('stpr_gen_hess_at:geom_create failed?',1)
         if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     $        ('stpr_gen_hess_at:geom_rtdb_load failed?',2)
c     
c     get the number of atoms
c     
         if (.not. geom_ncent(geom,nat)) call errquit
     $        ('stpr_gen_hess_at:geom_ncent failed?',3)
c     
c     copy reference geometry and store on rtdb
c     
         if (.not.
     $        geom_rtdb_store(rtdb,geom,'reference'))
     $        call errquit
     $        ('stpr_gen_hess_at: geom_rtdb_store failed',4)
         if (.not.geom_destroy(geom))
     $        call errquit
     $        ('stpr_gen_hess_at: geom_destroy failed?',5)
      endif
      call ga_sync()
c     
c     all nodes working here
c     
      ignore_status = rtdb_parallel(.true.)
      if (.not.rtdb_ma_get(rtdb,'gen_hess:active atoms',MT_INT,
     $     nat_active,h_active)) then
         all_atom = .true.
      else
         all_atom = .false.
         if (.not.MA_get_index(h_active,i_active)) call errquit
     $        ('stpr_gen_hess_at: ma_get_index error',6)
      endif
c     set no xyz options
      if (.not.rtdb_get(rtdb,'vib:no_x',mt_log,1,no_x))
     $     no_x = .false.
      if (.not.rtdb_get(rtdb,'vib:no_y',mt_log,1,no_y))
     $     no_y = .false.
      if (.not.rtdb_get(rtdb,'vib:no_z',mt_log,1,no_z))
     $     no_z = .false.
      if ( whoiam .eq. master ) then
        if (util_print('directions', print_medium)) then
          if (no_x) write(luout,*)' stpr_gen_hess_at: no X steps taken'
          if (no_y) write(luout,*)' stpr_gen_hess_at: no Y steps taken'
          if (no_z) write(luout,*)' stpr_gen_hess_at: no Z steps taken'
        endif
      endif
c     all nodes work here 
      restart = .false.
      restart = stpr_check_genat_restart(rtdb,iatom_start,ixyz_start)
c
      if (.not. restart) then
c
c     Perform an energy calculation at the central point and
c     save wavefunction information for rapid convergence
c
         if (.not. task_energy(rtdb)) call errquit
     $        ('stpr_gen_hess_at: reference point energy failed',7)
         if (.not. task_save_state(rtdb,'hess')) call errquit
     $        ('stpr_gen_hess_at: failed to save central point',8)
c
c     save reference geometry
c     
         if (whoiam .eq. master) then
            ignore_status = rtdb_parallel(.false.)
            if (.not.geom_create(geom,'geometry')) call errquit
     $           ('stpr_gen_hess_at:geom_create failed?',9)
            if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     $           ('stpr_gen_hess_at:geom_rtdb_load failed?',10)
            if (.not.geom_rtdb_store(rtdb,geom,'reference'))call errquit
     $           ('stpr_gen_hess_at: geom_rtdb_store failed',11)
            if (.not.geom_destroy(geom)) call errquit
     $           ('stpr_gen_hess_at: geom_destroy failed?',12)
            ignore_status = rtdb_parallel(.true.)
         endif
         call ga_sync()
      endif
c     
c     broadcast number of atoms
c     
      call ga_brdcst(Msg_gen_at_nat+MSGINT, nat, mitob(1), 0)
c     
      call ga_sync()
c     
      if (whoiam.eq.master) then
         rank_hess = 3*nat
         if (.not. MA_Push_Get(MT_DBL,rank_hess,
     $        'minus delta nuclear gradient vector',
     $        h_gradm,k_gradm)) call errquit
     $        ('stpr_gen_hess_at: ma_push_get failed?',13)
         if (.not. MA_Push_Get(MT_DBL,rank_hess,
     $        'plus delta nuclear gradient vector',
     $        h_gradp,k_gradp)) call errquit
     $        ('stpr_gen_hess_at: ma_push_get failed?',14)
         if (.not. MA_Push_Get(MT_DBL,3*rank_hess,
     &       'derivative dipole moment',
     &       h_ddipole, k_ddipole)) call errquit
     &       ('stpr_gen_hess_at: ma_push_get failed',15)
         if (.not. MA_Push_Get(MT_DBL,(rank_hess*rank_hess),
     $        'nuclear hessian matrix',
     $        h_hess,k_hess)) call errquit
     $        ('stpr_gen_hess_at: ma_push_get failed?',16)
c
         if (restart) then 
            call stpr_get_genat_restart(rank_hess,
     $           dbl_mb(k_hess),
     $           dbl_mb(k_gradm),.false.)
* no restart for dipole yet
            dipole_okay = .false.
            call dfill((3*rank_hess),0.0d00,dbl_mb(k_ddipole),1)
         else
            call dfill((rank_hess*rank_hess),0.0d00,
     $           dbl_mb(k_hess),1)
            call dfill((3*rank_hess),0.0d00,dbl_mb(k_ddipole),1)
         endif
      endif
c     
c     all nodes doing work here
c     
      call ga_sync()
      if(.not.restart) then
         iatom_start = 1
         xyz_pass    = .true.
      else
         xyz_pass = .false.
         if (whoiam.eq.master) then
           if (util_print('restart', print_low)) then
             write(luout,*)' **** gen_hess restart ****'
             write(luout,*)' iatom_start = ',iatom_start
             write(luout,*)' ixyz_start  = ',ixyz_start
             call util_flush(luout)
           endif
         endif
      endif
      do 00100 iatom = iatom_start,nat
         if (.not.all_atom) then
            skip_atom = .true.
            do iat_test = 0,(nat_active-1)
               if (iatom.eq.int_mb(i_active+iat_test))
     $              skip_atom = .false.
            enddo
         else
            skip_atom = .false.
         endif
         if (skip_atom) goto 00201
         if (xyz_pass) then
            ixyz_start = 1
         else
            xyz_pass = .true.
         endif
         do 00200 ixyz = ixyz_start,3
            if (ixyz.eq.1 .and. no_x) goto 00199
            if (ixyz.eq.2 .and. no_y) goto 00199
            if (ixyz.eq.3 .and. no_z) goto 00199
            ignore_status = rtdb_parallel(.false.)
            if(whoiam.eq.master) then
              if (util_print('directions', print_medium)) then
                write(luout,*)' iatom = ',iatom, ' ixyz = ',
     $               ixyz,' plus '
                call util_flush(luout)
              endif
c     
c     master node only
c     
               if (.not.geom_create(geom,'reference')) call errquit
     $              ('stpr_gen_hess_at:geom_create failed?',17)
               if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_rtdb_load failed?',18)
               if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_cent_get failed?',19)
               xyz(ixyz) = xyz(ixyz) + delta
               if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_cent_get failed?',20)
c     
c     copy modified geometry and store on rtdb
c     
               write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'plus'
               new_geom_name = 'reference '//tmpname(1:11)
               if (.not.
     $              geom_rtdb_store(rtdb,geom,new_geom_name))
     $              call errquit
     $              ('stpr_gen_hess_at: geom_rtdb_store failed',21)
               if (.not.geom_destroy(geom))
     $              call errquit
     $              ('stpr_gen_hess_at: geom_destroy failed?',22)
c     
c     set "geometry" to "reference atom/xyz"
c     
               if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     $              call errquit
     $              ('stpr_gen_hess_at: rtdb_cput failed ?',23)
            endif
            ignore_status = rtdb_parallel(.true.)
            call ga_sync()
c
c     restore starting information for central point
c     
            if (.not. task_restore_state(rtdb,'hess')) call errquit
     $          ('stpr_gen_hess_at: failed to restore central point',
     &          24)
c     
            if (.not. task_gradient(rtdb)) call errquit
     $           ('stpr_gen_hess_at: energy+grad failed',25)
c     
            call ga_sync()
            ignore_status = rtdb_parallel(.false.)
            if (whoiam.eq.master) then
               if (.not. geom_rtdb_delete(rtdb,new_geom_name))
     $              call errquit
     $              ('stpr_gen_hess_at: geom_rtdb_delete failed',26)
               if(.not. rtdb_get(rtdb,'task:gradient',MT_DBL,
     $              rank_hess,dbl_mb(k_gradp)))
     $              call errquit('stpr_gen_hess_at: rtdb_get failed',27)
               if(.not. rtdb_get(rtdb,'task:dipole',mt_dbl,
     &             3,pdipole)) dipole_okay = .false.
c     
               if (util_print('directions', print_medium)) then
                 write(luout,*)
     $                ' iatom = ',iatom, ' ixyz = ',ixyz,' minus'
                 call util_flush(luout)
               endif
c
c     master node only
c     
               if (.not.geom_create(geom,'reference')) call errquit
     $              ('stpr_gen_hess_at:geom_create failed?',28)
               if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_rtdb_load failed?',29)
               if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_cent_get failed?',30)
               xyz(ixyz) = xyz(ixyz) - delta
               if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     $              call errquit
     $              ('stpr_gen_hess_at:geom_cent_get failed?',31)
c     
c     copy modified geometry and store on rtdb
c     
               write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'minus'
               new_geom_name = 'reference '//tmpname(1:11)
               if (.not.
     $              geom_rtdb_store(rtdb,geom,new_geom_name))
     $              call errquit
     $              ('stpr_gen_hess_at: geom_rtdb_store failed',32)
               if (.not.geom_destroy(geom))
     $              call errquit
     $              ('stpr_gen_hess_at: geom_destroy failed?',33)
c     
c     set "geometry" to "reference atom/xyz"
c     
               if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     $              call errquit
     $              ('stpr_gen_hess_at: rtdb_cput failed ?',34)
            endif
            ignore_status = rtdb_parallel(.true.)
            call ga_sync()
c
c     restore starting information from central point
c     
            if (.not. task_restore_state(rtdb,'hess')) call errquit
     $          ('stpr_gen_hess_at: failed to restore central point',
     &          35)
c     
            if (.not. task_gradient(rtdb))
     $          call errquit('stpr_gen_hess_at: energy+grad failed',
     &          36)
c     
            call ga_sync()
            ignore_status = rtdb_parallel(.false.)
            if (whoiam.eq.master) then
               if (.not. geom_rtdb_delete(rtdb,new_geom_name))
     $            call errquit
     $            ('stpr_gen_hess_at: geom_rtdb_delete (2) failed',37)
               if(.not. rtdb_get(rtdb,'task:gradient',MT_DBL,
     $             rank_hess,dbl_mb(k_gradm)))
     $             call errquit('stpr_gen_hess_at: rtdb_get failed',38)
c
               if(.not. rtdb_get(rtdb,'task:dipole',mt_dbl,
     &             3,mdipole)) dipole_okay = .false.
c     
               call stpr_fd_upd_hess(dbl_mb(k_hess),dbl_mb(k_gradm),
     $              dbl_mb(k_gradp),2.0d00,delta,nat,iatom,ixyz)
               if (dipole_okay)
     &             call stpr_fd_upd_dipole(dbl_mb(k_ddipole),
     &             mdipole,pdipole,2.0d00,delta,nat,iatom,ixyz)
            endif
00199       continue
            if (whoiam.eq.master) then
               call stpr_put_genat_restart(rank_hess,dbl_mb(k_hess),
     $              dbl_mb(k_gradm),iatom,ixyz,nat,.false.)
               ignore_status = rtdb_parallel(.true.)
            endif
00200    continue
00201    continue
00100 continue
c     
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
         if (all_atom) then
c     
c     average contributions from finite diff hess since they "should"
c     be the same
c     
            call stpr_gen_hess_foldave(dbl_mb(k_hess),rank_hess)
         else
            call stpr_gen_hess_fold(dbl_mb(k_hess),rank_hess)
c     
c     assumes zeroed hessian prior to fd computations
c     
            call stpr_gen_set_diag(dbl_mb(k_hess),rank_hess)  
         endif
         if (util_print('hessian', print_low)) then
           write(luout,*)' '
           write(luout,*)' finite difference hessian delta = ',delta
           call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     $          rank_hess,rank_hess,1)
           write(luout,*)' '
         endif
         if (dipole_okay.and.
     &       util_print('derivative dipole',print_low)) then
           write(tmpname,*)
     &         ' finite difference derivative dipole; delta = ',delta
           call stpr_print_ddipole(dbl_mb(k_ddipole),
     &         tmpname, nat)
         endif
         call stpr_wrt_fd_from_sq(dbl_mb(k_hess),rank_hess,FILEHESS)
         if (dipole_okay)
     &       call stpr_wrt_fd_dipole(dbl_mb(k_ddipole),nat,FILEDDIPOLE)
         if (util_print('inforamation', print_medium)) then
           write(luout,*)' triangle hessian written to ',
     $          FILEHESS(1:inp_strlen(FILEHESS))
           if (dipole_okay)
     &         write(luout,*)' derivative dipole written to ',
     &         FILEDDIPOLE(1:inp_strlen(FILEDDIPOLE))
         endif
c     
c     set "geometry" to "geometry"
c     
         if(.not.rtdb_cput(rtdb,'geometry',1,'geometry'))
     $        call errquit
     $        ('stpr_gen_hess_at: rtdb_cput failed ?',39)
c
c     restore and delete save state
c
         if (.not. task_restore_state(rtdb, 'hess')) call errquit
     $        ('stpr_gen_hess_at: failed to restore central point',40)
         if (.not. task_delete_state(rtdb, 'hess')) call errquit
     $        ('stpr_gen_hess_at: failed to restore central point',41)
c     
c     free memory
c     
         
         if(.not.ma_pop_stack(h_hess)) call errquit
     $        ('stpr_gen_hess_at: ma_pop_stack(h_hess) failed?',42)
         if(.not.ma_pop_stack(h_ddipole)) call errquit
     $        ('stpr_gen_hess_at: ma_pop_stack(h_ddipole) failed?',43)
         if(.not.ma_pop_stack(h_gradp)) call errquit
     $        ('stpr_gen_hess_at: ma_pop_stack(h_gradp) failed?',44)
         if(.not.ma_pop_stack(h_gradm)) call errquit
     $        ('stpr_gen_hess_at: ma_pop_stack(h_gradm) failed?',45)
      endif
      if(.not.all_atom) then
         if (.not.ma_free_heap(h_active)) call errquit
     $        ('stpr_gen_hess_at: ma_free_heap(h_active) failed?',46)
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
      stpr_gen_hess_at = .true.
      end
