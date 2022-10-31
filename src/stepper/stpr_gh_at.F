C> \ingroup stpr_api
C> @{
C>
C> \brief Generate the Hessian by finite differences
C>
C> This routine generates the Hessian by finite differences at the
C> current geometry and using the default basis set.
C>
      logical function stpr_gen_hess_at(rtdb,delta)
      implicit none
#include "errquit.fh"
*     
*     $Id$
*     routine to generate hessian at current geometry using default basis set
*     ra kendall 7/95 PNL
*     
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
#include "sym.fh"
c::   functions
      logical numerical_single
      logical stpr_check_genat_restart
      external stpr_check_genat_restart
c::   passed
      integer rtdb           !< [Input] The RTDB handle
      double precision delta !< [Input] Step size of finite displacements
c::   local
      integer i
      integer ma_type
      integer whoiam, master
      integer geom
      integer nat, rank_hess
      integer h_gradm, k_gradm ! handle/index for central gradient
      integer h_gradp, k_gradp ! handle/index for delta gradient
      integer h_hess, k_hess   ! handle/index for hessian matrix
      integer h_hessc,k_hessc
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
      logical o_stepper_debug
      logical o_debug_restart
      logical o_directions
      logical ogeom_aliased
      logical last_pass
      logical oprint
      double precision time_taskg1, time_taskg2
      integer time_taskg, print_level
      double precision xyz(3), chg
      double precision pdipole(3), mdipole(3)
      integer date_len
      character*40 date_string
      character*255 tmpname
      character*16 tag_name
      character*40 new_geom_name
      character*40 first_geom_name
      double precision q1
      logical task_gradient, task_save_state, task_restore_state,
     &    task_delete_state, task_energy,hesscompress
      external task_gradient, task_save_state, task_restore_state,
     &    task_delete_state, task_energy
c     
      oprint = ga_nodeid().eq.0
      o_stepper_debug = util_print('debug_stepper',print_debug)
      o_directions    = util_print('directions', print_medium)
      o_debug_restart = util_print('debug_stepper_restart',print_debug)
      q1=1d0

c     Check if single differencing is requested 
      if (.not. rtdb_get(rtdb,'hessian:numerical_single',mt_log,1,
     &    numerical_single)) then
         numerical_single = .false.
      else
         if (oprint) then
         write (6,*)
         write (6,*) 'Single differencing is requested'
         write (6,*)
         end if
      endif
c
      call util_print_get_level(print_level) ! Save for later
c
      dipole_okay = .true.
      stpr_gen_hess_at = .false.
      whoiam = ga_nodeid()
      master = 0
c
      call stpr_filenames ! Initialize common with all stepper filename
c
      restart = .false.
      restart = stpr_check_genat_restart(rtdb,iatom_start,ixyz_start)
*
      if (o_stepper_debug) then
        call ga_sync()
        write(6,*)' going false 1',whoiam
        call util_flush(6)
      endif
*
      call ga_sync()
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
c     
c     create/load reference geometry
c     
        
        if (restart) then 
          if (.not.geom_create(geom,'reference')) call errquit
     &        ('stpr_gen_hess_at:geom_create failed?',1, GEOM_ERR)
          if (o_debug_restart)
     &        write(6,*)' load 1 reference ',
     &        '   node=',whoiam
          if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &        call errquit
     &        ('stpr_gen_hess_at:geom_rtdb_load failed?',2, RTDB_ERR)
          if (rtdb_cget(rtdb,'stpr:gh:at:alias:name',1,tmpname)) then
            ogeom_aliased = .true.
          else
            ogeom_aliased = .false.
          endif
        else
          if (rtdb_cget(rtdb,'geometry',1,first_geom_name)) then
            ogeom_aliased = .true.
            if (.not.
     &          rtdb_cput
     &          (rtdb,'stpr:gh:at:alias:name',1,first_geom_name)
     &          ) call errquit
     &          ('stpr_gen_hess_at:rtdb_cput of alias name failed?',1,
     &       RTDB_ERR)
          else
            first_geom_name = 'geometry'          
            ogeom_aliased = .false.
          endif
          if (.not.geom_create(geom,first_geom_name)) call errquit
     &        ('stpr_gen_hess_at:geom_create failed?',1, GEOM_ERR)
          if (o_debug_restart)
     &        write(6,*)' load 2 first_geom_name: ',first_geom_name,
     &        '   node=',whoiam
          if (.not.geom_rtdb_load(rtdb,geom,first_geom_name))
     &        call errquit
     &        ('stpr_gen_hess_at:geom_rtdb_load failed?',2, RTDB_ERR)
        endif
c     
c     get the number of atoms
c     
        if (.not. geom_ncent(geom,nat)) call errquit
     &      ('stpr_gen_hess_at:geom_ncent failed?',3, GEOM_ERR)
c     
c     copy reference geometry and store on rtdb
c     
        if (.not.restart) then
          if (.not.
     &        geom_rtdb_store(rtdb,geom,'reference'))
     &        call errquit
     &        ('stpr_gen_hess_at: geom_rtdb_store failed',4, RTDB_ERR)
        endif
        if (.not.geom_destroy(geom))
     &      call errquit
     &      ('stpr_gen_hess_at: geom_destroy failed?',5, GEOM_ERR)
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
*
      if (o_stepper_debug) then
        write(6,*)' gone true 1',whoiam
        call util_flush(6)
        call ga_sync()
      endif
*
c     
c     all nodes working here
c     
c      if (.not.rtdb_ma_get(rtdb,'gen_hess:active atoms',MT_INT,
      if (.not.rtdb_ma_get(rtdb,'gen_hess:actlist',ma_type,
     &    nat_active,h_active)) then
        all_atom = .true.
      else
        all_atom = .false.
        if (.not.MA_get_index(h_active,i_active)) call errquit
     &      ('stpr_gen_hess_at: ma_get_index error',6, MA_ERR)
        if(oprint) then
        write(*,*) "active atoms for hessian generation"
        do i=1,nat_active
             write(*,*) int_mb(i_active+i-1)
        end do
        end if
      endif
c     set no xyz options
      if (.not.rtdb_get(rtdb,'vib:no_x',mt_log,1,no_x))
     &    no_x = .false.
      if (.not.rtdb_get(rtdb,'vib:no_y',mt_log,1,no_y))
     &    no_y = .false.
      if (.not.rtdb_get(rtdb,'vib:no_z',mt_log,1,no_z))
     &    no_z = .false.
*
      if (o_stepper_debug) then
        call ga_sync()
        write(6,*)' going false 2',whoiam
        call util_flush(6)
      endif
*
      call ga_sync()
      ignore_status = rtdb_parallel(.false.)
      if ( whoiam .eq. master ) then
        if ((o_directions).and.(oprint)) then
          if (no_x) write(luout,*)' stpr_gen_hess_at: no X steps taken'
          if (no_y) write(luout,*)' stpr_gen_hess_at: no Y steps taken'
          if (no_z) write(luout,*)' stpr_gen_hess_at: no Z steps taken'
        endif
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
*
      if (o_stepper_debug) then
        write(6,*)' gone true 2',whoiam
        call util_flush(6)
        call ga_sync()
      endif
*
c     all nodes work here 
c
      if (.not. restart) then

      if (numerical_single) then
c
c     Perform an energy+grad calculation at the central point and
c     save wavefunction information for rapid convergence
c
        if (.not. task_gradient(rtdb)) call errquit
     &      ('stpr_gen_hess_at: reference point energy+grad failed',7)
      else
c
c     Perform an energy calculation at the central point and
c     save wavefunction information for rapid convergence
c
        if (.not. task_energy(rtdb)) call errquit
     &      ('stpr_gen_hess_at: reference point energy failed',7,
     &       CALC_ERR)
      end if
        if (.not. task_save_state(rtdb,'hess')) call errquit
     &      ('stpr_gen_hess_at: failed to save central point',8,
     &       RTDB_ERR)
c
c     save reference geometry
c     
*
        if (o_stepper_debug) then
          call ga_sync()
          write(6,*)' going false 3',whoiam
          call util_flush(6)
        endif
*
        call ga_sync()
        ignore_status = rtdb_parallel(.false.)
        if (whoiam .eq. master) then
          if (.not.geom_create(geom,'geometry')) call errquit
     &        ('stpr_gen_hess_at:geom_create failed?',9, GEOM_ERR)
          if (o_debug_restart)
     &        write(6,*)' load 3 geometry ',
     &        '   node=',whoiam
          if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &        ('stpr_gen_hess_at:geom_rtdb_load failed?',10, RTDB_ERR)
          if (.not.geom_rtdb_store(rtdb,geom,'reference'))call errquit
     &        ('stpr_gen_hess_at: geom_rtdb_store failed',11, RTDB_ERR)
          if (.not.geom_destroy(geom)) call errquit
     &        ('stpr_gen_hess_at: geom_destroy failed?',12, GEOM_ERR)
        endif
        ignore_status = rtdb_parallel(.true.)
        call ga_sync()
*
        if (o_stepper_debug) then
          write(6,*)' gone true 3',whoiam
          call util_flush(6)
          call ga_sync()
        endif
*
      endif
c     
c     broadcast number of atoms
c     
      call ga_brdcst(Msg_gen_at_nat+MSGINT, nat,
     $     MA_sizeof(MT_INT,1,MT_BYTE), 0)
c     
*
      if (o_stepper_debug) then
        call ga_sync()
        write(6,*)' going false 4',whoiam
        call util_flush(6)
      endif
*
      call ga_sync()
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
        rank_hess = 3*nat
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'minus delta nuclear gradient vector',
     &      h_gradm,k_gradm)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',13, MA_ERR)
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'plus delta nuclear gradient vector',
     &      h_gradp,k_gradp)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',14, MA_ERR)
        if (.not. MA_Push_Get(MT_DBL,3*rank_hess,
     &      'derivative dipole moment',
     &      h_ddipole, k_ddipole)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed',15, MA_ERR)
        if (.not. MA_Push_Get(MT_DBL,(rank_hess*rank_hess),
     &      'nuclear hessian matrix',
     &      h_hess,k_hess)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',16, MA_ERR)
c
        if (restart) then 
          call stpr_get_genat_restart(rank_hess,
     &        dbl_mb(k_hess),
     &        dbl_mb(k_gradm),.false.,
     &        dipole_okay,dbl_mb(k_ddipole))
        else
          call dfill((rank_hess*rank_hess),0.0d00,
     &        dbl_mb(k_hess),1)
          call dfill((3*rank_hess),0.0d00,dbl_mb(k_ddipole),1)
        endif
          write(luout,*)' '
          write(luout,*)'initial hessian' 
          call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     &        rank_hess,rank_hess,1)
          write(luout,*)' '
 
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
*
      if (o_stepper_debug) then
        write(6,*)' gone true 4',whoiam
        call util_flush(6)
        call ga_sync()
      endif
*
c     
c     all nodes doing work here
c     
      if (numerical_single) then
c
c     Read gradient and dipole at the central point
c
      call ga_brdcst(Msg_gen_at_nat+MSGINT, rank_hess,
     $        MA_sizeof(MT_INT,1,MT_BYTE), 0)
c
c     me.ne.0 need gradient
c
      if(ga_nodeid().ne.0) then
         if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &        'minus delta nuclear gradient vector',
     &        h_gradm,k_gradm)) call errquit
     &        ('stpr_gen_hess_at: ma_push_get failed?',13, MA_ERR)
      endif
            if(.not. rtdb_get(rtdb,'task:gradient',MT_DBL,
     &          rank_hess,dbl_mb(k_gradm)))
     &          call errquit('stpr_gen_hess_at: rtdb_get failed',38,
     R     RTDB_ERR)
            if(.not. rtdb_get(rtdb,'task:dipole',mt_dbl,
     &          3,mdipole)) dipole_okay = .false.
      endif

      call ga_sync()
      if(.not.restart) then
        iatom_start = 1
        xyz_pass    = .true.
      else
        xyz_pass = .false.
*
        if (o_stepper_debug) then
          call ga_sync()
          write(6,*)' going false 5',whoiam
          call util_flush(6)
        endif
*
        call ga_sync()
        ignore_status = rtdb_parallel(.false.)
        if (whoiam.eq.master) then
          if (util_print('restart', print_low)) then
            write(luout,*)' **** gen_hess restart ****'
            write(luout,*)' iatom_start = ',iatom_start
            write(luout,*)' ixyz_start  = ',ixyz_start
            call util_flush(luout)
          endif
        endif
        ignore_status = rtdb_parallel(.true.)
        call ga_sync()
*
        if (o_stepper_debug) then
          write(6,*)' gone true 5',whoiam
          call util_flush(6)
          call ga_sync()
        endif
*
      endif
      call ga_sync()
      do 00100 iatom = iatom_start,nat
        if (.not.all_atom) then
          skip_atom = .true.
          do iat_test = 0,(nat_active-1)
            if (iatom.eq.int_mb(i_active+iat_test))
     &          skip_atom = .false.
          enddo
        else
          if (.not.geom_create(geom,'reference')) call errquit
     &        ('stpr_gen_hess_at:geom_create failed?',9, GEOM_ERR)
          if (o_debug_restart)
     &        write(6,*)' load 4 reference:',
     &        '  iatom=',iatom,
     &        '   node=',whoiam
          if (.not.geom_rtdb_load(rtdb,geom,'reference')) call errquit
     &        ('stpr_gen_hess_at:geom_rtdb_load failed?',10, RTDB_ERR)
          skip_atom = .not.(sym_atom(geom,iatom,q1))
          if (.not.geom_destroy(geom)) call errquit
     &        ('stpr_gen_hess_at: geom_destroy failed?',12, GEOM_ERR)
        endif
        if (skip_atom) goto 00201
        call stpr_ghat_set_active(rtdb,iatom)
        if (xyz_pass) then
          ixyz_start = 1
        else
          xyz_pass = .true.
        endif
        do 00200 ixyz = ixyz_start,3
          last_pass = (ixyz.eq.3).and.(iatom.eq.nat)
          if (ixyz.eq.1 .and. no_x) goto 00199
          if (ixyz.eq.2 .and. no_y) goto 00199
          if (ixyz.eq.3 .and. no_z) goto 00199
*
          if (o_stepper_debug) then
            call ga_sync()
            write(6,*)' going false 6',whoiam
            call util_flush(6)
          endif
*
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if(whoiam.eq.master) then
            if (o_directions) then
              call util_date(date_string)
              date_len = inp_strlen(date_string)
              write(luout,
     &            '(1x,a5,i4,1x,a4,i2,a3,1x,a10,f10.1,6x,a5,a27)')
     &            'atom:',iatom,
     &            'xyz:',ixyz,'(+)',
     &            'wall time:',util_wallsec(),
     &            'date:',date_string(1:date_len)
              call util_flush(luout)
            endif
c     
c     master node only
c     
            if (.not.geom_create(geom,'reference')) call errquit
     &          ('stpr_gen_hess_at:geom_create failed?',17, GEOM_ERR)
          if (o_stepper_debug)
     &          write(6,*)' load 5 reference'
            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_rtdb_load failed?',18, RTDB_ERR)
            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',19, GEOM_ERR)
            xyz(ixyz) = xyz(ixyz) + delta
            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',20, GEOM_ERR)
c     
c     copy modified geometry and store on rtdb
c     
            write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'plus'
            new_geom_name = 'reference '//tmpname(1:11)
            if (.not.geom_strip_sym(geom)) call errquit
     &          ('stpr_gen_hess_at:1: geom_strip_sym failed',911,
     &       GEOM_ERR)
            if (.not.
     &          geom_rtdb_store(rtdb,geom,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_store failed',21,
     &       RTDB_ERR)
            if (.not.geom_destroy(geom))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_destroy failed?',22, GEOM_ERR)
c     
c     set "geometry" to "reference atom/xyz"
c     
            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: rtdb_cput failed ?',23, RTDB_ERR)
          endif
          ignore_status = rtdb_parallel(.true.)
          call ga_sync()
*
          if (o_stepper_debug) then
            write(6,*)' gone true 6',whoiam
            call util_flush(6)
            call ga_sync()
          endif
*
c
c     restore starting information for central point
c     
          if (.not. task_restore_state(rtdb,'hess')) call errquit
     &        ('stpr_gen_hess_at: failed to restore central point',
     &        24, RTDB_ERR)
c     
          if (print_level.gt.print_low)
     &          call util_print_set_level(print_low) ! Minimize output
          time_taskg1 = (-1.0d00)*util_wallsec()
          if (.not. task_gradient(rtdb)) call errquit
     &        ('stpr_gen_hess_at: energy+grad failed',25, CALC_ERR)
          time_taskg1 = time_taskg1 + util_wallsec()
          if (print_level.gt.print_low)
     &          call util_print_set_level(print_level) ! Restore current print
c     
*
          if (o_stepper_debug) then
            call ga_sync()
            write(6,*)' going false 7',whoiam
            call util_flush(6)
          endif
*
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if (whoiam.eq.master) then
            if (.not. geom_rtdb_delete(rtdb,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_delete failed',26,
     &       RTDB_ERR)
            if(.not. rtdb_get(rtdb,'task:gradient',MT_DBL,
     &          rank_hess,dbl_mb(k_gradp)))
     &          call errquit('stpr_gen_hess_at: rtdb_get failed',27,
     &       RTDB_ERR)
c     zero non-active components
            call stpr_grad_act(rtdb,nat,dbl_mb(k_gradp))
            if(.not. rtdb_get(rtdb,'task:dipole',mt_dbl,
     &          3,pdipole)) dipole_okay = .false.
c
          endif

         call ga_sync()
         if (numerical_single) then
c     Skip calculations            
         else
          if (whoiam.eq.master) then
            if (o_directions) then
              call util_date(date_string)
              date_len = inp_strlen(date_string)
              write(luout,
     &            '(1x,a5,i4,1x,a4,i2,a3,1x,a10,f10.1,6x,a5,a27)')
     &            'atom:',iatom,
     &            'xyz:',ixyz,'(-)',
     &            'wall time:',util_wallsec(),
     &            'date:',date_string(1:date_len)
              call util_flush(luout)
            endif
c
c     master node only
c     
            if (.not.geom_create(geom,'reference')) call errquit
     &          ('stpr_gen_hess_at:geom_create failed?',28, GEOM_ERR)
          if (o_debug_restart)
     &          write(6,*)' load 6 reference',
     &          '   node=',whoiam
            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_rtdb_load failed?',29, GEOM_ERR)
            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',30, GEOM_ERR)
            xyz(ixyz) = xyz(ixyz) - delta
            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',31, GEOM_ERR)
c     
c     copy modified geometry and store on rtdb
c     
            write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'minus'
            new_geom_name = 'reference '//tmpname(1:11)
            if (.not.geom_strip_sym(geom)) call errquit
     &          ('stpr_gen_hess_at:2: geom_strip_sym failed',911,
     &       GEOM_ERR)
            if (.not.
     &          geom_rtdb_store(rtdb,geom,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_store failed',32,
     &       RTDB_ERR)
            if (.not.geom_destroy(geom))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_destroy failed?',33, GEOM_ERR)
c     
c     set "geometry" to "reference atom/xyz"
c     
            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: rtdb_cput failed ?',34, RTDB_ERR)
          endif
          ignore_status = rtdb_parallel(.true.)
          call ga_sync()
*
          if (o_stepper_debug) then
            write(6,*)' gone true 7',whoiam
            call util_flush(6)
            call ga_sync()
          endif
*
c
c     restore starting information from central point
c     
          if (.not. task_restore_state(rtdb,'hess')) call errquit
     &        ('stpr_gen_hess_at: failed to restore central point',
     &        35, RTDB_ERR)
c     
*          call util_print_push
          if (print_level.gt.print_low)
     &          call util_print_set_level(print_low) ! Minimize output
          time_taskg2 = (-1.0d00)*util_wallsec()

          if (.not. task_gradient(rtdb))
     &        call errquit('stpr_gen_hess_at: energy+grad failed',
     &        36, CALC_ERR)
          time_taskg2 = time_taskg2 + util_wallsec()
*          call util_print_pop
          if (print_level.gt.print_none)
     &          call util_print_set_level(print_level) ! Restore current print
c     
*
          if (o_stepper_debug) then
            call ga_sync()
            write(6,*)' going false 8',whoiam
            call util_flush(6)
          endif

c      endif corresponds to skipped calculations
       endif
*
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if (whoiam.eq.master) then
            if (.not. geom_rtdb_delete(rtdb,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_delete (2) failed',37,
     &       RTDB_ERR)

c      Read gradient at negative offset from disk
         if (numerical_single) then
c      Gradient is in k_gradm
         else
            if(.not. rtdb_get(rtdb,'task:gradient',MT_DBL,
     &          rank_hess,dbl_mb(k_gradm)))
     &          call errquit('stpr_gen_hess_at: rtdb_get failed',38,
     &       RTDB_ERR)
c
c     zero non-active components
            call stpr_grad_act(rtdb,nat,dbl_mb(k_gradm))
            if(.not. rtdb_get(rtdb,'task:dipole',mt_dbl,
     &          3,mdipole)) dipole_okay = .false.
         endif
c      Compute hessian
         if (numerical_single) then
            call stpr_fd_upd_hess(rtdb,
     &          dbl_mb(k_hess),dbl_mb(k_gradm),
     &          dbl_mb(k_gradp),1.0d00,delta,nat,iatom,ixyz)
            if (dipole_okay)
     &          call stpr_fd_upd_dipole(dbl_mb(k_ddipole),
     &          mdipole,pdipole,1.0d00,delta,nat,iatom,ixyz,q1)
         else
            call stpr_fd_upd_hess(rtdb,
     &          dbl_mb(k_hess),dbl_mb(k_gradm),
     &          dbl_mb(k_gradp),2.0d00,delta,nat,iatom,ixyz)
            if (dipole_okay)
     &          call stpr_fd_upd_dipole(dbl_mb(k_ddipole),
     &          mdipole,pdipole,2.0d00,delta,nat,iatom,ixyz,q1)
         endif
          endif
          ignore_status = rtdb_parallel(.true.)
          call ga_sync()
*
          if (o_stepper_debug) then
            write(6,*)' gone true 8',whoiam
            call util_flush(6)
            call ga_sync()
          endif
*
00199     continue
*
          if (o_stepper_debug) then
            call ga_sync()
            write(6,*)' going false 9',whoiam
            call util_flush(6)
          endif
*
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if (whoiam.eq.master) then
            call stpr_put_genat_restart(rank_hess,dbl_mb(k_hess),
     &          dbl_mb(k_gradm),iatom,ixyz,nat,.false.,
     &          dipole_okay,dbl_mb(k_ddipole))
          endif
          ignore_status = rtdb_parallel(.true.)
          call ga_sync()

*
          if (o_stepper_debug) then
            write(6,*)' gone true 9',whoiam
            call util_flush(6)
            call ga_sync()
            write(6,*)'@----------------------- iatom,ixyz,nodeid ',
     &          iatom,ixyz,whoiam
            call util_flush(6)
            call ga_sync()
          endif
*
          if (.not.last_pass) then
            time_taskg = int(1.2d00*(time_taskg1 + time_taskg2))
            if (util_test_time_remaining(rtdb,time_taskg)) then
              if (o_stepper_debug) write(luout,*)
     &            'remaining time enough for next step'
            else
              if (whoiam.eq.master) write(luout,*)
     &            'remaining time not enough for next step'
              stpr_gen_hess_at = .false.
              return
            endif
          endif
*
00200   continue
00201   continue
*
        if (o_stepper_debug) then
          call ga_sync()
          write(6,*)'@----------------------- iatom,nodeid ',
     &        iatom,whoiam
          call util_flush(6)
          call ga_sync()
        endif
*
00100 continue
c     
*
      if (o_stepper_debug) then
        call ga_sync()
        write(6,*)' going false 10',whoiam
        call util_flush(6)
      endif
*
      call ga_sync()
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
        if (.not.geom_create(geom,'reference')) call errquit
     &      ('stpr_gen_hess_at:geom_create failed?',1, GEOM_ERR)
          if (o_debug_restart)
     &      write(6,*)' load 7 reference',
     &      '   node=',whoiam
        if (.not.geom_rtdb_load(rtdb,geom,'reference')) call errquit
     &      ('stpr_gen_hess_at:geom_rtdb_load failed?',2, RTDB_ERR)
        call sym_hess_symmetrize(geom,dbl_mb(k_hess))
        if (dipole_okay) 
     &      call sym_der_dip_symmetrize(geom, dbl_mb(k_ddipole))
        if (.not.geom_destroy(geom))
     &      call errquit
     &      ('stpr_gen_hess_at: geom_destroy failed?',33, GEOM_ERR)
        if (all_atom) then
c     
c     average contributions from finite diff hess since they "should"
c     be the same
c     
          call stpr_gen_hess_foldave(dbl_mb(k_hess),rank_hess)
        else
          call stpr_gen_hess_foldave(dbl_mb(k_hess),rank_hess)
cedo          call stpr_gen_hess_fold(dbl_mb(k_hess),rank_hess)
c     
c     assumes zeroed hessian prior to fd computations
c     
          call stpr_gen_set_diag(dbl_mb(k_hess),rank_hess)  
        endif
        if (util_print('hessian', print_low)) then
          write(luout,*)' '
          write(luout,*)' finite difference hessian delta = ',delta
          call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     &        rank_hess,rank_hess,1)
          write(luout,*)' '
        endif
        if (dipole_okay.and.
     &      util_print('derivative dipole',print_low)) then
          write(tmpname,*)
     &        ' finite difference derivative dipole; delta = ',delta
          call stpr_print_ddipole(dbl_mb(k_ddipole),
     &        tmpname, nat, -1.0d00)
        endif
c
c     compress hessian if actlist is on
c
        if(.not. rtdb_get(rtdb,'hessian:compress',mt_log,
     &       1,hesscompress)) hesscompress = .false.
        if(hesscompress)then
           if (.not. MA_Push_Get(MT_DBL,9*nat_active*nat_active,
     &        'hesscompress',
     &        h_hessc,k_hessc)) call errquit
     &        ('stpr_gen_hess_at: ma_push_get failed?',13, MA_ERR)
          call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     &        rank_hess,rank_hess,1)
           call stpr_hesscompr(nat,nat_active,
     A          int_mb(i_active),dbl_mb(k_hess),dbl_mb(k_hessc))
           rank_hess=3*nat_active
          call output(dbl_mb(k_hessc),1,rank_hess,1,rank_hess,
     &        rank_hess,rank_hess,1)
           call dcopy(9*nat_active*nat_active,dbl_mb(k_hessc),1,
     C          dbl_mb(k_hess),1)
           if (.not.ma_pop_stack(h_hessc)) call errquit
     &      ('stpr_ghat_set_active: ma_pop k_hessc failed?',22, MA_ERR)
        endif
        call stpr_wrt_fd_from_sq(dbl_mb(k_hess),rank_hess,FILEHESS)
        if (dipole_okay) then
          call stpr_wrt_fd_dipole(dbl_mb(k_ddipole),nat,FILEDDIPOLE)
        else
          call util_file_unlink(FILEDDIPOLE)
        endif
        if (util_print('information', print_medium)) then
          write(luout,*)' triangle hessian written to ',
     &        FILEHESS(1:inp_strlen(FILEHESS))
          if (dipole_okay)
     &        write(luout,*)' derivative dipole written to ',
     &        FILEDDIPOLE(1:inp_strlen(FILEDDIPOLE))
        endif
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
*
      if (o_stepper_debug) then
        write(6,*)' gone true 10',whoiam
        call util_flush(6)
        call ga_sync()
      endif
*
c
c     restore and delete save state
c
      if (.not. task_restore_state(rtdb, 'hess')) call errquit
     &    ('stpr_gen_hess_at: failed to restore central point',40,
     &       RTDB_ERR)
      if (.not. task_delete_state(rtdb, 'hess')) call errquit
     &    ('stpr_gen_hess_at: failed to restore central point',41,
     &       RTDB_ERR)
*
      if (o_stepper_debug) then
        call ga_sync()
        write(6,*)' going false 11',whoiam
        call util_flush(6)
      endif
*
      call ga_sync()
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
c     
c     reset "geometry" to appropriate designation
c     
        ignore_status = rtdb_delete(rtdb,'geometry:actlist')
        ignore_status = rtdb_delete(rtdb,'geometry')  ! always try to delete geometry alias
        if (ogeom_aliased) then
          if(.not.rtdb_cget(rtdb,'stpr:gh:at:alias:name',1,tmpname))
     &        call errquit
     &        ('stpr_gh_at: get of <stpr:gh:at:alias:name> failed ?',
     &        39, RTDB_ERR)
          if(.not.rtdb_cput(rtdb,'geometry',1,tmpname))
     &        call errquit
     &        ('stpr_gen_hess_at: rtdb_cput failed ?',39, RTDB_ERR)
          ignore_status = rtdb_delete(rtdb,'stpr:gh:at:alias:name') ! always delete reference alias
        endif
c     
c     free memory
c     
        
        if(.not.ma_pop_stack(h_hess)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_hess) failed?',42,
     &       MA_ERR)
        if(.not.ma_pop_stack(h_ddipole)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_ddipole) failed?',43,
     &       MA_ERR)
        if(.not.ma_pop_stack(h_gradp)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_gradp) failed?',44,
     &       MA_ERR)
        if(.not.ma_pop_stack(h_gradm)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_gradm) failed?',45,
     &       MA_ERR)
      endif
      if (numerical_single.and.whoiam.ne.0) then
        if(.not.ma_pop_stack(h_gradm)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_gradm) failed?',145,
     &       MA_ERR)
      endif
      ignore_status = rtdb_parallel(.true.)
      call ga_sync()
*
      if (o_stepper_debug) then
        write(6,*)' gone true 11',whoiam
        call util_flush(6)
        call ga_sync()
      endif
*
      if(.not.all_atom) then
        if (.not.ma_free_heap(h_active)) call errquit
     &      ('stpr_gen_hess_at: ma_free_heap(h_active) failed?',46,
     &       MA_ERR)
      endif
*
      if (o_stepper_debug) then
        write(6,*)' out of stpr_gen_hess_at ',whoiam
      endif
*
      stpr_gen_hess_at = .true.
      call ga_sync()
      end
C>
C> @}
C> \ingroup stpr_priv
C> @{
      subroutine stpr_ghat_set_active(rtdb,iatom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "sym.fh"
#include "geom.fh"
#include "util.fh"
*::passed
      integer rtdb  ! [input] the rtdb handle
      integer iatom ! [input] the current loop atom
*::local
      logical current_rtdb_status, ignore_status
      double precision q2
      integer sga_act
      integer geom
      integer iiiat
      integer nat
      integer h_sga, k_sga
      logical odebug_print
*::code
      odebug_print = util_print('information',print_debug)
      odebug_print = odebug_print.or.
     &    util_print('debug_stepper',print_debug)
      odebug_print = odebug_print.or.
     &    util_print('debug_stepper_restart',print_debug)
      odebug_print = odebug_print .and. (ga_nodeid().eq.0)
*
      current_rtdb_status = rtdb_parallel(.false.)
      call ga_sync()
*
      if (ga_nodeid().eq.0) then
*
* generate active atom list based on symmetry for theory+gradients 
* this is orthoginal to user defined list for hessian generation.
*
        if (.not.geom_create(geom,'reference')) call errquit
     &      ('stpr_ghat_set_active:geom_create failed?',17, GEOM_ERR)
          if (odebug_print)
     &      write(6,*)' load 8 reference'
        if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &      call errquit
     &      ('stpr_ghat_set_active:geom_rtdb_load failed?',18, RTDB_ERR)
        if (.not. geom_ncent(geom,nat)) call errquit
     &      ('stpr_ghat_set_active:geom_ncent failed?',3, GEOM_ERR)
        if (.not.ma_push_get
     &      (mt_int,nat,'stpr:geom:actlist',h_sga,k_sga))
     &      call errquit
     &      ('stpr_ghat_set_active:ma for geom actlist failed',911,
     &       MA_ERR)
        sga_act = 0
        do iiiat=1,iatom
          if (sym_atom_pair(geom,iatom,iiiat,q2)) then
            int_mb(k_sga+sga_act) = iiiat
            sga_act = sga_act + 1
            if (odebug_print)
     &          write(6,*)'debug:T:',iatom,iiiat
          else
            if (odebug_print)
     &          write(6,*)'debug:F:',iatom,iiiat
          endif
        enddo
        if (sga_act.gt.0) then
          if (.not.rtdb_put(rtdb,'geometry:actlist',mt_int,
     &        sga_act, int_mb(k_sga))) call errquit
     &        ('stpr_ghat_set_active: rtdb put of actlist failed',911,
     &       RTDB_ERR)
        endif
        if (odebug_print)
     &      write(6,*)'debug: active atom list: ',sga_act,':',
     &      (int_mb(k_sga+iiiat),iiiat=0,(sga_act-1))
        if (.not.geom_destroy(geom))
     &      call errquit
     &      ('stpr_ghat_set_active: geom_destroy failed?',22, GEOM_ERR)
        if (.not.ma_pop_stack(h_sga)) call errquit
     &      ('stpr_ghat_set_active: ma_pop k_sga failed?',22, MA_ERR)
      endif
      call ga_sync()
      ignore_status = rtdb_parallel(current_rtdb_status)
      end
c     zero non-active components of grad
      subroutine stpr_grad_act(rtdb,nat,grad)
      implicit none
      integer rtdb,nat
      double precision grad(3,*)
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer nactive   ! number of active atoms
      integer h_actlist ! ma handle for actlist
      integer k_actlist ! ma index for actlist
      integer ma_type
      integer ii,i0,iact
      logical dozero
      integer ga_nodeid
      external ga_nodeid
c
       if (.not.rtdb_ma_get(rtdb,'gen_hess:actlist', ma_type, 
     $     nactive,h_actlist)) return
C       if(.true.)return
       if (nactive.le.0 .or. nactive.gt.nat) call errquit
     &    ('stpr_grr: invalid number of active atoms',nactive,
     &       INPUT_ERR)
      if (.not.MA_get_index(h_actlist,k_actlist)) call errquit
     &    ('stpr_grr: ma_get_index failed for actlist',911, MA_ERR)
      i0=1
      do ii=1,nat
         dozero=.true.
         do iact=i0,nactive
            if(ii.eq.int_mb(k_actlist+iact-1)) then
               dozero=.false.
               i0=iact+1
            endif
         enddo
         if(dozero) then
            grad(1,ii)=0d0
            grad(2,ii)=0d0
            grad(3,ii)=0d0
         endif
      enddo
      if (.not.MA_free_heap(h_actlist)) call errquit
     C     ('stpr_grr: mafreeheap failed',0,MA_ERR)
      return
      end
      subroutine stpr_hesscompr(nat,nactive,
     A          i_active,hessin,hessout)
      implicit none
      integer nat,nactive
      integer i_active(*)
      double precision hessin(3,nat,3,nat)
      double precision hessout(3,nactive,3,nactive)
c
      integer iatom,jatom
      integer icount,jcount,ix,jx
c
      do iatom=1,nactive
         icount=i_active(iatom)
         do jatom=1,nactive
            jcount=i_active(jatom)
            do ix=1,3
               do jx=1,3
                  hessout(jx,jatom,ix,iatom)=hessin(jx,jcount,ix,icount)
               enddo
            enddo
         enddo
      enddo
      return
      end
C>
C> @}
