      logical function stpr_gen_hess_at(rtdb,delta,flag)
      implicit none
* $Id: stpr_gh_at.F,v 1.2 1996-02-16 00:47:35 d3e129 Exp $
* routine to generate hessian at current geometry using default basis set
* ra kendall 7/95 PNL
*
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "context.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "inp.fh"
#include "stdio.fh"
#include "stpr_optimize.fh"
c::functions
      logical stpr_check_genat_restart
      logical stpr_genat_mode_context_push
      logical stpr_genat_mode_context_pop
      external stpr_check_genat_restart
      external stpr_genat_mode_context_push
      external stpr_genat_mode_context_pop
c::passed
      integer rtdb
      double precision delta
      integer flag
c::local
      integer whoiam, master
      integer geom
      integer nat, rank_hess
      integer h_gradm, k_gradm ! handle/index for central gradient
      integer h_gradp, k_gradp ! handle/index for delta gradient
      integer h_hess, k_hess   ! handle/index for hessian matrix
      integer iatom, ixyz
      integer iatom_start, ixyz_start
      integer h_active, i_active, iat_test, nat_active
      logical all_atom, skip_atom
      logical ignore_status
      logical restart
      logical xyz_pass
      integer vecin_len
      double precision xyz(3), chg
      character*255 tmpname
      character*255 vecin, vecout, vecin_file, vecout_file
      character*16 tag_name
      character*40 new_geom_name
      logical dft_mode, scf_mode
c
      stpr_gen_hess_at = .false.
      dft_mode = flag.eq.stpr_dft_flag
      scf_mode = flag.eq.stpr_scf_flag
      whoiam = ga_nodeid()
      master = 0
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
c
c.. create/load reference geometry
        if (.not.geom_create(geom,'geometry')) call errquit
     &      ('stpr_gen_hess_at:geom_create failed?',1)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &      ('stpr_gen_hess_at:geom_rtdb_load failed?',2)
c..  get the number of atoms
        if (.not. geom_ncent(geom,nat)) call errquit
     &      ('stpr_gen_hess_at:geom_ncent failed?',3)
c
c.. copy reference geometry and store on rtdb
        if (.not.
     &      geom_rtdb_store(rtdb,geom,'reference'))
     &      call errquit
     &      ('stpr_gen_hess_at: geom_rtdb_store failed',5)
        if (.not.geom_destroy(geom))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_destroy failed?',20)
c.. set vector scope for computation (assumes central scf point has been computed)
        if (scf_mode) then
          vecin  = 'scf:input vectors'
          vecout = 'scf:output vectors'
          if (.not.rtdb_cget(rtdb,vecout,1,vecin_file)) call errquit
     &        ('stpr_gen_hess_at: output vectors not set',911)
          vecin_len = inp_strlen(vecin_file)
          vecout_file = vecin_file(1:vecin_len)//'.step'
        else if (dft_mode) then
          vecin  = 'dft:input vectors'
          vecout = 'dft:output vectors'
* vecout=t -> vecin_file is okay
* vecout=f -> vecin=t -> vecin_file okay
*          -> vecin=f -> vecin_file movecs        
          if (rtdb_cget(rtdb,vecout,1,vecin_file)) then
            continue
          else
            if (.not.rtdb_cget(rtdb,vecin,1,vecin_file))
     &        vecin_file = 'movecs'
          endif
          vecin_len = inp_strlen(vecin_file)
          vecout_file = vecin_file(1:vecin_len)//'.step'
        else
            call errquit
     &          ('stpr_gen_hess_at: not scf or dft mode ',flag)
        endif
      endif
      ignore_status = rtdb_parallel(.true.)
      if (.not.rtdb_ma_get(rtdb,'gen_hess:active atoms',MT_INT,
     &       nat_active,h_active)) then
        all_atom = .true.
      else
        all_atom = .false.
        if (.not.MA_get_index(h_active,i_active)) call errquit
     &         ('stpr_gen_hess_at: ma_get_index error',911)
      endif
c all nodes work here 
      restart = .false.
      restart = stpr_check_genat_restart(iatom_start,ixyz_start)
c
c..      broadcast number of atoms
c
      call ga_brdcst(Msg_gen_at_nat+MSGINT, nat, mitob(1), 0)
c
      call ga_sync()
c
      if (whoiam.eq.master) then
        rank_hess = 3*nat
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'cent nuclear gradient vector',
     &      h_gradm,k_gradm)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',10)
        if (.not. MA_Push_Get(MT_DBL,rank_hess,
     &      'delta nuclear gradient vector',
     &      h_gradp,k_gradp)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',11)
        if (.not. MA_Push_Get(MT_DBL,(rank_hess*rank_hess),
     &      'nuclear hessian matrix',
     &      h_hess,k_hess)) call errquit
     &      ('stpr_gen_hess_at: ma_push_get failed?',12)
c
        if (restart) then 
          call stpr_get_genat_restart(rank_hess,
     &        dbl_mb(k_hess),
     &        dbl_mb(k_gradm),.false.)
        else
          call dfill((rank_hess*rank_hess),0.0d00,
     &        dbl_mb(k_hess),1)
c
        endif
c
      endif
c
c... all nodes doing work here
      call ga_sync()
      if(.not.restart) then
        iatom_start = 1
        xyz_pass    = .true.
      else
        xyz_pass = .false.
        if (whoiam.eq.master) then
          write(luout,*)' **** gen_hess restart ****'
          write(luout,*)' iatom_start = ',iatom_start
          write(luout,*)' ixyz_start  = ',ixyz_start
          call util_flush(luout)
        endif
      endif
      do 00100 iatom = iatom_start,nat
        if (.not.all_atom) then
          skip_atom = .true.
          do iat_test = 0,(nat_active-1)
            if (iatom.eq.int_mb(i_active+iat_test)) skip_atom = .false.
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
          ignore_status = rtdb_parallel(.false.)
          if(whoiam.eq.master) then
            write(luout,*)' iatom = ',iatom, ' ixyz = ',ixyz,' plus '
c... master node only
            if (.not.geom_create(geom,'reference')) call errquit
     &          ('stpr_gen_hess_at:geom_create failed?',14)
            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_rtdb_load failed?',15)
            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',16)
            xyz(ixyz) = xyz(ixyz) + delta
            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',17)
c.. copy modified geometry and store on rtdb
            write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'plus '
            new_geom_name = 'reference '//tmpname(1:11)
            if (.not.
     &          geom_rtdb_store(rtdb,geom,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_store failed',19)
            if (.not.geom_destroy(geom))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_destroy failed?',20)
c set "geometry" to "reference atom/xyz"
            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: rtdb_cput failed ?',911)
c set input and output vectors
            if(.not.rtdb_cput(rtdb,vecin,1,vecin_file)) call errquit
     &          ('stpr_gen_hess_at: failed to set input vectors',911)
            if(.not.rtdb_cput(rtdb,vecout,1,vecout_file)) call errquit
     &          ('stpr_gen_hess_at: failed to set output vectors',911)
          endif
          call ga_sync()
          ignore_status = rtdb_parallel(.true.)
          if (scf_mode) then
            call scf_gradient(rtdb)
          else if (dft_mode) then
            call dft_grad_top(rtdb)
          else
            call errquit
     &          ('stpr_gen_hess_at: not scf or dft mode ',flag)
          endif
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if (whoiam.eq.master) then
*... set context
            if (.not.stpr_genat_mode_context_push(flag))
     &          call errquit
     &          ('stpr_gen_hess_at: mode_context failed',911)
            if(.not.context_prefix('gradients',tmpname))
     &          call errquit
     &          ('stpr_gen_hess_at:context_prefix failed ',911)
            if(.not. rtdb_get(rtdb,tmpname,MT_DBL,
     &          rank_hess,dbl_mb(k_gradp)))
     &          call errquit('stpr_gen_hess_at: rtdb_get failed',21)
c... remove Prefix context
            if (.not.stpr_genat_mode_context_pop(flag))
     &          call errquit
     &          ('stpr_gen_hess_at: context_pop failed',911)
c...
            write(luout,*)' iatom = ',iatom, ' ixyz = ',ixyz,' minus'
c... master node only
            if (.not.geom_create(geom,'reference')) call errquit
     &          ('stpr_gen_hess_at:geom_create failed?',14)
            if (.not.geom_rtdb_load(rtdb,geom,'reference'))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_rtdb_load failed?',15)
            if (.not.geom_cent_get(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',16)
            xyz(ixyz) = xyz(ixyz) - delta
            if (.not.geom_cent_set(geom,iatom,tag_name,xyz,chg))
     &          call errquit
     &          ('stpr_gen_hess_at:geom_cent_get failed?',17)
c.. copy modified geometry and store on rtdb
            write(tmpname,'(i3,1x,i1,1x,a5)')iatom,ixyz,'minus'
            new_geom_name = 'reference '//tmpname(1:11)
            if (.not.
     &          geom_rtdb_store(rtdb,geom,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_rtdb_store failed',19)
            if (.not.geom_destroy(geom))
     &          call errquit
     &          ('stpr_gen_hess_at: geom_destroy failed?',20)
c set "geometry" to "reference atom/xyz"
            if(.not.rtdb_cput(rtdb,'geometry',1,new_geom_name))
     &          call errquit
     &          ('stpr_gen_hess_at: rtdb_cput failed ?',911)
c set input and output vectors
            if(.not.rtdb_cput(rtdb,vecin,1,vecin_file)) call errquit
     &          ('stpr_gen_hess_at: failed to set input vectors',911)
            if(.not.rtdb_cput(rtdb,vecout,1,vecout_file)) call errquit
     &          ('stpr_gen_hess_at: failed to set output vectors',911)
          endif
          call ga_sync()
          ignore_status = rtdb_parallel(.true.)
          if (scf_mode) then
            call scf_gradient(rtdb)
          else if (dft_mode) then
            call dft_grad_top(rtdb)
          else
            call errquit
     &          ('stpr_gen_hess_at: not scf or dft mode ',flag)
          endif
          call ga_sync()
          ignore_status = rtdb_parallel(.false.)
          if (whoiam.eq.master) then
*... set context
            if (.not.stpr_genat_mode_context_push(flag))
     &          call errquit
     &          ('stpr_gen_hess_at: mode_context failed',911)
            if(.not.context_prefix('gradients',tmpname))
     &          call errquit
     &          ('stpr_gen_hess_at:context_prefix failed ',911)
            if(.not. rtdb_get(rtdb,tmpname,MT_DBL,
     &          rank_hess,dbl_mb(k_gradm)))
     &          call errquit('stpr_gen_hess_at: rtdb_get failed',21)
c... remove Prefix context
            if (.not.stpr_genat_mode_context_pop(flag))
     &          call errquit
     &          ('stpr_gen_hess_at: context_pop failed',911)
c...
            call stpr_fd_upd_hess(dbl_mb(k_hess),dbl_mb(k_gradm),
     &          dbl_mb(k_gradp),2.0d00,delta,nat,iatom,ixyz)
            call stpr_put_genat_restart(rank_hess,dbl_mb(k_hess),
     &          dbl_mb(k_gradm),iatom,ixyz,nat,.false.)
          endif
          ignore_status = rtdb_parallel(.true.)
00200   continue
00201   continue
00100 continue
c
      ignore_status = rtdb_parallel(.false.)
      if (whoiam.eq.master) then
        if (all_atom) then
* average contributions from finite diff hess since they "should" be the same
          call stpr_gen_hess_foldave(dbl_mb(k_hess),rank_hess)
        else
          call stpr_gen_hess_fold(dbl_mb(k_hess),rank_hess)
*.............................  assumes zeroed hessian prior to fd computations
          call stpr_gen_set_diag(dbl_mb(k_hess),rank_hess)  
        endif
        write(luout,*)' finite difference hessian delta = ',delta
        call output(dbl_mb(k_hess),1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
        call stpr_wrt_fd_from_sq(dbl_mb(k_hess),rank_hess,'hess_ascii')
        write(luout,*)' triangle hessian written to hess_ascii'
        
c set "geometry" to "geometry"
        if(.not.rtdb_cput(rtdb,'geometry',1,'geometry'))
     &      call errquit
     &      ('stpr_gen_hess_at: rtdb_cput failed ?',911)
c
c... free memory
        if(.not.ma_pop_stack(h_hess)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_hess) failed?',27)
        if(.not.ma_pop_stack(h_gradp)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_gradp) failed?',27)
        if(.not.ma_pop_stack(h_gradm)) call errquit
     &      ('stpr_gen_hess_at: ma_pop_stack(h_gradm) failed?',27)
      endif
      if(.not.all_atom) then
        if (.not.ma_free_heap(h_active)) call errquit
     &         ('stpr_gen_hess_at: ma_free_heap(h_active) failed?',27)
      endif
      call ga_sync()
      ignore_status = rtdb_parallel(.true.)
      stpr_gen_hess_at = .true.
      end
