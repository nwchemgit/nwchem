*
* $Id: stpr_gen_aux.F,v 1.6 1997-04-22 20:01:45 d3e129 Exp $
* 
* this file contains auxilary routines to stpr_gen functions
*
* current routines: 
* stpr_fd_upd_hess             ! computes finite difference either central or forward
* stpr_wrt_fd_from_sq          ! writes hessian to file
* stpr_check_genat_restart     ! check for restart "is info available?"
* stpr_get_genat_restart       ! get restart info
* stpr_put_genat_restart       ! put restart info out to restart file
* stpr_gen_hess_foldave        ! averages off diaginal contributions
* stpr_gen_hess_fold           ! sums off diaginal contributions (partial computations)
* stpr_gen_set_diag            ! sets inactive atom diagonal contribs to large value
*
c
      subroutine stpr_fd_upd_hess(hess,gradm,gradp,s_delta,delta,nat,
     &    iatom_t,ixyz_t)
      implicit none
#include "stdio.fh"
c::passed
      integer nat
      integer iatom_t, ixyz_t
      double precision hess(3,nat,3,nat)
      double precision gradm(3,nat)
      double precision gradp(3,nat)
      double precision delta, s_delta
c::local
      integer iatom,ixyz
      double precision rdelta, value
c
c 
c finite difference  [g(x+delta) - g(x-delta)]/(s_delta*delta) (s_delta = 2.0)
c central difference [g(x+delta) - g(x)]/(s_delta*delta) (s_delta = 1.0)
c 
c
      rdelta = 1.0d00/(s_delta*delta)
      do 00100 iatom = 1,nat
        do 00200 ixyz = 1,3
          value = rdelta*(gradp(ixyz,iatom)-gradm(ixyz,iatom))
          hess(ixyz_t,iatom_t,ixyz,iatom) = value
          hess(ixyz,iatom,ixyz_t,iatom_t) = value
00200   continue
00100 continue
c
      end
      subroutine stpr_wrt_fd_from_sq(hess,rank_hess,filename)
      implicit none
c
#include "stdio.fh"
#include "inp.fh"
c
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
      character*(*) filename
c
      logical does_it_exist
c
      integer i, j, lu
c
      lu = 66
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if (does_it_exist)
     &    write(luout,*)
     &    ' stpr_wrt_fd_from_sq: overwrite of existing file:',
     &    filename(1:inp_strlen(filename))
      open(unit=lu,file=filename,
     &    form='formatted',
     &    access='sequential',
     &    status='unknown')
c
      do 00100 i = 1,rank_hess
        do 00200 j = 1,i
          write(lu,10000)hess(i,j)
00200   continue
00100 continue
c
10000 format(1x,1pd20.10)
c
      close(unit=lu,status='keep')
c
      end
      logical function stpr_check_genat_restart(
     &    rtdb, iatom_start,ixyz_start)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "cstprfiles.fh"
      integer rtdb
      integer iatom_start
      integer ixyz_start
c
      integer int_restart
      integer rank, ijunk
      logical does_it_exist
c
      logical ostart, orestart, ocontinue
c
c assume not a restart
      iatom_start = 1
      ixyz_start  = 1
      int_restart = 0
      call util_get_rtdb_state(rtdb,ostart,ocontinue,orestart)
      if (ostart) then
        call util_file_unlink(FILEATR)
      else if (ocontinue.or.orestart) then
        if (ga_nodeid().eq.0) then
          does_it_exist = .false.
          inquire(file=FILEATR,exist=does_it_exist)
          if (does_it_exist) then
            open(unit=69,file=FILEATR,
     &          form='unformatted',
     &          access='sequential',
     &          status='old')
            read(69)iatom_start,ixyz_start,rank, ijunk
            close(unit=69,status='keep')
            int_restart = 1
          else
            write(luout,*)'*** Warning restart called for but no   ***'
            write(luout,*)'*** fd restart file for nuclear hessian ***'
            write(luout,*)'*** starting from scratch so to speak   ***'
            iatom_start = 1
            ixyz_start  = 1
            int_restart = 0
          endif
        endif
      else
        call errquit
     &      ('stpr_check_genat_restart: error with rtdb state',911)
      endif
c
      call ga_brdcst(Msg_gen_at_iatom  +MSGINT,iatom_start,mitob(1),0)
      call ga_brdcst(Msg_gen_at_ixyz   +MSGINT,ixyz_start, mitob(1),0)
      call ga_brdcst(Msg_gen_at_restart+MSGINT,int_restart,mitob(1),0)
c
      if (int_restart.eq.1) then
        stpr_check_genat_restart = .true.
      else if (int_restart.eq.0) then
        stpr_check_genat_restart = .false.
      else
        write(luout,*)' invalid int_restart value ', ga_nodeid()
        call errquit(' stpr_check_genat_restart: error ',int_restart)
      endif
      end
      subroutine stpr_get_genat_restart(rank_in,hess,grad0,get_grad0)
      implicit none
#include "global.fh"
#include "stdio.fh"
#include "cstprfiles.fh"
      integer rank_in
      double precision hess(rank_in,rank_in)
      double precision grad0(rank_in)
      logical get_grad0
c
      logical does_it_exist
      integer ijunk1, ijunk2, rank, iflag_grad0
c
      if (ga_nodeid().ne.0) then
        write(luout,*)' non-master node called me ',ga_nodeid()
        call errquit('stpr_get_genat_restart: error ',911)
      endif
c
      inquire(file='gen_at_restart',exist=does_it_exist)
      if (does_it_exist) then
        open(unit=69,file=FILEATR,
     &      form='unformatted',
     &      access='sequential',
     &      status='old')
        read(69)ijunk1,ijunk2,rank,iflag_grad0
        if (rank.ne.rank_in) then
          write(luout,*)'rank not the same as rank_in '
          write(luout,*)' rank    :',rank
          write(luout,*)' rank_in :',rank_in
          close(unit=69,status='keep')
          call errquit('stpr_get_genat_restart: error ',911)
        endif
        if (get_grad0.and.iflag_grad0.ne.1) then
          write(luout,*)' grad 0 not written but requested '
          call errquit(' stpr_get_genat_restart: error',911)
        endif
        if ((.not.get_grad0).and.iflag_grad0.eq.1) then
          write(luout,*)' grad 0 written but not requested '
          call errquit(' stpr_get_genat_restart: error',911)
        endif
        if (get_grad0) read(69) grad0
        read(69) hess
      endif
      end
      subroutine stpr_put_genat_restart(rank,hess,grad0,
     &    iatom_in,ixyz_in,nat,put_grad0)
      implicit none
#include "global.fh"
#include "stdio.fh"
#include "cstprfiles.fh"
      integer rank
      integer iatom_in, ixyz_in, nat
      double precision hess(rank,rank)
      double precision grad0(rank)
      logical put_grad0
c
      integer iatom_start, ixyz_start
      integer iflag_grad0
c
      if (ga_nodeid().ne.0) then
        write(luout,*)' non-master node called me ',ga_nodeid()
        call errquit('stpr_put_genat_restart: error ',911)
      endif
c
      call util_file_unlink(fileatr)
      if(iatom_in.eq.nat.and.ixyz_in.eq.3) return
      iatom_start = iatom_in
      ixyz_start  = ixyz_in + 1
      if(ixyz_in.eq.3) then
        iatom_start = iatom_start + 1
        ixyz_start  = 1
      endif
c
      if (put_grad0) then
        iflag_grad0 = 1
      else
        iflag_grad0 = 0
      endif
      open(unit=69,file=FILEATR,
     &    form='unformatted',
     &    access='sequential',
     &    status='new')
      write(69)iatom_start,ixyz_start,rank,iflag_grad0
      if (put_grad0) write(69)grad0
      write(69)hess
      close(unit=69,status='keep')
      end
      subroutine stpr_gen_hess_foldave(hess,rank_hess)
*! averages off diaginal contributions
      implicit none
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
*
      integer i,j
      double precision dbl_tmp 
*
      do i = 1,rank_hess
        do j = 1,(i-1)
          dbl_tmp   = hess(i,j) + hess(j,i)
          dbl_tmp   = dbl_tmp/2.0d00
          hess(i,j) = dbl_tmp
          hess(j,i) = dbl_tmp
        enddo
      enddo
      end
      subroutine stpr_gen_hess_fold(hess,rank_hess)
*! sums off diaginal contributions assuming a partial computation
      implicit none
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
*
      integer i,j
      double precision dbl_tmp 
#define DEBUG
#if defined(DEBUG)      
      integer icount
      double precision dbl_diff, max_dbl_diff
c
      write(6,*)' hessian before fold operaton'
      call output(hess,1,rank_hess,1,rank_hess,
     &    rank_hess,rank_hess,1)
      icount = 0
      max_dbl_diff = -1.0d00
#endif
*
      do i = 1,rank_hess
        do j = 1,(i-1)
#if defined(DEBUG)      
          dbl_diff = abs(hess(i,j)) - abs(hess(j,i))
          max_dbl_diff = max(max_dbl_diff,dbl_diff)
          icount = icount + 1
          write(6,12345)icount,dbl_diff,max_dbl_diff
12345     format('<',i2,'> <diff=',f14.8,'> <diff_max=',f14.8)
#endif
          dbl_tmp   = hess(i,j) + hess(j,i)
          hess(i,j) = dbl_tmp
          hess(j,i) = dbl_tmp
        enddo
      enddo
#if defined(DEBUG)      
      write(6,*)' hessian after fold operaton'
      call output(hess,1,rank_hess,1,rank_hess,
     &    rank_hess,rank_hess,1)
#endif
      end
      subroutine stpr_gen_set_diag(hess,rank_hess)
*! sets diag to default value for stiff frequency analysis
*! e.g., active atom computation
      implicit none
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
*
      integer i
      double precision dbl_tmp 
#define DEBUG
#if defined(DEBUG)      
      write(6,*)' hessian before diag set operaton'
      call output(hess,1,rank_hess,1,rank_hess,
     &    rank_hess,rank_hess,1)
#endif
*
      dbl_tmp = 1.0d00
      do i = 1,rank_hess
        if (hess(i,i).eq.0.0d00) hess(i,i) = dbl_tmp
      enddo
#if defined(DEBUG)      
      write(6,*)' hessian after diag set operaton'
      call output(hess,1,rank_hess,1,rank_hess,
     &    rank_hess,rank_hess,1)
#endif
      end
