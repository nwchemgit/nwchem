*
* $Id: stpr_gen_aux.F,v 1.14 2006-02-08 16:30:51 windus Exp $
* 
* this file contains auxilary routines to stpr_gen functions
*
* current routines: 
* stpr_fd_upd_dipole           ! computes finite difference dipole moment
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
      subroutine stpr_fd_upd_dipole(ddipole,mdipole,pdipole,
     &    s_delta,delta,nat,iatom,ixyz,q1)
      implicit none
c
c::passed
      integer nat
      integer iatom
      integer ixyz
      double precision s_delta
      double precision delta
      double precision mdipole(3)
      double precision pdipole(3)
      double precision ddipole(3,3,nat)
      double precision q1
c:local
      integer moment
      double precision rdelta
      double precision value
c
      rdelta = 1.0d00/(s_delta*delta)
      do moment = 1,3
        value = rdelta*(pdipole(moment)-mdipole(moment)) * q1
        ddipole(moment,ixyz,iatom) = value
      enddo
c
      end
      subroutine stpr_fd_upd_hess(rtdb,
     &    hess,gradm,gradp,s_delta,delta,nat,
     &    iatom_t,ixyz_t)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "sym.fh"
c::passed
      integer rtdb
      integer nat
      integer iatom_t, ixyz_t
      double precision hess(3,nat,3,nat)
      double precision gradm(3,nat)
      double precision gradp(3,nat)
      double precision delta, s_delta
c::local
      integer geom
      integer iatom,ixyz
      double precision rdelta, value
      double precision q2
c
      if (.not.geom_create(geom,'reference')) call errquit
     $    ('stpr_fd_upd_hess:geom_create failed?',1, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'reference')) call errquit
     $        ('stpr_fd_upd_hess:geom_rtdb_load failed?',2, RTDB_ERR)
c 
c finite difference  [g(x+delta) - g(x-delta)]/(s_delta*delta) (s_delta = 2.0)
c central difference [g(x+delta) - g(x)]/(s_delta*delta) (s_delta = 1.0)
c 
c
      rdelta = 1.0d00/(s_delta*delta)
      do 00100 iatom = 1,iatom_t
        if (sym_atom_pair(geom,iatom_t,iatom,q2)) then
**           write(6,*) ' iatom_t iatom q2 ', iatom_t, iatom, q2
           if (iatom.ne.iatom_t) q2 = q2 + q2
          do 00200 ixyz = 1,3
            value = rdelta*(gradp(ixyz,iatom)-gradm(ixyz,iatom))
            value = q2*value
            hess(ixyz_t,iatom_t,ixyz,iatom) = value
**            hess(ixyz,iatom,ixyz_t,iatom_t) = value
00200     continue
        endif
00100 continue
c
      if (.not.geom_destroy(geom))
     $    call errquit
     $    ('stpr_fd_upd_hess: geom_destroy failed?',33, GEOM_ERR)
      end
      subroutine stpr_wrt_fd_dipole(ddipole,nat,filename)
      implicit none
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
      integer nat
      double precision ddipole(3,3,nat)
      character*(*) filename
c
      integer lu
      integer print_level
      integer atom, xyz, moment
      logical does_it_exist
c
      call util_print_get_level(print_level)
      lu = 67 
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if ((does_it_exist).and.(print_level.gt.print_none))
     &    write(luout,*)
     &    'stpr_wrt_fd_dipole: overwrite of existing file',
     &    filename(1:inp_strlen(filename))
      open(unit=lu,file=filename,
     &    form='formatted',
     &    access='sequential',
     &    status='unknown')
c
      do atom = 1,nat
        do  xyz = 1,3
          do moment = 1,3
            write(lu,10000)ddipole(moment,xyz,atom)
          enddo
        enddo
      enddo
c
10000 format(1x,1pd20.10)
c
      close(unit=lu,status='keep')
c
      end
      subroutine stpr_wrt_fd_from_sq(hess,rank_hess,filename)
      implicit none
c
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
c
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
      character*(*) filename
c
      logical does_it_exist
c
      integer i, j, lu
      integer print_level
c
      call util_print_get_level(print_level)
      lu = 66
      does_it_exist = .false.
      inquire(file=filename,exist=does_it_exist)
      if ((does_it_exist).and.(print_level.gt.print_none))
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
#include "errquit.fh"
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
      integer rank, ijunk1, ijunk2
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
        if (ga_nodeid().eq.0)
     &    call util_file_unlink(FILEATR)
      else if (ocontinue.or.orestart) then
        if (ga_nodeid().eq.0) then
          does_it_exist = .false.
          inquire(file=FILEATR,exist=does_it_exist)
          if (does_it_exist) then
            open(unit=69,file=FILEATR,
     &          form='unformatted',
     &          access='sequential',
     &          status='old')
            read(69)iatom_start,ixyz_start,rank, ijunk1, ijunk2
            close(unit=69,status='keep')
            int_restart = 1
          else
            write(luout,*)'*** Warning continue called for but no   ***'
            write(luout,*)'*** fd restart file for nuclear hessian ***'
            write(luout,*)'*** starting from scratch so to speak   ***'
          endif
        endif
      else
        call errquit
     &      ('stpr_check_genat_restart: error with rtdb state',911,
     &       RTDB_ERR)
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
        call errquit(' stpr_check_genat_restart: fatal error ',
     &      int_restart, INPUT_ERR)
      endif
      end
      subroutine stpr_get_genat_restart(rank_in,hess,grad0,get_grad0,
     &    dipole_okay,ddipole)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
      integer rank_in
      double precision hess(rank_in,rank_in)
      double precision grad0(rank_in)
      double precision ddipole(3*rank_in)
      logical get_grad0
      logical dipole_okay
c
      logical does_it_exist
      integer ijunk1, ijunk2, rank, iflag_grad0
      integer dipole_there
c
      if (ga_nodeid().ne.0) then
        write(luout,*)' non-master node called me ',ga_nodeid()
        call errquit('stpr_get_genat_restart: error ',911, INPUT_ERR)
      endif
c
      inquire(file=FILEATR,exist=does_it_exist)
      if (does_it_exist) then
        open(unit=69,file=FILEATR,
     &      form='unformatted',
     &      access='sequential',
     &      status='old')
        read(69)ijunk1,ijunk2,rank,iflag_grad0,dipole_there
        if (dipole_there.eq.1) then
          dipole_okay = .true.
        else
          dipole_okay = .false.
        endif
        if (rank.ne.rank_in) then
          write(luout,*)'rank not the same as rank_in '
          write(luout,*)' rank    :',rank
          write(luout,*)' rank_in :',rank_in
          close(unit=69,status='keep')
          call errquit('stpr_get_genat_restart: error ',911, INPUT_ERR)
        endif
        if (get_grad0.and.iflag_grad0.ne.1) then
          write(luout,*)' grad 0 not written but requested '
          call errquit(' stpr_get_genat_restart: error',911, INPUT_ERR)
        endif
        if ((.not.get_grad0).and.iflag_grad0.eq.1) then
          write(luout,*)' grad 0 written but not requested '
          call errquit(' stpr_get_genat_restart: error',911, INPUT_ERR)
        endif
        if (get_grad0) read(69) grad0
        read(69) hess
        if (dipole_okay) read(69) ddipole
        if (util_print('debug_stepper_restart',print_debug)
     &      .or.
     &      util_print('debug_stepper',print_debug)) then
          write(6,*)'hessian read from restart file '
          call output(hess,1,rank,1,rank,rank,rank,1)
          call stpr_print_ddipole(ddipole,
     &        'dipole derivative read from restart file',
     &        (rank/3),
     &        1.0d-07)
        endif
      else
        write(6,*)' no finite difference hessian restart ',
     &      'information read '
      endif
      end
      subroutine stpr_put_genat_restart(rank,hess,grad0,
     &    iatom_in,ixyz_in,nat,put_grad0,
     &    dipole_okay,ddipole)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
      integer rank
      integer iatom_in, ixyz_in, nat
      double precision hess(rank,rank)
      double precision grad0(rank)
      double precision ddipole(3*rank)
      logical put_grad0
      logical dipole_okay
c
      integer iatom_start, ixyz_start
      integer iflag_grad0
      integer dipole_there
c
      if (ga_nodeid().ne.0) then
        write(luout,*)' non-master node called me ',ga_nodeid()
        call errquit('stpr_put_genat_restart: error ',911, INPUT_ERR)
      endif
c
      if(iatom_in.eq.nat.and.ixyz_in.eq.3) then
        call util_file_unlink(fileatr)
        return
      endif
      iatom_start = iatom_in
      ixyz_start  = ixyz_in + 1
      if(ixyz_in.eq.3) then
        iatom_start = iatom_start + 1
        ixyz_start  = 1
      endif
c
      dipole_there = 0
      if (dipole_okay) dipole_there = 1
c
      if (put_grad0) then
        iflag_grad0 = 1
      else
        iflag_grad0 = 0
      endif
      call util_file_unlink(fileatr)
      open(unit=69,file=FILEATR,
     &    form='unformatted',
     &    access='sequential',
     &    status='new')
      write(69)iatom_start,ixyz_start,rank,iflag_grad0,dipole_there
      if (put_grad0) write(69)grad0
      write(69)hess
      if (dipole_okay)write(69) ddipole
      close(unit=69,status='keep')
      if (util_print('debug_stepper_restart',print_debug)
     &    .or.
     &    util_print('debug_stepper',print_debug)) then
        write(6,*)'hessian put to restart file '
        call output(hess,1,rank,1,rank,rank,rank,1)
        call stpr_print_ddipole(ddipole,
     &      'dipole derivative put to restart file',
     &      (rank/3),
     &      1.0d-07)
      endif
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
#include "util.fh"
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
*
      integer i,j
      double precision dbl_tmp 
      integer icount
      double precision dbl_diff, max_dbl_diff
      logical o_debug
c
      o_debug = util_print('debug_stepper_restart',print_debug)
      o_debug = o_debug .or.
     &    util_print('debug_stepper',print_debug)
      if (o_debug) then
        write(6,*)' hessian before fold operaton'
        call output(hess,1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
        icount = 0
        max_dbl_diff = -1.0d00
      endif
*
      do i = 1,rank_hess
        do j = 1,(i-1)
          if (o_debug) then
            dbl_diff = abs(hess(i,j)) - abs(hess(j,i))
            max_dbl_diff = max(max_dbl_diff,dbl_diff)
            icount = icount + 1
            write(6,12345)icount,dbl_diff,max_dbl_diff
          endif
          dbl_tmp   = hess(i,j) + hess(j,i)
          hess(i,j) = dbl_tmp
          hess(j,i) = dbl_tmp
        enddo
      enddo
      if (o_debug) then
        write(6,*)' hessian after fold operaton'
        call output(hess,1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
      endif
12345 format('<',i2,'> <diff=',f14.8,'> <diff_max=',f14.8)
      end
      subroutine stpr_gen_set_diag(hess,rank_hess)
*! sets diag to default value for stiff frequency analysis
*! e.g., active atom computation
      implicit none
#include "util.fh"
      integer rank_hess
      double precision hess(rank_hess,rank_hess)
*
      integer i
      double precision dbl_tmp 
      logical o_debug
      o_debug = util_print('debug_stepper_restart',print_debug)
      o_debug = o_debug .or.
     &    util_print('debug_stepper',print_debug)
      if (o_debug) then
        write(6,*)' hessian before diag set operaton'
        call output(hess,1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
      endif
*
      dbl_tmp = 1.0d00
      do i = 1,rank_hess
        if (hess(i,i).eq.0.0d00) hess(i,i) = dbl_tmp
      enddo
      if (o_debug) then
        write(6,*)' hessian after diag set operaton'
        call output(hess,1,rank_hess,1,rank_hess,
     &      rank_hess,rank_hess,1)
      endif
      end
      subroutine stpr_print_ddipole(ddipole,msg,nat,thresh)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
      integer nat
      double precision ddipole(3,3,nat)
      double precision thresh
      character*(*) msg
      character*1 cname(3)
      character*1 kapname
c
      double precision val
      integer moment
      integer xyz, atom
c
      cname(1) = 'x'
      cname(2) = 'y'
      cname(3) = 'z'
c
      write(luout,'(/a/)') msg(1:inp_strlen(msg))
      write(luout,*)' '
      do moment = 1,3 
        kapname = cname(moment)
        call inp_ucase(kapname)
        write(luout,*)' '
        write(luout,10001) kapname
        do atom = 1,nat
          do xyz = 1,3
            if (moment.ge.1.and.moment.le.3) then
              val = ddipole(moment,xyz,atom)
*
* from A Physicists Desk Reference, The Second Edition of Physics Vade Mecum
*      Herbert L. Anderson, Editor in Chief
*      Copyright (C) 1989 American Institute of Physics
*      335 East 45th Street, New York, NY 10017
*
*1 debye = 10**(-18) esu cm * [1 e/4.8032068 x 10**(-10) esu]*[1 m /100cm]*[a0/5.29177249 x 10**(-11) m]
*1 debye = (1.0/4.8032068/5.29177249) * 10**(-18 + 10 - 2 + 11) e a0
*1 debye = (1.0/4.8032068/5.29177249) * 10**(1) e a0
*1 e a0  = (4.8032068*5.29177249) * 10**(-1) debye
*1 e a0  = 25.417477608 * 10**(-1) debye
*1 e a0  = 2.5417477608 debye
*
*use 1 e a0 = 2.541 7478 debye
*
              val = val*2.5417478d00
              val = val/0.529177249d00  ! bohr->angstrom matches current geom data
              if (abs(val).gt.thresh) then
                write(luout,10002)
     &              cname(moment),atom,cname(xyz),
     &              ddipole(moment,xyz,atom),val
              endif
            else
              write(luout,10000)moment
              call errquit('stpr_print_ddipole: fatal error',911,
     &       UNKNOWN_ERR)
            endif
          enddo
        enddo
      enddo
      write(luout,*)' '
      write(luout,*)' '
10000 format('invalid moment value',i10)
10001 format(1x,a1,1x,
     &    'vector of derivative dipole (au) [debye/angstrom]')
10002 format(1x,'d_dipole_',a1,'/<atom=',i4,',',a1,'> = ',
     &    f10.4,5x,'[',f10.4,']')
      end
