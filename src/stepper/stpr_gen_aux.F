*
* $Id$
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
C> \ingroup stpr_priv
C> @{
C>
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
C>
C> \brief Update the Hessian from new gradient data
C>
C> Apply the finite difference formula to the gradient data provided
C> to calculate a contribution to and update the Hessian.
C>
C> Depending on the values of the parameters `gradm`, `gradp` and
C> `s_delta` the subroutine will evaluate either the forward difference
C> equation
C> \f{eqnarray*}{
C>    f'(x) &=& \frac{f(x+h)-f(x)}{h}
C> \f}
C> or the central difference formula
C> \f{eqnarray*}{
C>    f'(x) &=& \frac{f(x+h)-f(x-h)}{2h}
C> \f}
C>
      subroutine stpr_fd_upd_hess(rtdb,
     &    hess,gradm,gradp,s_delta,delta,nat,
     &    iatom_t,ixyz_t)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
#include "sym.fh"
c::passed
      integer rtdb                       !< [Input] The RTDB handle
      integer nat                        !< [Input] The number of atoms
      integer iatom_t                    !< [Input] The displaced atom
      integer ixyz_t                     !< [Input] The displaced coordinate
      double precision hess(3,nat,3,nat) !< [In/Output] The Hessian
      double precision gradm(3,nat)      !< [Input] The gradient at 
                                         !< negative (or no) displacement
      double precision gradp(3,nat)      !< [Input] The gradient at 
                                         !< positive displacement
      double precision delta             !< [Input] The step size
      double precision s_delta           !< [Input] The scale factor
                                         !< for total displacement
                                         !< * 1 - for forward difference
                                         !< * 2 - for central difference
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
C>
C> \brief Check and retrieve restart data
C>
C> Check the runtime database for the state of the calculation
C> and report back whether this calculation needs to be restarted
C> or not. The function also provides the rank of the initial
C> atom and atomic coordinate to start from.
C>
C> The initial atom and atomic coordinate are read from the
C> restart file. For this purpose the function may open the restart
C> file if available. But the restart file must be closed before
C> the function returns.
C>
C> \return Returns .TRUE. if the calculation needs to be restarted,
C> and .FALSE. otherwise.
C>
      logical function stpr_check_genat_restart(
     &    rtdb, iatom_start,ixyz_start)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "cstprfiles.fh"
      integer rtdb        !< [Input] The RTDB handle
      integer iatom_start !< [Output] The atom to start from
      integer ixyz_start  !< [Output] The atomic coordinate to start
                          !< from
c
      integer int_restart
      integer rank, ijunk1, ijunk2
      integer mitob1
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
            write(luout,*)'*** Warning continue called for but no  ***'
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
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(Msg_gen_at_iatom  +MSGINT,iatom_start,mitob1,0)
      call ga_brdcst(Msg_gen_at_ixyz   +MSGINT,ixyz_start, mitob1,0)
      call ga_brdcst(Msg_gen_at_restart+MSGINT,int_restart,mitob1,0)
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
C>
C> \brief Load the restart information
C>
C> Load the data needed for restarting a Hessian calculation. 
C>
C> Obviously the restart data is read from the restart file. Hence,
C> this subroutine must open the restart file, but the restart file
C> must be closed before the routine returns.
C>
      subroutine stpr_get_genat_restart(rank_in,hess,grad0,get_grad0,
     &    dipole_okay,ddipole)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
      integer rank_in    !< [Input] The dimension of the Hessian
      double precision hess(rank_in,rank_in) !< [Output] The Hessian
      double precision grad0(rank_in)        !< [Output] The gradient
      double precision ddipole(3*rank_in)    !< [Output] The dipole
                                             !< derivative
      logical get_grad0    !< [Input] Should the gradient be loaded
                           !< * `.TRUE.` - load the gradient
                           !< * `.FALSE.` - do not load the gradient
      logical dipole_okay  !< [Output] Is the dipole derivative returned
                           !< * `.TRUE.` - The dipole derivative is 
                           !< returned
                           !< * `.FALSE.` - The dipole derivative is 
                           !< not returned
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
        close(unit=69,status='keep')
      else
        write(6,*)' no finite difference hessian restart ',
     &      'information read '
      endif
      end
C>
C> \brief Store the restart data
C>
C> Write the restart data to a file. The restart file is opened and
C> overwritten with new restart data. The restart file must be closed
C> before the subroutine returns.
C>
      subroutine stpr_put_genat_restart(rank,hess,grad0,
     &    iatom_in,ixyz_in,nat,put_grad0,
     &    dipole_okay,ddipole)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "cstprfiles.fh"
      integer rank  !< [Input] The dimension of the Hessian
      integer iatom_in !< [Input] The current atom
      integer ixyz_in  !< [Input] The current atomic coordinate
      integer nat      !< [Input] The number of atoms
      double precision hess(rank,rank) !< [Input] The Hessian
      double precision grad0(rank)     !< [Input] The gradient
      double precision ddipole(3*rank) !< [Input] The dipole derivative
      logical put_grad0    !< [Input] Should the gradient be stored?
      logical dipole_okay  !< [Input] Should the dipole derivative be
                           !< stored?
      logical lopen        !< Is the file open?
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
        call util_file_unlink(FILEATR)
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
c     If unit 69 is not closed then the OPEN statement will be ignored
c     in which case the restart data is appended to the existing file
c     rather than overwriting it.
c     Therefore we force the file to be closed. If it is not then there
c     is a bug somewhere else.
      inquire(unit=69,opened=lopen)
      if (lopen) call errquit("stpr_put_genat_restart: unit 69 should "
     &                      //"be closed at this point",0,UERR)
      open(unit=69,file=FILEATR,
     &    form='unformatted',
     &    access='sequential',
     &    status='unknown')
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
C>
C> @}
