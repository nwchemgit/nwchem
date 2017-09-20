      subroutine ccsd_lambda(
     1           d_f1,d_v2,d_d1,d_q1,d_o1,
     2           d_t1,d_t2,d_lambda1,d_lambda2,
     4           k_f1_offset,k_v2_offset,
     5           k_d1_offset,k_q1_offset,k_o1_offset,
     6           k_t1_offset,k_t2_offset,
     7           k_l1_offset,k_l2_offset,
     8           size_l1,size_l2,
     9           dipole_scf,    dipole_cor,
     1           quadrupole_scf,quadrupole_cor,
     2           octupole_scf,  octupole_cor)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_diis.fh"
#include "tce_restart.fh"
c
      integer iter
      integer d_f1,d_v2,d_t1,d_t2
      integer d_lambda1,d_lambda2
      integer d_lr1,d_lr2
      integer k_f1_offset,k_v2_offset
      integer k_t1_offset,k_t2_offset
      integer k_l1_offset,k_l2_offset
      integer size_l1,size_l2
c DM
      integer l_hh_offset,k_hh_offset,size_hh
      integer l_hp_offset,k_hp_offset,size_hp
      integer l_ph_offset,k_ph_offset,size_ph
      integer l_pp_offset,k_pp_offset,size_pp
      integer d_hh,d_hp,d_ph,d_pp
c
      double precision ref,corr
      double precision r1,r2
      double precision residual
      double precision ddotfile
      integer d_d1(3)           ! MO dipole integral file
      integer d_q1(6)           ! MO quadrupole integral file
      integer d_o1(10)          ! MO octupole integral file
      integer k_d1_offset(3)    ! Offset for d1 integral file
      integer k_q1_offset(6)    ! Offset for q1 integral file
      integer k_o1_offset(10)   ! Offset for o1 integral file
      double precision dipole_scf(3)     ! Dipole moment - SCF XYZ
      double precision dipole_cor(3)     ! Dipole moment - CC XYZ
      double precision quadrupole_scf(6) ! Quadrupole moment
      double precision quadrupole_cor(6) ! Quadrupole moment
      double precision octupole_scf(10)  ! Octupole moment
      double precision octupole_cor(10)  ! Octupole moment
      double precision cpu, wall
      external ddotfile
      logical nodezero
      integer dummy
      character*255 filename
c
      nodezero = (ga_nodeid().eq.0)
      if (skip_l) then
        if (nodezero) write(LuOut,9700)
        return
      endif
c
      nodezero=(ga_nodeid().eq.0)
      irrep_y=0
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1)) write(LuOut,9400) "CCSD Lambda"
        call tce_filename('lr1',filename)
        call createfile(filename,d_lr1,size_l1)
        call ccsd_lambda1(d_f1,d_lr1,d_t1,d_t2,d_v2,d_lambda1,
     1     d_lambda2,k_f1_offset,k_l1_offset,k_t1_offset,
     2     k_t2_offset,k_v2_offset,k_l1_offset,k_l2_offset)
        call reconcilefile(d_lr1,size_l1)
        call tce_filename('lr2',filename)
        call createfile(filename,d_lr2,size_l2)
        call ccsd_lambda2(d_f1,d_lr2,d_t1,d_t2,d_v2,d_lambda1,
     1     d_lambda2,k_f1_offset,k_l2_offset,k_t1_offset,
     2     k_t2_offset,k_v2_offset,k_l1_offset,k_l2_offset)
        call reconcilefile(d_lr2,size_l2)
        r1 = ddotfile(d_lr1,d_lr1,size_l1)
        r2 = ddotfile(d_lr2,d_lr2,size_l2)
        residual = max(r1,r2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) then
            write(LuOut,9410)
          endif
c
c  density matrix calculations: 1-body density matrix
c 
c
c     ====================
c     1PDM Hole-Hole Block
c     ====================
      if(nodezero) then 
       write(LuOut,9531)
       call util_flush(LuOut)
      endif
      call tce_filename('hh',filename)
      call tce_dens_hh_offset(l_hh_offset,k_hh_offset,size_hh)
      call createfile(filename,d_hh,size_hh)
      call ccsd_1pdm_hh_mo(d_hh,d_t1,d_t2,
     1     d_lambda1,d_lambda2,
     2     k_hh_offset,k_t1_offset,k_t2_offset,
     1     k_l1_offset,k_l2_offset)
      call reconcilefile(d_hh,size_hh)
      call tce_dens_hh_print(d_hh,k_hh_offset,printtol,irrep_y)
      call deletefile(d_hh)
      if (.not.ma_pop_stack(l_hh_offset))
     1    call errquit("ccsd_lambda: MApopstack failed",29,MA_ERR)
c     ============================
c     1PDM Particle-Particle Block
c     ============================
      if(nodezero) then
       write(LuOut,9532)
       call util_flush(LuOut)
      endif
      call tce_filename('pp',filename)
      call tce_dens_pp_offset(l_pp_offset,k_pp_offset,size_pp)
      call createfile(filename,d_pp,size_pp)
      call ccsd_1pdm_pp_mo(d_pp,d_t1,d_t2,
     1     d_lambda1,d_lambda2,
     2     k_pp_offset,k_t1_offset,k_t2_offset,
     3     k_l1_offset,k_l2_offset)
      call reconcilefile(d_pp,size_pp)
      call tce_dens_pp_print(d_pp,k_pp_offset,printtol,irrep_y)
      call deletefile(d_pp)
      if (.not.ma_pop_stack(l_pp_offset))
     1    call errquit("ccsd_lambda: MApopstack failed",30,MA_ERR)
c     ========================
c     1PDM Hole-Particle Block
c     ========================
      if(nodezero) then
       write(LuOut,9533)
       call util_flush(LuOut)
      endif
      call tce_filename('hp',filename)
      call tce_dens_hp_offset(l_hp_offset,k_hp_offset,size_hp)
      call createfile(filename,d_hp,size_hp)
      call ccsd_1pdm_hp_mo(d_hp,d_lambda1,k_hp_offset,k_l1_offset)
      call reconcilefile(d_hp,size_hp)
      call tce_dens_hp_print(d_hp,k_hp_offset,printtol,irrep_y)
      call deletefile(d_hp)
      if (.not.ma_pop_stack(l_hp_offset))
     1    call errquit("ccsd_lambda: MApopstack failed",32,MA_ERR)
c     ========================
c     1PDM Particle-Hole Block
c     ========================
      if(nodezero) then
       write(LuOut,9534)
       call util_flush(LuOut)
      endif
      call tce_filename('ph',filename)
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,size_ph)
      call createfile(filename,d_ph,size_ph)
      call ccsd_1pdm_ph_mo(d_ph,d_t1,d_t2,d_lambda1,d_lambda2,
     1     k_ph_offset,k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset)
      call reconcilefile(d_ph,size_ph)
      call tce_dens_ph_print(d_ph,k_ph_offset,printtol,irrep_y)
      call deletefile(d_ph)
      if (.not.ma_pop_stack(l_ph_offset))
     1    call errquit("ccsd_lambda: MApopstack failed",33,MA_ERR)
c
      if(nodezero) then
       write(LuOut,9535)
       call util_flush(6)
      endif
c
c
          if (ampnorms) then
            r1 = ddotfile(d_lambda1,d_lambda1,size_l1)
            r2 = ddotfile(d_lambda2,d_lambda2,size_l2)
            if (nodezero) then
              write(LuOut,9082) "Lambda singles",r1
              write(LuOut,9082) "Lambda doubles",r2
            endif
          endif
          call deletefile(d_lr2)
          call deletefile(d_lr1)
          call tce_diis_tidy()
c        endif
          if(nodezero) then 
           write(LuOut,9002)
           call util_flush(6)
          end if
          call tce_print_y1(d_lambda1,k_l1_offset,printtol,irrep_y)
          call tce_print_y2(d_lambda2,k_l2_offset,printtol,irrep_y)
          if (save_l(1)) then
            if(nodezero) then
            write(LuOut,*) 'Saving L1 now...'
            endif
            call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                           handle_l1,irrep_y)
          endif
          if (save_l(2)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L2 now...'
            endif
            call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                           handle_l2,irrep_y)
          endif
          if((perturbative.eq.'lambda(t)').or.
     1       (perturbative.eq.'lambda[t]').or.
     2       (perturbative.eq.'2_t').or.
     3       (perturbative.eq.'2_tq')) go to 4511
          call ccsd_dipole(
     1           d_f1,d_v2,d_d1,d_q1,d_o1,
     2           d_t1,d_t2,d_lambda1,d_lambda2,
     4           k_f1_offset,k_v2_offset,
     5           k_d1_offset,k_q1_offset,k_o1_offset,
     6           k_t1_offset,k_t2_offset,
     7           k_l1_offset,k_l2_offset,
     8           size_l1,size_l2,
     9           dipole_scf,    dipole_cor,
     1           quadrupole_scf,quadrupole_cor,
     2           octupole_scf,  octupole_cor)
 4511     continue
          return
        endif
      if (save_l(1).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving L1 now...'
          endif
          call y1_restart_save(d_lambda1,k_l1_offset,size_l1,0,
     1                         handle_l1,irrep_y)
        endif
        if (save_l(2).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving L2 now...'
          endif
          call y2_restart_save(d_lambda2,k_l2_offset,size_l2,0,
     1                         handle_l2,irrep_y)
        endif
        call tce_diis(.true.,iter,.true.,.true.,.false.,.false.,
     1                d_lr1,d_lambda1,k_l1_offset,size_l1,
     2                d_lr2,d_lambda2,k_l2_offset,size_l2,
     3                dummy,dummy,dummy,dummy,
     4                dummy,dummy,dummy,dummy)
        call deletefile(d_lr2)
        call deletefile(d_lr1)
        if (nodezero) call util_flush(LuOut)
      enddo
      call errquit('ccsd_lambda: maxiter exceeded',iter,CALC_ERR)
      return
 9002 format('Lambda-CCSD leading excitations:')
 9400 format(/,1x,A,' iterations',/,
     1  1x,'---------------------------------------------',/
     2  1x,'Iter          Residuum            Cpu    Wall',/
     3  1x,'---------------------------------------------')
 9410 format(
     1  1x,'---------------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f25.13,2f8.1)
 9070 format(1x,A,' correlation energy / hartree = ',f25.15)
 9080 format(1x,A,' total energy / hartree       = ',f25.15)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9700 format(1x,'Skipping CCSD-Lambda altogether')
c
 9531 format('1-PDM CCSD : h-h block')
 9532 format('1-PDM CCSD : p-p block')
 9533 format('1-PDM CCSD : h-p block')
 9534 format('1-PDM CCSD : p-h block')
 9535 format('1-PDM CCSD : all done',//)
c
      end

      subroutine ccsd_dipole(
     1           d_f1,d_v2,d_d1,d_q1,d_o1,
     2           d_t1,d_t2,d_lambda1,d_lambda2,
     4           k_f1_offset,k_v2_offset,
     5           k_d1_offset,k_q1_offset,k_o1_offset,
     6           k_t1_offset,k_t2_offset,
     7           k_l1_offset,k_l2_offset,
     8           size_l1,size_l2,
     9           dipole_scf,    dipole_cor,
     1           quadrupole_scf,quadrupole_cor,
     2           octupole_scf,  octupole_cor)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_diis.fh"
c
      integer axis
      integer d_f1,d_v2,d_t1,d_t2
      integer d_lambda1,d_lambda2
      integer d_lr1,d_lr2
      integer k_f1_offset,k_v2_offset
      integer k_t1_offset,k_t2_offset
      integer k_l1_offset,k_l2_offset
      integer size_l1,size_l2
      integer k_s0_offset,l_s0_offset,size_s0 ! any scalar
      integer d_d0                       ! Dipole moment
      integer d_q0                       ! Quadrupole moment
      integer d_o0                       ! Octupole moment
      integer d_d1(3)                    ! MO dipole integral file
      integer d_q1(6)                    ! MO quadrupole integral file
      integer d_o1(10)                   ! MO octupole integral file
      integer k_d1_offset(3)             ! Offset for d1 integral file
      integer k_q1_offset(6)             ! Offset for q1 integral file
      integer k_o1_offset(10)            ! Offset for o1 integral file
      double precision dipole_cor_t      ! Dipole moment - CC total
      double precision dipole_exc_t      ! Dipole moment - EOM-CC total
      double precision dipole_scf(3)     ! Dipole moment - SCF XYZ
      double precision dipole_cor(3)     ! Dipole moment - CC XYZ
      double precision dipole_exc(3)     ! Dipole moment - EOM-CC
      double precision cdipole           ! Dipole moment correlation
      double precision quadrupole_scf(6) ! Quadrupole moment
      double precision quadrupole_cor(6) ! Quadrupole moment
      double precision quadrupole_exc(6) ! Quadrupole moment
      double precision cquadrupole       ! Quadrupole moment correlation
      double precision octupole_scf(10)  ! Octupole moment
      double precision octupole_cor(10)  ! Octupole moment
      double precision octupole_exc(10)  ! Octupole moment
      double precision coctupole         ! Octupole moment correlation
      double precision debye,buck,octo
      parameter (octo=0.711688d0)
      parameter (buck=1.344911d0)
      parameter (debye=2.54176568d0)
      double precision cpu, wall
      logical nodezero
      integer dummy
      integer sym_abelian_axis
      external sym_abelian_axis
      character*255 filename
      call tce_e_offset(l_s0_offset,k_s0_offset,size_s0)
      nodezero=(ga_nodeid().eq.0)
      irrep_y=0
      call tce_diis_init()
      do axis=1,3
        irrep_d=sym_abelian_axis(geom,axis)
        call tce_filename('d0',filename)
        call createfile(filename,d_d0,1) ! size_d0 = 1
        call ccsd_density1(d_d1(axis),d_d0,d_t1,d_t2,d_lambda1,
     1    d_lambda2,k_d1_offset(axis),k_s0_offset,k_t1_offset,
     2    k_t2_offset,k_l1_offset,k_l2_offset)
        call reconcilefile(d_d0,1)
        call get_block(d_d0,cdipole,1,0)
        dipole_cor(axis)=dipole_scf(axis)+cdipole
        if (nodezero.and.util_print('multipole',print_debug)) then
            write(LuOut,1001) axis,dipole_scf(axis),
     1                             dipole_cor(axis)
        endif
        call deletefile(d_d0)
      enddo
      dipole_cor_t = dsqrt( dipole_cor(1)*dipole_cor(1) +
     1                      dipole_cor(2)*dipole_cor(2) +
     2                      dipole_cor(3)*dipole_cor(3) )
      if (nodezero) write(LuOut,9430) "CCSD",
     1              dipole_cor(1),dipole_cor(1)*debye,
     2              dipole_cor(2),dipole_cor(2)*debye,
     3              dipole_cor(3),dipole_cor(3)*debye,
     4              dipole_cor_t ,dipole_cor_t*debye
      if (multipole.gt.1) then
        do axis=1,6
          irrep_d=sym_abelian_axis(geom,axis+3)
          call tce_filename('q0',filename)
          call createfile(filename,d_q0,1)
          call ccsd_density1(d_q1(axis),d_q0,d_t1,d_t2,d_lambda1,
     1      d_lambda2,k_q1_offset(axis),k_s0_offset,k_t1_offset,
     2      k_t2_offset,k_l1_offset,k_l2_offset)
          call reconcilefile(d_q0,1)
          call get_block(d_q0,cquadrupole,1,0)
          quadrupole_cor(axis)=quadrupole_scf(axis)-cquadrupole
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,1002) axis,quadrupole_scf(axis),
     1                               quadrupole_cor(axis)
          endif
          call deletefile(d_q0)
        enddo
        if (nodezero) write(LuOut,9432) 'CCSD','second',
     1                quadrupole_cor(1),quadrupole_cor(1)*buck,
     2                quadrupole_cor(2),quadrupole_cor(2)*buck,
     3                quadrupole_cor(3),quadrupole_cor(3)*buck,
     4                quadrupole_cor(4),quadrupole_cor(4)*buck,
     5                quadrupole_cor(5),quadrupole_cor(5)*buck,
     6                quadrupole_cor(6),quadrupole_cor(6)*buck
      endif
      if (multipole.gt.2) then
        do axis=1,10
          irrep_d=sym_abelian_axis(geom,axis+9)
          call tce_filename('o0',filename)
          call createfile(filename,d_o0,1)
          call ccsd_density1(d_o1(axis),d_o0,d_t1,d_t2,d_lambda1,
     1      d_lambda2,k_o1_offset(axis),k_s0_offset,k_t1_offset,
     2      k_t2_offset,k_l1_offset,k_l2_offset)
          call reconcilefile(d_o0,1)
          call get_block(d_o0,coctupole,1,0)
          octupole_cor(axis)=octupole_scf(axis)+coctupole
          if (nodezero.and.util_print('multipole',print_debug)) then
              write(LuOut,1003) axis,octupole_scf(axis),
     1                               octupole_cor(axis)
          endif
          call deletefile(d_o0)
        enddo
        if (nodezero) write(LuOut,9433) 'CCSD','third',
     1                octupole_cor(1), octupole_cor(1)*octo,
     2                octupole_cor(2), octupole_cor(2)*octo,
     3                octupole_cor(3), octupole_cor(3)*octo,
     4                octupole_cor(4), octupole_cor(4)*octo,
     5                octupole_cor(5), octupole_cor(5)*octo,
     2                octupole_cor(6), octupole_cor(6)*octo,
     3                octupole_cor(7), octupole_cor(7)*octo,
     4                octupole_cor(8), octupole_cor(8)*octo,
     5                octupole_cor(9), octupole_cor(9)*octo,
     6                octupole_cor(10),octupole_cor(10)*octo
      endif
      if (.not.ma_pop_stack(l_s0_offset))
     1     call errquit("ccsd_lambda: MA problem",307,MA_ERR)
      return
 1001 format(1x,'dipole:   axis=',i1,' SCF=',f15.7,' CCSD=',f15.7)
 1002 format(1x,'quadpole: axis=',i1,' SCF=',f15.7,' CCSD=',f15.7)
 1003 format(1x,'octpole:  axis=',i2,' SCF=',f15.7,' CCSD=',f15.7)
 9430 format(/,1x,A,' dipole moments / hartree & Debye',/,
     1  1x,'------------------------------------',/
     2  1x,'  X   ',2f15.7,/
     3  1x,'  Y   ',2f15.7,/
     4  1x,'  Z   ',2f15.7,/
     5  1x,'Total ',2f15.7,/
     6  1x,'------------------------------------')
 9432 format(/,1x,A,' ',A,' moments / hartree & Buckingham',/,
     1  1x,'--------------------------------',/
     2  1x,'XX ',2f15.7,/
     3  1x,'XY ',2f15.7,/
     4  1x,'XZ ',2f15.7,/
     5  1x,'YY ',2f15.7,/
     6  1x,'YZ ',2f15.7,/
     7  1x,'ZZ ',2f15.7,/
     8  1x,'--------------------------------')
 9433 format(/,1x,A,' ',A,' moments / hartree & Octos',/,
     1  1x,'--------------------------------',/
     2  1x,'XXX ',2f15.7,/
     3  1x,'XXY ',2f15.7,/
     4  1x,'XXZ ',2f15.7,/
     5  1x,'XYY ',2f15.7,/
     6  1x,'XYZ ',2f15.7,/
     7  1x,'XZZ ',2f15.7,/
     8  1x,'YYY ',2f15.7,/
     9  1x,'YYZ ',2f15.7,/
     1  1x,'YZZ ',2f15.7,/
     2  1x,'ZZZ ',2f15.7,/
     3  1x,'--------------------------------')
 9531 format('1-body CC density matrix: h-h block: correlated part',/)
 9532 format('1-body CC density matrix: p-p block: correlated part',/)
 9533 format('1-body CC density matrix: h-p block: correlated part',/)
 9534 format('1-body CC density matrix: p-h block: correlated part',/)
 9535 format('1-body CC density matrix: all done',/)
      end
c $Id$
