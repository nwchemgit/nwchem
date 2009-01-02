      subroutine ccsd_lambda(d_d0,d_f1,d_v2,d_d1,d_t1,d_t2,
     1           d_lambda1,d_lambda2,
     2           k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,
     4           k_l1_offset,k_l2_offset,
     5           size_l1,size_l2,dipole_scf,dipole_cor)
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
      integer d_d0,d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_lambda1,d_lambda2
      integer d_lr1,d_lr2
      integer k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_l1_offset,k_l2_offset
      integer size_l1,size_l2
      double precision ref,corr
      double precision cpu, wall
      double precision r1,r2
      double precision residual
      double precision ddotfile
      double precision dipole_scf(3),dipole_cor(3)
      external ddotfile
      logical nodezero
      integer dummy
      character*255 filename
c
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
          call ccsd_dipole(d_d0,d_f1,d_v2,d_d1,d_t1,d_t2,
     1         d_lambda1,d_lambda2,
     2         k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3         k_t1_offset,k_t2_offset,
     4         k_l1_offset,k_l2_offset,
     5         size_l1,size_l2,dipole_scf,dipole_cor)
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
      end

      subroutine ccsd_dipole(d_d0,d_f1,d_v2,d_d1,d_t1,d_t2,
     1           d_lambda1,d_lambda2,
     2           k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,
     4           k_l1_offset,k_l2_offset,
     5           size_l1,size_l2,dipole_scf,dipole_cor)
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
      integer d_d0,d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_lambda1,d_lambda2
      integer d_lr1,d_lr2
      integer k_d0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_l1_offset,k_l2_offset
      integer size_l1,size_l2
      double precision cpu, wall
      double precision dipole_scf(3),dipole_cor(3)
      double precision dipole_scf_t ,dipole_cor_t
      double precision cdipole
      double precision debye,buck,octo
      parameter (octo=0.711688d0)
      parameter (buck=1.344911d0)
      parameter (debye=2.54176568d0)
      logical nodezero
      integer dummy
      integer sym_abelian_axis
      external sym_abelian_axis
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      irrep_y=0
      call tce_diis_init()
      do axis=1,3
        irrep_d=sym_abelian_axis(geom,axis)
        call tce_filename('d0',filename)
        call createfile(filename,d_d0,1) ! size_d0 = 1
        call ccsd_density1(d_d1(axis),d_d0,d_t1,d_t2,d_lambda1,
     1    d_lambda2,k_d1_offset(axis),k_d0_offset,k_t1_offset,
     2    k_t2_offset,k_l1_offset,k_l2_offset)
        call reconcilefile(d_d0,1)
        call get_block(d_d0,cdipole,1,0)
        dipole_cor(axis)=dipole_scf(axis)+cdipole
        call deletefile(d_d0)
      enddo
#if defined(NEWDIPOLE)
      dipole_scf_t = dsqrt( dipole_scf(1)*dipole_scf(1) +
     1                      dipole_scf(2)*dipole_scf(2) +
     2                      dipole_scf(3)*dipole_scf(3) )
      if (nodezero) write(LuOut,9430) "SCF",
     1              dipole_scf(1),dipole_scf(1)*debye,
     2              dipole_scf(2),dipole_scf(2)*debye,
     3              dipole_scf(3),dipole_scf(3)*debye,
     4              dipole_scf_t ,dipole_scf_t*debye
#endif
      dipole_cor_t = dsqrt( dipole_cor(1)*dipole_cor(1) +
     1                      dipole_cor(2)*dipole_cor(2) +
     2                      dipole_cor(3)*dipole_cor(3) )
      if (nodezero) write(LuOut,9430) "CCSD",
     1              dipole_cor(1),dipole_cor(1)*debye,
     2              dipole_cor(2),dipole_cor(2)*debye,
     3              dipole_cor(3),dipole_cor(3)*debye,
     4              dipole_cor_t ,dipole_cor_t*debye
#if defined(NEWDIPOLE)
      if (multipole.ge.2) then
        do axis=1,6
          irrep_d=sym_abelian_axis(geom,axis)
          call tce_filename('q0',filename)
          call createfile(filename,d_q0,1)
          call ccsd_density1(d_q1(axis),d_q0,d_t1,d_t2,d_lambda1,
     1      d_lambda2,k_q1_offset(axis),k_q0_offset,k_t1_offset,
     2      k_t2_offset,k_l1_offset,k_l2_offset)
          call reconcilefile(d_q0,1)
          call get_block(d_q0,cquadrupole,1,0)
          quadrupole_cor(axis)=quadrupole_scf(axis)-cquadrupole
          call deletefile(d_q0)
        enddo
        if (nodezero) write(LuOut,9432) "SCF",
     1                quadrupole_scf(1),quadrupole_scf(1)*buck,
     2                quadrupole_scf(2),quadrupole_scf(2)*buck,
     3                quadrupole_scf(3),quadrupole_scf(3)*buck,
     4                quadrupole_scf(4),quadrupole_scf(4)*buck,
     5                quadrupole_scf(5),quadrupole_scf(5)*buck,
     6                quadrupole_scf(6),quadrupole_scf(6)*buck
        if (nodezero) write(LuOut,9432) "CCSD",
     1                quadrupole_cor(1),quadrupole_cor(1)*buck,
     2                quadrupole_cor(2),quadrupole_cor(2)*buck,
     3                quadrupole_cor(3),quadrupole_cor(3)*buck,
     4                quadrupole_cor(4),quadrupole_cor(4)*buck,
     5                quadrupole_cor(5),quadrupole_cor(5)*buck,
     6                quadrupole_cor(6),quadrupole_cor(6)*buck
      endif
      if (multipole.ge.3) then
        do axis=1,10
          irrep_d=sym_abelian_axis(geom,axis)
          call tce_filename('o0',filename)
          call createfile(filename,d_o0,1)
          call ccsd_density1(d_o1(axis),d_o0,d_t1,d_t2,d_lambda1,
     1      d_lambda2,k_o1_offset(axis),k_o0_offset,k_t1_offset,
     2      k_t2_offset,k_l1_offset,k_l2_offset)
          call reconcilefile(d_o0,1)
          call get_block(d_o0,coctupole,1,0)
          octupole_cor(axis)=octupole_scf(axis)+coctupole
          call deletefile(d_o0)
        enddo
        if (nodezero) write(LuOut,9433) "SCF",
     1                octupole_scf(1), octupole_scf(1)*octo,
     2                octupole_scf(7), octupole_scf(7)*octo,
     3                octupole_scf(10),octupole_scf(10)*octo,
     4                octupole_scf(2), octupole_scf(2)*octo,
     5                octupole_scf(3), octupole_scf(3)*octo,
     2                octupole_scf(4), octupole_scf(4)*octo,
     3                octupole_scf(8), octupole_scf(8)*octo,
     4                octupole_scf(6), octupole_scf(6)*octo,
     5                octupole_scf(9), octupole_scf(9)*octo,
     6                octupole_scf(5), octupole_scf(5)*octo
        if (nodezero) write(LuOut,9433) "CCSD",
     1                octupole_cor(1), octupole_cor(1)*octo,
     2                octupole_cor(7), octupole_cor(7)*octo,
     3                octupole_cor(10),octupole_cor(10)*octo,
     4                octupole_cor(2), octupole_cor(2)*octo,
     5                octupole_cor(3), octupole_cor(3)*octo,
     2                octupole_cor(4), octupole_cor(4)*octo,
     3                octupole_cor(8), octupole_cor(8)*octo,
     4                octupole_cor(6), octupole_cor(6)*octo,
     5                octupole_cor(9), octupole_cor(9)*octo,
     6                octupole_cor(5), octupole_cor(5)*octo
      endif
#endif
      return
 9430 format(/,1x,A,' dipole moments / hartree & Debye',/,
     1  1x,'------------------------------------',/
     2  1x,'  X   ',2f15.7,/
     3  1x,'  Y   ',2f15.7,/
     4  1x,'  Z   ',2f15.7,/
     5  1x,'Total ',2f15.7,/
     6  1x,'------------------------------------')
 9432 format(/,1x,A,' quadrupole moments / hartree & Buckingham',/,
     1  1x,'--------------------------------',/
     2  1x,'XX ',2f15.7,/
     3  1x,'XY ',2f15.7,/
     4  1x,'XZ ',2f15.7,/
     5  1x,'YY ',2f15.7,/
     6  1x,'YZ ',2f15.7,/
     7  1x,'ZZ ',2f15.7,/
     8  1x,'--------------------------------')
c
c     Ordering in mptval (from defNxyz):
c     xxx, xxy, xxz, yyx, xyz, zzx, yyy, zzy, zzz
 9433 format(/,1x,A,' octupole moments / hartree & Octos',/,
     1  1x,'--------------------------------',/
     2  1x,'XXX ',2f15.7,/
     3  1x,'YYY ',2f15.7,/
     4  1x,'ZZZ ',2f15.7,/
     5  1x,'XXY ',2f15.7,/
     6  1x,'XXZ ',2f15.7,/
     7  1x,'YYX ',2f15.7,/
     8  1x,'YYZ ',2f15.7,/
     9  1x,'ZZX ',2f15.7,/
     1  1x,'ZZY ',2f15.7,/
     2  1x,'XYZ ',2f15.7,/
     3  1x,'--------------------------------')
      end
