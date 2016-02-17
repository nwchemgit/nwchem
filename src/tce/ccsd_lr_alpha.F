      subroutine ccsd_lr_alpha(rtdb,d_d0,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2           k_d0_offset,k_a0_offset,k_f1_offset,k_v2_offset,
     3           k_d1_offset,k_t1_offset,k_t2_offset,
     4           k_l1_offset,k_l2_offset,k_tr1_offset,k_tr2_offset,
     7           size_tr1,size_tr2,dipole_scf,dipole_cor)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer rtdb
      integer i,j,dummy,axis
      integer omegacount,omegasign,dynaxis,dynfreq
      integer irrep_g
      parameter (irrep_g=0)
      integer d_d0,d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_lambda1,d_lambda2
      integer d_tr1(6),d_tr2(6)
      integer k_d0_offset,k_a0_offset
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision cpu, wall
      double precision omega,cdipole
      double precision dipole_scf(3),dipole_cor(3)
      double precision debye
      parameter (debye=2.54176568d0)
      logical nodezero
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
      do omegacount=1,anumfreq
        omega = afreq(omegacount)
        if (omega.eq.(0.0d0)) then
          dynfreq = 1
        else
          dynfreq = 2
        endif
        do axis = 1, 3
          if (respaxis(axis)) then
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
            irrep_o=irrep_d
            irrep_x=irrep_d
c
            if (skip_tr) then
              if (nodezero) write(LuOut,9700)
            else
              do omegasign=1,dynfreq
                if (omegasign.eq.1) then
                  omega = (1.0d0)*afreq(omegacount)
                  dynaxis = 0
                elseif (omegasign.eq.2) then
                  omega = (-1.0d0)*afreq(omegacount)
                  dynaxis = 3
                endif
c
                if (nodezero) write(LuOut,9431) omega
c
                if ((omegacount.eq.1).and.guess_lr_real) then
                  call ccsd_lr_guess(d_f1,d_v2,d_d1,d_t1,d_t2,
     1                 d_tr1,d_tr2,k_f1_offset,k_v2_offset,
     2                 k_d1_offset,k_t1_offset,k_t2_offset,
     3                 k_tr1_offset,k_tr2_offset,
     4                 size_tr1,size_tr2,axis,dynaxis)
                endif
c
                call ccsd_lr_iter(axis,dynaxis,omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     4               size_tr1,size_tr2)
c
              enddo ! omegasign loop
            endif ! skip_tr
c
            if (omega.eq.(0.0d0)) then
              irrep_tr=irrep_d
              call tce_filename('d0',filename)
              call createfile(filename,d_d0,1)
              call tce_zero(d_d0,1)
              call dip_r(d_d0,d_d1(axis),d_t1,
     1             d_tr1(axis),d_tr2(axis),d_v2,
     2             k_d0_offset,k_d1_offset(axis),k_t1_offset,
     3             k_tr1_offset(axis),k_tr2_offset(axis),k_v2_offset)
              call reconcilefile(d_d0,1)
              cdipole=0.0d0
              call get_block(d_d0,cdipole,1,0)
              dipole_cor(axis)=dipole_scf(axis)+cdipole
              call deletefile(d_d0)
            endif
c
          endif ! respaxis(axis)
        enddo ! axis loop
c
        if (nodezero) write(LuOut,9430) "CCSD-LR",
     1                dipole_cor(1),dipole_cor(1)*debye,
     2                dipole_cor(2),dipole_cor(2)*debye,
     3                dipole_cor(3),dipole_cor(3)*debye
c
        if (nopolar) then
          if (nodezero) write(LuOut,9701)
        else
          call ccsd_lr_alpha_eval(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2         k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     6         k_tr1_offset,k_tr2_offset)
        endif
c
      enddo ! omegacount loop
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9420 format(1x,i4,f25.13,2f8.1)
 9430 format(/,1x,A,' dipole moments / hartree & Debye',/,
     1  1x,'--------------------------------',/
     2  1x,'X ',2f15.7,/
     3  1x,'Y ',2f15.7,/
     4  1x,'Z ',2f15.7,/
     5  1x,'--------------------------------')
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
 9700 format(1x,'Skipping CCSD-LR iterations')
 9701 format(1x,'Skipping evaluation of polarizability by request')
      return
      end


      subroutine ccsd_lr_iter(axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     4           size_tr1,size_tr2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis,iter
      integer omegacount,omegasign,dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(6),d_tr2(6),d_rr1(3),d_rr2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision residual,rr1,rr2
      double precision ddotfile
      external ddotfile
      character*255 filename
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1)) write(LuOut,9400) "CCSD-LR"
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
c        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call daxfile(1,(-1.0d0)*omega,d_tr1(axis+dynaxis),
     1       d_rr1(axis),size_tr1(axis))
        call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),k_f1_offset,
     2       k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
c
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
c        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call daxfile(1,(-1.0d0)*omega,d_tr2(axis+dynaxis),
     1       d_rr2(axis),size_tr2(axis))
        call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
c
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
        residual = max(rr1,rr2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) then
            write(LuOut,9410)
          endif
          if (ampnorms) then
            call tce_residual_tr1(d_tr1(axis+dynaxis),
     1                            k_tr1_offset(axis),rr1)
            call tce_residual_tr2(d_tr2(axis+dynaxis),
     1                            k_tr2_offset(axis),rr2)
            if (nodezero) then
              write(LuOut,9082) "T(1) singles",rr1
              write(LuOut,9082) "T(1) doubles",rr2
            endif
          endif
          call deletefile(d_rr2(axis))
          call deletefile(d_rr1(axis))
          call tce_diis_tidy()
          if (save_tr(1)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1(1) now...'
            endif
            call tr1_restart_save(d_tr1(axis+dynaxis),
     1           k_tr1_offset(axis),size_tr1(axis),
     2           axis+dynaxis,handle_tr1(axis),irrep_x)
          endif
          if (save_tr(2)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2(1) now...'
            endif
            call tr2_restart_save(d_tr2(axis+dynaxis),
     1           k_tr2_offset(axis),size_tr2(axis),
     2           axis+dynaxis,handle_tr2(axis),irrep_x)
          endif
          return
        endif
        if (save_tr(1).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving T1(1) now...'
          endif
          call tr1_restart_save(d_tr1(axis+dynaxis),
     1         k_tr1_offset(axis),size_tr1(axis),
     2         axis+dynaxis,handle_tr1(axis),irrep_x)
        endif
        if (save_tr(2).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving T2(1) now...'
          endif
          call tr2_restart_save(d_tr2(axis+dynaxis),
     1         k_tr2_offset(axis),size_tr2(axis),
     2         axis+dynaxis,handle_tr2(axis),irrep_x)
        endif
        call tce_diis2(.false.,iter,.true.,.true.,.false.,.false.,
     1       d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2       size_tr1(axis),
     3       d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4       size_tr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,omega)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_lr_iter: maxiter exceeded',iter,CALC_ERR)
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9400 format(/,1x,A,' iterations',/,
     1  1x,'---------------------------------------------',/
     2  1x,'Iter          Residuum            Cpu    Wall',/
     3  1x,'---------------------------------------------')
 9410 format(
     1  1x,'---------------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f25.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end


      subroutine ccsd_lr_alpha_eval(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     6           k_tr1_offset,k_tr2_offset)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
c
      integer rtdb
      integer i,j,dummy,axis
      integer dynfreq,dynaxis,axisA,axisB,omegasign
      integer d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_lambda1,d_lambda2
      integer d_tr1(6),d_tr2(6)
      integer k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer irrep_g
      parameter (irrep_g=0)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision alpha1,alpha2
      double precision alpha(3,3)       ! Dipole polarizability tensor
      double precision alphacopy(3,3)   ! Dipole polarizability tensor copy
      double precision alphaiso         ! Isotropic dipole polarizability
      double precision alphaani         ! Anisotropic dipole polarizability
      double precision alphaevr(3)      ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)      ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3)     ! dummy dgeev eigenvectors
      double precision auXnm            ! Conversion factor from a.u. (frequency) to nm (wavelength)
      parameter (auXnm=45.563353d0)
      double precision au2ang           ! Conversion factor from bohr to Angstrom
      parameter (au2ang=5.29177249d-01)
      double precision au2ang3          ! Conversion factor from bohr^3 to Angstrom^3
      parameter (au2ang3=au2ang*au2ang*au2ang)
      double precision tol_omega
      parameter (tol_omega=1.0d-50)
      double precision ddotfile
      external ddotfile
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c CCSD-LR evaluation step
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
      do axisA = 1, 3
        do axisB = 1, axisA
c
c         initialize alpha here to ensure all elements are sensible
c         as the body of the if-statement is not always executed
c
          alpha(axisA,axisB)=0.0d0
c
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,1)
            call tce_zero(d_a0,1)
c
#ifdef DEBUG_PRINT
            if (nodezero) then
              write(LuOut,*) "axisA = ",axisA
              write(LuOut,*) "axisB = ",axisB
              call util_flush(LuOut)
            endif
#endif
c
            do omegasign=1, 2
c
              if ((omegasign.eq.2).and.(omega.ne.0.0d0)) then
                dynaxis = 3
              else
                dynaxis = 0
              endif
c
              irrep_d=irrep_a
              irrep_tr=irrep_b
              call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1             d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     2             d_lambda1,d_lambda2,
     3             k_d1_offset(axisA),k_a0_offset,
     4             k_t1_offset,k_t2_offset,
     5             k_tr1_offset(axisB),k_tr2_offset(axisB),
     6             k_l1_offset,k_l2_offset)
c
              alpha1=0.0d0
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,1)
              if (axisA.eq.axisB) then
                alpha2=alpha1
              else
                irrep_d=irrep_b
                irrep_tr=irrep_a
                call alpha_1(d_d1(axisB),d_a0,d_t1,d_t2,
     1               d_tr1(axisA+dynaxis),d_tr2(axisA+dynaxis),
     2               d_lambda1,d_lambda2,
     3               k_d1_offset(axisB),k_a0_offset,
     4               k_t1_offset,k_t2_offset,
     5               k_tr1_offset(axisA),k_tr2_offset(axisA),
     6               k_l1_offset,k_l2_offset)
c
                alpha2=0.0d0
                call reconcilefile(d_a0,1)
                call get_block(d_a0,alpha2,1,0)
                call tce_zero(d_a0,1)
              endif ! axisA.eq.axisB
c
#ifdef DEBUG_PRINT
              if (nodezero) then
                write(LuOut,*) "alpha1 = ",alpha1
                write(LuOut,*) "alpha2 = ",alpha2
                call util_flush(LuOut)
              endif
#endif
c
              alpha(axisA,axisB)=alpha(axisA,axisB)+alpha1+alpha2
c
            enddo ! omegasign loop
c
            alpha1=0.0d0
            alpha2=0.0d0
c
            irrep_tra=irrep_a
            irrep_trb=irrep_b
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1           d_tr1(axisA),d_tr2(axisA),
     2           d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
c
            alpha1=0.0d0
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,1)
c
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1           d_tr1(axisA+dynaxis),d_tr2(axisA+dynaxis),
     2           d_tr1(axisB),d_tr2(axisB),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
c
            alpha2=0.0d0
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha2,1,0)
            call deletefile(d_a0)
            alpha(axisA,axisB)=alpha(axisA,axisB)+alpha1+alpha2
            alpha(axisA,axisB)=-0.5d0*alpha(axisA,axisB)
c
#ifdef DEBUG_PRINT
            if (nodezero) then
              write(LuOut,*) "alpha1 = ",alpha1
              write(LuOut,*) "alpha2 = ",alpha2
              write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
              call util_flush(LuOut)
            endif
#endif
c
          endif ! respaxis(axis)
        enddo ! axisB loop
      enddo ! axisA loop
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
c
c     Take a copy of alpha for the diagonalization, remember though
c     that only a triangle of alpha was calculated, the rest contains
c     junk.
c
      do i = 1, 3
        do j = 1, i
          alphacopy(i,j)=alpha(i,j)
          alphacopy(j,i)=alpha(i,j)
        enddo
      enddo
c
      call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
c
      alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
      alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1         + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2         + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
      alphaani = dsqrt(0.5d0*alphaani)
c
      if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1    .and.nodezero) write(LuOut,9911)
c
      if (nodezero) then
        if (dabs(omega).gt.tol_omega) then
          write(LuOut,9434) "CCSD Linear Response",
     1      dabs(omega),dabs(auXnm/omega),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,alphaani,au2ang3*alphaani
        else
          write(LuOut,9436) "CCSD Linear Response",
     1      dabs(omega),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,alphaani,au2ang3*alphaani
        endif
      endif
c
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
c
      if (.not. rtdb_put(rtdb, 'raman:alfare ', mt_dbl, 9, alphacopy))
     &  call errquit('ccsd_lr_alpha_eval: rtdb_put fail', 0, RTDB_ERR)
c     zero this out and call it imaginary part, since that is not computed
      do i = 1, 3
        do j = 1, 3
          alphacopy(i,j)=0
        enddo
      enddo
      if (.not. rtdb_put(rtdb, 'raman:alfaim ', mt_dbl, 9, alphacopy))
     &  call errquit('ccsd_lr_alpha_eval: rtdb_put fail', 1, RTDB_ERR)

 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9400 format(/,1x,A,' iterations',/,
     1  1x,'---------------------------------------------',/
     2  1x,'Iter          Residuum            Cpu    Wall',/
     3  1x,'---------------------------------------------')
 9410 format(
     1  1x,'---------------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f25.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9434 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength = ',f15.7,' / nm',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     2  1x,'                    atomic units (bohr^3)       '
     2    ,'       |                   angstroms^3           ',/
     2  1x,'                 X              Y              Z',
     2  1x,'      |             X              Y              Z',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     4  1x,'X      ',3f15.7,3x,'|',3f15.7,/
     5  1x,'Y      ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Z      ',3f15.7,3x,'|',3f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     6  1x,'Eigs = ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Isotropic   = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     6  1x,'Anisotropic = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------')
 9436 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength =        Infinity / nm',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     2  1x,'                    atomic units (bohr^3)       '
     2    ,'       |                   angstroms^3           ',/
     2  1x,'                 X              Y              Z',
     2  1x,'      |             X              Y              Z',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     4  1x,'X      ',3f15.7,3x,'|',3f15.7,/
     5  1x,'Y      ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Z      ',3f15.7,3x,'|',3f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     6  1x,'Eigs = ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Isotropic   = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     6  1x,'Anisotropic = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------')
 9435 format(/,1x,A,' C6 coefficients ',/
     1  1x,'--------------------------------',/
     2  1x,'C6(XX)  ',f15.7,/
     3  1x,'C6(YY)  ',f15.7,/
     4  1x,'C6(ZZ)  ',f15.7,/
     5  1x,'C6(XY)  ',f15.7,/
     6  1x,'C6(XZ)  ',f15.7,/
     7  1x,'C6(YZ)  ',f15.7,/
     8  1x,'C6(AVG) ',f15.7,/
     9  1x,'C6(ANI) ',f15.7,/
     1  1x,'--------------------------------')
 9911 format(/,1x,'Warning: you have not solved ',
     1            'the response equations for all axes.  ',
     2            'Please analyze the results carefully as ',
     3            'the average and anisotropic polarizabilities ',
     4            'are surely wrong.',/)
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end




      subroutine ccsd_lr_guess(d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     1           k_f1_offset,k_v2_offset,k_d1_offset,
     2           k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     3           size_tr1,size_tr2,axis,dynaxis)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
c
      integer axis,dynaxis
      integer d_f1,d_v2,d_d1(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer d_t1
      integer k_t1_offset
      integer d_t2
      integer k_t2_offset
      integer d_tr1(6)
      integer k_tr1_offset(3)
      integer size_tr1(3)
      integer d_tr2(6)
      integer k_tr2_offset(3)
      integer size_tr2(3)
      integer d_rr1(3)
      integer d_rr2(3)
      double precision omega
      double precision cpu, wall
      logical nodezero
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
      if (nodezero) write(6,9121) 'Initial guess x = b/Adiag'
c
c     0. Create dummies
c
      call tce_filename(rr1filename(axis),filename)
      call createfile(filename,d_rr1(axis),size_tr1(axis))
      call tce_zero(d_rr1(axis),size_tr1(axis))
c
      call tce_filename(rr2filename(axis),filename)
      call createfile(filename,d_rr2(axis),size_tr2(axis))
      call tce_zero(d_rr2(axis),size_tr2(axis))
c
c     1. Form b
c
      call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1     k_tr1_offset(axis),k_d1_offset(axis),
     2     k_t1_offset,k_t2_offset)
      call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1     k_tr2_offset(axis),k_d1_offset(axis),
     2     k_t1_offset,k_t2_offset)
c
c     3. Hit with preconditioner
c
      call tce_jacobi_tr1(d_rr1(axis),d_tr1(axis+dynaxis),
     1                    k_tr1_offset(axis),0.0d0)
      call tce_jacobi_tr2(d_rr2(axis),d_tr2(axis+dynaxis),
     1                    k_tr2_offset(axis),0.0d0)
c
c     4. Delete dummies
c
      call deletefile(d_rr1(axis))
      call deletefile(d_rr2(axis))
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9420 format(1x,i4,f25.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end

#ifdef DENOM_EXPER
c====================================================================
c             DENOMINATOR EXPERIMENT
c====================================================================
          call tce_filename('denom1',filename)
          call createfile(filename,d_denom1,size_t1)
          call tce_zero(d_denom1,size_t1)
c
          call tce_filename('denom2',filename)
          call createfile(filename,d_denom2,size_t2)
          call tce_zero(d_denom2,size_t2)
c
          call tce_make_denom1(d_denom1,k_t1_offset,irrep_t,
     1                         1,0.0d0,0.0d0) ! denom_power, omega, shift
          call tce_make_denom2(d_denom2,k_t2_offset,irrep_t,
     1                         1,0.0d0,0.0d0) ! denom_power, omega, shift
c
          call tce_print_x1(d_denom1,k_t1_offset,1.0d-6,irrep_t)
          call tce_print_x2(d_denom2,k_t2_offset,1.0d-6,irrep_t)
c
          call deletefile(d_denom1)
          call deletefile(d_denom2)
c
          do axis = 1, 3
          do i = 1, 2
          if (respaxis(axis)) then
c
            if(nodezero) write(6,*) '================================='
            if(nodezero) write(6,*) 'denom_power = ',i
c
            irrep_d=sym_abelian_axis(geom,axis)
            call sym_irrepname(geom,irrep_d+1,irrepname)
            if (nodezero.and.util_print('mod1',print_default)) then
              write(LuOut,*)
              write(LuOut,9440) axisname(axis),irrepname
            endif
c
            call tce_filename('denom1',filename)
            call createfile(filename,d_denom1,size_tr1(axis))
            call tce_zero(d_denom1,size_tr1(axis))
c
            call tce_filename('denom2',filename)
            call createfile(filename,d_denom2,size_tr2(axis))
            call tce_zero(d_denom2,size_tr2(axis))
c
            call tce_make_denom1(d_denom1,k_tr1_offset(axis),irrep_d,
     1                             i,0.0d0,0.0d0) ! denom_power, omega, shift
            call tce_make_denom2(d_denom2,k_tr2_offset(axis),irrep_d,
     1                             i,0.0d0,0.0d0) ! denom_power, omega, shift
c
            call tce_print_x1(d_denom1,k_tr1_offset(axis),1.0d-6,
     1                        irrep_d)
            call tce_print_x2(d_denom2,k_tr2_offset(axis),1.0d-6,
     2                        irrep_d)
c
            call deletefile(d_denom1)
            call deletefile(d_denom2)
c
          endif
          enddo
          enddo
c
          call errquit('tce_energy: MANUAL STOP',0,CALC_ERR)
c====================================================================
#endif
c
c $Id$
