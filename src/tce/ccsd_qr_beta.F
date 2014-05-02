      subroutine ccsd_qr_beta(rtdb,d_a0,d_b0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_y1,d_y2,d_yr1,d_yr2,d_tr1,d_tr2,
     2           k_a0_offset,k_b0_offset,
     3           k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     5           k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset,
     6           size_yr1,size_yr2)
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
      integer i,j,dummy,axis,dynaxis,dynfreq
      integer omegacount,omegasign
      integer irrep_g
      parameter (irrep_g=0)
      integer d_a0,d_b0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_y1,d_y2
      integer d_tr1(9),d_tr2(9),d_yr1(6),d_yr2(6)
      integer k_a0_offset,k_b0_offset
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer k_yr1_offset(3),k_yr2_offset(3)
      integer size_yr1(3),size_yr2(3)
      integer d_rr1(3),d_rr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      logical nodezero,guess
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
c
      nodezero=(ga_nodeid().eq.0)
c
#ifdef DEBUG_PRINT
      if (nodezero) write(6,*) 'top of ccsd_qr_beta'
      call util_flush(6)
#endif
c
      do omegacount=1,bnumfreq
        omega = bfreq(omegacount)
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
                  omega = (1.0d0)*bfreq(omegacount)
                  dynaxis = 0
                elseif (omegasign.eq.2) then
                  omega = (-1.0d0)*bfreq(omegacount)
                  dynaxis = 3
                endif
c
                if (nodezero) write(LuOut,9431) omega
c
                guess = ((omegacount.eq.1).and.guess_qr_real)
                call ccsd_qr_iter(guess,axis,dynaxis,omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_y1,d_y2,
     2               d_tr1,d_tr2,d_yr1,d_yr2,k_f1_offset,k_v2_offset,
     3               k_d1_offset,k_t1_offset,k_t2_offset,
     4               k_y1_offset,k_y2_offset,k_tr1_offset,k_tr2_offset,
     5               k_yr1_offset,k_yr2_offset,size_yr1,size_yr2)
c
              enddo ! omegasign loop
            endif ! skip_tr
          endif ! respaxis(axis)
        enddo ! axis loop
c
#ifdef DEBUG_PRINT
        if (nodezero) write(6,*) 'calling ccsd_qr_alpha_eval'
        call util_flush(6)
#endif
c
        call ccsd_qr_alpha_eval(omega,d_a0,d_f1,d_v2,d_d1,
     1       d_t1,d_t2,d_y1,d_y2,d_tr1,d_tr2,d_yr1,d_yr2,
     2       k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4       k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     5       k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset)
c
#ifdef DEBUG_PRINT
        if (nodezero) write(6,*) 'calling ccsd_qr_beta_eval'
        call util_flush(6)
#endif
c
        if (nopolar) then
          if (nodezero) write(LuOut,9701)
        else
          call ccsd_qr_beta_eval(omega,d_b0,d_f1,d_v2,d_d1,
     1         d_t1,d_t2,d_y1,d_y2,d_tr1,d_tr2,d_yr1,d_yr2,
     2         k_b0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4         k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     5         k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset)
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
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
 9700 format(1x,'Skipping CCSD-QR iterations')
 9701 format(1x,'Skipping evaluation of hyperpolarizability by request')
      return
      end


      subroutine ccsd_qr_iter(guess,axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_y1,d_y2,
     2           d_tr1,d_tr2,d_yr1,d_yr2,
     3           k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     5           k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset,
     6           size_yr1,size_yr2)
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
      integer i,j,dummy,axis,dynaxis,iter
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_y1,d_y2,d_b1(3),d_b2(3)
      integer d_tr1(9),d_tr2(9),d_yr1(6),d_yr2(6)
      integer d_yrr1(3),d_yrr2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_y1_offset,k_y2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer k_yr1_offset(3),k_yr2_offset(3)
      integer size_yr1(3),size_yr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision residual
      double precision rr1            ! Residual tr1
      double precision rr2            ! Residual tr2
      double precision ddotfile
      external ddotfile
      logical nodezero,guess
      character*4 irrepname
      character*255 filename
      character*5 yrr1filename(3) ! File name stub
      data yrr1filename/'yrr1x','yrr1y','yrr1z'/
      character*5 yrr2filename(3) ! File name stub
      data yrr2filename/'yrr2x','yrr2y','yrr2z'/
      character*3 b1filename(3) ! File name stub
      data b1filename/'b1x','b1y','b1z'/
      character*3 b2filename(3) ! File name stub
      data b2filename/'b2x','b2y','b2z'/
c
      nodezero=(ga_nodeid().eq.0)
c
c     form and store constant RHS of linear equations
c
      if (nodezero) write(LuOut,9120) 'Forming RHS'
      cpu=-util_cpusec()
      wall=-util_wallsec()
c
      call tce_filename(b1filename(axis),filename)
      call createfile(filename,d_b1(axis),size_yr1(axis))
      call tce_zero(d_b1(axis),size_yr1(axis))
c
      call tce_filename(b2filename(axis),filename)
      call createfile(filename,d_b2(axis),size_yr2(axis))
      call tce_zero(d_b2(axis),size_yr2(axis))
c
      if (save_lr(1)) then
        if(nodezero) then
          write(LuOut,*) 'Saving QR-B1 now...'
        endif
        call yr1_restart_save(d_b1(axis),k_yr1_offset(axis),
     1       size_yr1(axis),axis,handle_qb1(axis),irrep_d)
      endif
      if (save_lr(2)) then
        if(nodezero) then
          write(LuOut,*) 'Saving QR-B2 now...'
        endif
        call yr2_restart_save(d_b2(axis),k_yr2_offset(axis),
     1       size_yr2(axis),axis,handle_qb2(axis),irrep_d)
      endif
c
      irrep_y=irrep_g ! using Lambda(0) amplitudes
      call ccsd_y_o1(d_b1(axis),d_d1(axis),d_t1,d_t2,
     1     d_y1,d_y2,k_yr1_offset(axis),
     2     k_d1_offset(axis),k_t1_offset,k_t2_offset,
     3     k_y1_offset,k_y2_offset)
      call ccsd_y_tr1(d_f1,d_b1(axis),d_t1,d_t2,d_v2,
     1     d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),d_y1,
     2     d_y2,k_f1_offset,k_yr1_offset(axis),
     3     k_t1_offset,k_t2_offset,k_v2_offset,
     4     k_tr1_offset(axis),k_tr2_offset(axis),
     5     k_y1_offset,k_y2_offset)
      call ccsd_y_o2(d_b2(axis),d_d1(axis),d_t1,
     1     d_y1,d_y2,k_yr2_offset(axis),
     2     k_d1_offset(axis),k_t1_offset,
     3     k_y1_offset,k_y2_offset)
      call ccsd_y_tr2(d_f1,d_b2(axis),d_t1,d_v2,
     1     d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),d_y1,
     2     d_y2,k_f1_offset,k_yr2_offset(axis),
     3     k_t1_offset,k_v2_offset,
     4     k_tr1_offset(axis),k_tr2_offset(axis),
     5     k_y1_offset,k_y2_offset)
c
      call reconcilefile(d_b1(axis),size_yr1(axis))
      call reconcilefile(d_b2(axis),size_yr2(axis))
c
      if (guess) then
        if (nodezero) write(6,9121) 'Initial guess y = b/Adiag'
c
        call tce_zero(d_yr1(axis),size_yr1(axis))
        call tce_zero(d_yr2(axis),size_yr2(axis))
c
        irrep_y = irrep_d
        call tce_jacobi_lr1(d_b1(axis),d_yr1(axis),k_yr1_offset(axis),
     1                      0.0d0,0.0d0)
        call tce_jacobi_lr2(d_b2(axis),d_yr2(axis),k_yr2_offset(axis),
     1                      0.0d0,0.0d0)
      endif
c
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu,wall
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1)) write(LuOut,9400) "CCSD-YR"
c
        call tce_filename(yrr1filename(axis),filename)
        call createfile(filename,d_yrr1(axis),size_yr1(axis))
        call copyfile(d_b1(axis),d_yrr1(axis),size_yr1(axis))
        call daxpyfile(1,(1.0d0)*omega,d_yr1(axis+dynaxis),
     1       d_yrr1(axis),size_yr1(axis))
c
        irrep_y=irrep_d ! using Lambda(1) amplitudes
        call eomccsd_y1(d_f1,d_yrr1(axis),d_t1,d_t2,d_v2,
     1       d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),k_f1_offset,
     2       k_yr1_offset(axis),k_t1_offset,k_t2_offset,
     3       k_v2_offset,k_yr1_offset(axis),k_yr2_offset(axis))
c
        call tce_filename(yrr2filename(axis),filename)
        call createfile(filename,d_yrr2(axis),size_yr2(axis))
        call copyfile(d_b2(axis),d_yrr2(axis),size_yr2(axis))
        call daxpyfile(1,(1.0d0)*omega,d_yr2(axis+dynaxis),
     1       d_yrr2(axis),size_yr2(axis))
c
        irrep_y=irrep_d ! using Lambda(1) amplitudes
        call eomccsd_y2(d_f1,d_yrr2(axis),d_t1,d_t2,d_v2,
     1       d_yr1(axis+dynaxis),d_yr2(axis+dynaxis),k_f1_offset,
     2       k_yr2_offset(axis),k_t1_offset,k_t2_offset,
     3       k_v2_offset,k_yr1_offset(axis),k_yr2_offset(axis))
c
        call reconcilefile(d_yrr1(axis),size_yr1(axis))
        call reconcilefile(d_yrr2(axis),size_yr2(axis))
c
        rr1 = ddotfile(d_yrr1(axis),d_yrr1(axis),size_yr1(axis))
        rr2 = ddotfile(d_yrr2(axis),d_yrr2(axis),size_yr2(axis))
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
            rr1 = ddotfile(d_yr1(axis+dynaxis),d_yr1(axis+dynaxis),
     1                     size_yr1(axis))
            rr1 = ddotfile(d_yr2(axis+dynaxis),d_yr2(axis+dynaxis),
     1                     size_yr2(axis))
            if (nodezero) then
              write(LuOut,9082) "Lambda(1) singles",rr1
              write(LuOut,9082) "Lambda(1) doubles",rr2
            endif
          endif
c
          call deletefile(d_yrr2(axis))
          call deletefile(d_yrr1(axis))
          call tce_diis_tidy()
          if (save_lr(1)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L1(1) now...'
            endif
            call yr1_restart_save(d_yr1(axis),k_yr1_offset(axis),
     1           size_yr1(axis),axis,handle_yr1(axis),irrep_d)
          endif
          if (save_lr(2)) then
            if(nodezero) then
              write(LuOut,*) 'Saving L2(1) now...'
            endif
            call yr2_restart_save(d_yr2(axis),k_yr2_offset(axis),
     1           size_yr2(axis),axis,handle_yr2(axis),irrep_d)
          endif
          call deletefile(d_b2(axis))
          call deletefile(d_b1(axis))
          return
        endif
        if (save_lr(1).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving L1(1) now...'
          endif
          call yr1_restart_save(d_yr1(axis),k_yr1_offset(axis),
     1         size_yr1(axis),axis,handle_yr1(axis),irrep_d)
        endif
        if (save_lr(2).and.(mod(iter,save_interval).eq.0)) then
          if(nodezero) then
            write(LuOut,*) 'Saving L2(1) now...'
          endif
          call yr2_restart_save(d_yr2(axis),k_yr2_offset(axis),
     1   size_yr2(axis),axis,handle_yr2(axis),irrep_d)
        endif
c        write(LuOut,*) "before tce_diis2"
        irrep_y=irrep_d
        call tce_diis2(.true.,iter,.true.,.true.,.false.,.false.,
     1       d_yrr1(axis),d_yr1(axis+dynaxis),k_yr1_offset(axis),
     2       size_yr1(axis),
     3       d_yrr2(axis),d_yr2(axis+dynaxis),k_yr2_offset(axis),
     4       size_yr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,omega)
c       write(LuOut,*) "after tce_diis2"
        call deletefile(d_yrr2(axis))
        call deletefile(d_yrr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call deletefile(d_b2(axis))
      call deletefile(d_b1(axis))
      call errquit('ccsd_qr_iter: maxiter exceeded',iter,CALC_ERR)
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(1x,A)
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


      subroutine ccsd_qr_alpha_eval(omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_y1,d_y2,d_tr1,d_tr2,d_yr1,d_yr2,
     2           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     4           k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset)
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
      integer i,j,dummy,axis,axisA,axisB
      integer irrep_g
      parameter (irrep_g=0)
      integer d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_y1,d_y2
      integer k_a0_offset
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset
      integer d_tr1(9),d_tr2(9),d_yr1(6),d_yr2(6)
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer k_yr1_offset(3),k_yr2_offset(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision residual
      double precision alpha1,alpha2,alpha3,alpha4
      double precision au2ang   ! Conversion factor from bohr to Angstrom
      double precision au2ang3  ! Conversion factor from bohr^3 to Angstrom^3
      double precision auXnm    ! Conversion factor from a.u. (frequency) to nm (wavelength)
      double precision alpha(3,3)       ! Dipole polarizability tensor
      double precision alphacopy(3,3)   ! Dipole polarizability tensor copy
      double precision alphaiso         ! Isotropic dipole polarizability
      double precision alphaani         ! Anisotropic dipole polarizability
      double precision alphaevr(3)      ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)      ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3)     ! dummy dgeev eigenvectors
      double precision ddotfile
      external ddotfile
      parameter (auXnm=45.563353d0)
      parameter (au2ang=5.29177249d-01)
      parameter (au2ang3=au2ang*au2ang*au2ang)
      logical nodezero
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
#ifdef DEBUG_PRINT
      if (nodezero) write(6,*) 'top of ccsd_qr_alpha_eval'
      call util_flush(6)
#endif
c
c CCSD-QR evaluation step
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
      do axisA = 1, 3
        irrep_a=sym_abelian_axis(geom,axisA)
        do axisB = 1, axisA
          irrep_b=sym_abelian_axis(geom,axisB)
          alpha(axisA,axisB)=0.0d0
          if ( respaxis(axisA).and.respaxis(axisB) ) then
c
#ifdef DEBUG_PRINT
            if (nodezero) then
              write(LuOut,*) "axisA = ",axisA
              write(LuOut,*) "axisB = ",axisB
              call util_flush(6)
            endif
#endif
c
            alpha1=0.0d0
            alpha2=0.0d0
            alpha3=0.0d0
            alpha4=0.0d0
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,1)
            call tce_zero(d_a0,1)
c
            irrep_y=irrep_g
            irrep_d=irrep_a
            irrep_tr=irrep_b
c
            call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1           d_tr1(axisB),d_tr2(axisB),
     2           d_y1,d_y2,k_d1_offset(axisA),k_a0_offset,
     4           k_t1_offset,k_t2_offset,k_tr1_offset(axisB),
     5           k_tr2_offset(axisB),k_y1_offset,k_y2_offset)
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,1)
c
            if (omega.eq.0.0d0) then
              alpha2=alpha1
            else
              call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1             d_tr1(axisB+3),d_tr2(axisB+3),
     2             d_y1,d_y2,k_d1_offset(axisA),k_a0_offset,
     4             k_t1_offset,k_t2_offset,k_tr1_offset(axisB),
     5             k_tr2_offset(axisB),k_y1_offset,k_y2_offset)
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,1)
            endif
c
            irrep_o=irrep_a
            irrep_y=irrep_b
c
            call alpha_2_yr(d_a0,d_d1(axisA),d_t1,d_t2,
     1           d_yr1(axisB),d_yr2(axisB),
     2           k_a0_offset,k_d1_offset(axisA),k_t1_offset,
     3           k_t2_offset,k_yr1_offset(axisB),k_yr2_offset(axisB))
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha3,1,0)
            call tce_zero(d_a0,1)
c
            if (omega.eq.0.0d0) then
              alpha4=alpha3
            else
              call alpha_2_yr(d_a0,d_d1(axisA),d_t1,d_t2,
     1             d_yr1(axisB+3),d_yr2(axisB+3),
     2             k_a0_offset,k_d1_offset(axisA),k_t1_offset,
     3             k_t2_offset,k_yr1_offset(axisB),k_yr2_offset(axisB))
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha4,1,0)
              call tce_zero(d_a0,1)
            endif
c
            alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2+alpha3+alpha4)
            call deletefile(d_a0)
c
#ifdef DEBUG_PRINT
            if (nodezero) then
              write(LuOut,*) "alpha1 = ",alpha1
              write(LuOut,*) "alpha2 = ",alpha2
              write(LuOut,*) "alpha3 = ",alpha3
              write(LuOut,*) "alpha4 = ",alpha4
              write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
              call util_flush(LuOut)
            endif
#endif
c
          endif ! respaxis(axisA,axisB)
        enddo ! axisB loop
      enddo ! axisA loop
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
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
      if (omega.ne.0.0d0) then
        if (nodezero) write(LuOut,9434) "CCSD Lambda Response",
     1    omega,auXnm/omega,
     2    alpha(1,1),alpha(2,1),alpha(3,1),
     3    au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4    alpha(2,1),alpha(2,2),alpha(3,2),
     5    au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6    alpha(3,1),alpha(3,2),alpha(3,3),
     7    au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8    alphaevr(1),alphaevr(2),alphaevr(3),
     9    au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1    alphaiso,au2ang3*alphaiso,alphaani,au2ang3*alphaani
      else
        if (nodezero) write(LuOut,9435) "CCSD Lambda Response",
     1    omega,
     2    alpha(1,1),alpha(2,1),alpha(3,1),
     3    au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4    alpha(2,1),alpha(2,2),alpha(3,2),
     5    au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6    alpha(3,1),alpha(3,2),alpha(3,3),
     7    au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8    alphaevr(1),alphaevr(2),alphaevr(3),
     9    au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1    alphaiso,au2ang3*alphaiso,alphaani,au2ang3*alphaani
      endif
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
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
 9435 format(/,1x,A,' polarizability / au ',/
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
 9911 format(/,1x,'Warning: you have not solved ',
     1            'the response equations for all axes.  ',
     2            'Please analyze the results carefully as ',
     3            'the average and anisotropic polarizabilities ',
     4            'are surely wrong.',/)
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end


      subroutine ccsd_qr_beta_eval(omega,d_b0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_y1,d_y2,d_tr1,d_tr2,d_yr1,d_yr2,
     2           k_b0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     4           k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset)
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
      integer i,j,dummy,axis
      integer axisA,axisB,axisC
      integer irrep_g
      parameter (irrep_g=0)
      integer d_b0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_y1,d_y2
      integer k_b0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset
      integer d_tr1(9),d_tr2(9),d_yr1(6),d_yr2(6)
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer k_yr1_offset(3),k_yr2_offset(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision beta1,beta2,beta3,beta4
      double precision beta(3,3,3)      ! Dipole hyperpolarizability tensor
      double precision au2ang   ! Conversion factor from bohr to Angstrom
      double precision au2ang3  ! Conversion factor from bohr^3 to Angstrom^3
      double precision auXnm    ! Conversion factor from a.u. (frequency) to nm (wavelength)
      double precision ddotfile
      external ddotfile
      parameter (auXnm=45.563353d0)
      parameter (au2ang=5.29177249d-01)
      parameter (au2ang3=au2ang*au2ang*au2ang)
      logical nodezero
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
#ifdef DEBUG_PRINT
      if (nodezero) write(6,*) 'top of ccsd_qr_beta_eval'
      call util_flush(6)
#endif
c
      do axisC = 1, 3
        do axisB = 1, 3
          do axisA = 1, 3
            beta(axisA,axisB,axisC)=0.0d0
          enddo
        enddo
      enddo
c
c CCSD-QR evaluation step
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
      do axisA = 1, 3
        irrep_a=sym_abelian_axis(geom,axisA)
        do axisB = 1, 3
          irrep_b=sym_abelian_axis(geom,axisB)
          do axisC = 1, 3
            irrep_c=sym_abelian_axis(geom,axisC)
            if ( ( respaxis(axisA).and.respaxis(axisB) )
     &                          .and.respaxis(axisC) ) then
c
#ifdef DEBUG_PRINT
              if (nodezero) then
                write(LuOut,*) "axisA = ",axisA
                write(LuOut,*) "axisB = ",axisB
                write(LuOut,*) "axisC = ",axisC
                call util_flush(6)
              endif
#endif
c
              beta1=0.0d0
              beta2=0.0d0
              beta3=0.0d0
              beta4=0.0d0
              call tce_filename('b0',filename)
              call createfile(filename,d_b0,1)
c
              irrep_y=irrep_g
              irrep_tra=irrep_a
              irrep_trb=irrep_b
              irrep_trc=irrep_c
c              write(LuOut,*) "beta_1"
              call beta_1(d_b0,d_t1,d_tr1(axisA),d_tr2(axisA),
     1             d_tr1(axisB),d_tr2(axisB),d_tr1(axisC),d_tr2(axisC),
     2             d_v2,d_y1,d_y2,k_b0_offset,k_t1_offset,
     3             k_tr1_offset(axisA),k_tr2_offset(axisA),
     4             k_tr1_offset(axisB),k_tr2_offset(axisB),
     5             k_tr1_offset(axisC),k_tr2_offset(axisC),
     6             k_v2_offset,k_y1_offset,k_y2_offset)
              call reconcilefile(d_b0,1)
              call get_block(d_b0,beta1,1,0)
              call tce_zero(d_b0,1)
c
              irrep_y=irrep_g
              irrep_d=irrep_a
              irrep_tra=irrep_b
              irrep_trb=irrep_c
c              write(LuOut,*) "beta_2"
              call beta_2(d_d1(axisA),d_b0,d_tr1(axisB),d_tr2(axisB),
     1             d_tr1(axisC),d_tr2(axisC),d_y1,d_y2,
     2             k_d1_offset(axisA),k_b0_offset,
     3             k_tr1_offset(axisB),k_tr2_offset(axisB),
     4             k_tr1_offset(axisC),k_tr2_offset(axisC),
     5             k_y1_offset,k_y2_offset)
              call reconcilefile(d_b0,1)
              call get_block(d_b0,beta2,1,0)
              call tce_zero(d_b0,1)
c
              irrep_y=irrep_a
              irrep_tra=irrep_b
              irrep_trb=irrep_c
c              write(LuOut,*) "beta_3"
              call beta_3(d_f1,d_b0,d_t1,d_t2,d_tr1(axisB),d_tr2(axisB),
     1             d_tr1(axisC),d_tr2(axisC),d_v2,
     2             d_yr1(axisA),d_yr2(axisA),
     3             k_f1_offset,k_b0_offset,k_t1_offset,k_t2_offset,
     4             k_tr1_offset(axisB),k_tr2_offset(axisB),
     5             k_tr1_offset(axisC),k_tr2_offset(axisC),k_v2_offset,
     6             k_yr1_offset(axisA),k_yr2_offset(axisA))
              call reconcilefile(d_b0,1)
              call get_block(d_b0,beta3,1,0)
              call tce_zero(d_b0,1)
c
              irrep_y=irrep_a
              irrep_d=irrep_b
              irrep_tr=irrep_c
c              write(LuOut,*) "beta_4"
              call beta_4(d_d1(axisB),d_b0,d_t1,d_t2,
     1             d_tr1(axisC),d_tr2(axisC),d_yr1(axisA),d_yr2(axisA),
     2             k_d1_offset(axisB),k_b0_offset,
     3             k_t1_offset,k_t2_offset,
     4             k_tr1_offset(axisC),k_tr2_offset(axisC),
     5             k_yr1_offset(axisA),k_yr2_offset(axisA))
              call reconcilefile(d_b0,1)
              call get_block(d_b0,beta4,1,0)
              call tce_zero(d_b0,1)
c
c              write(LuOut,*) "beta1 = ",beta1
c              write(LuOut,*) "beta2 = ",beta2
c              write(LuOut,*) "beta3 = ",beta3
c              write(LuOut,*) "beta4 = ",beta4
c
c             This term contributes to only one component
c
              beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta1
c
c             These terms both contribute to three components
c
              beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta2
              beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta2
              beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta2
c
              beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta3
              beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta3
              beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta3
c
c             This term contributes to six components
c
              beta(axisA,axisB,axisC)=beta(axisA,axisB,axisC)-beta4
              beta(axisA,axisC,axisB)=beta(axisA,axisC,axisB)-beta4
              beta(axisB,axisA,axisC)=beta(axisB,axisA,axisC)-beta4
              beta(axisB,axisC,axisA)=beta(axisB,axisC,axisA)-beta4
              beta(axisC,axisA,axisB)=beta(axisC,axisA,axisB)-beta4
              beta(axisC,axisB,axisA)=beta(axisC,axisB,axisA)-beta4
c
              call deletefile(d_b0)
c
            endif ! respaxis(axisA,axisB,axisC)
          enddo ! axisB loop
        enddo ! axisB loop
      enddo ! axisA loop
c
      if (nodezero) then
        write(LuOut,*)
        write(LuOut,9120) "CCSD Quadratic Response"
        write(LuOut,9504)
        write(LuOut,9502)
        do axisA = 1, 3
          do axisB = 1, 3
            do axisC = 1, 3
              if (dabs(beta(axisA,axisB,axisC)).gt.1.0d-6)
     1        write(LuOut,9503) axisname(axisA),axisname(axisB),
     2                          axisname(axisC),beta(axisA,axisB,axisC),
     3                          ' / au'
            enddo
          enddo
        enddo
        write(LuOut,9502)
      endif
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
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
 9502 format(1x,'-----------------------------------------------')
 9503 format(1x,'beta(',a1,',',a1,',',a1,') = ',f18.9,a)
 9504 format(1x,'Static Hyperpolarizability')
 9505 format(1x,'Second Harmonic Generation (SHG)')
 9506 format(1x,'Optical Rectification (OR)')
 9507 format(1x,'Electro-Optic Pockels Effect (EOPE)')
 9508 format(1x,'beta(',f1.5,',',f1.5,',',f1.5,')')
      return
      end



      subroutine ccsd_lr_d1_build(omega,d_rho1,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_y1,d_y2,d_tr1,d_tr2,d_yr1,d_yr2,
     2           k_rho1_offset,
     3           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset,
     5           k_tr1_offset,k_tr2_offset,k_yr1_offset,k_yr2_offset)
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
      integer i,j,dummy,axis,axisA
      integer irrep_g
      parameter (irrep_g=0)
      integer d_rho1(3)
      integer d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_y1,d_y2
      integer k_rho1_offset(3)
      integer k_a0_offset
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_y1_offset,k_y2_offset
      integer d_tr1(9),d_tr2(9),d_yr1(6),d_yr2(6)
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer k_yr1_offset(3),k_yr2_offset(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      logical nodezero
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
#ifdef DEBUG_PRINT
      if (nodezero) write(6,*) 'top of ccsd_lr_d1_build'
      call util_flush(6)
#endif
c
c CCSD-QR evaluation step
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
      do axisA = 1, 3
        irrep_a=sym_abelian_axis(geom,axisA)
        if (respaxis(axisA)) then
c
#ifdef DEBUG_PRINT
            if (nodezero) then
              write(LuOut,*) "axisA = ",axisA
              call util_flush(6)
            endif
#endif
c
c           irrep_t=irrep_g
            irrep_y=irrep_g
            irrep_tr=irrep_a
            irrep_yr=irrep_a
c
c           irrep_tr symmetry
c
            call ccsd_lr_d1_hh(d_rho1(axisA),
     1           d_t1,d_tr1(axisA),d_tr2(axisA),
     2           d_y1,d_y2,d_yr1(axisA),
     3           k_rho1_offset(axisA),
     4           k_t1_offset,k_tr1_offset(axisA),k_tr2_offset(axisA),
     5           k_y1_offset,k_y2_offset,k_yr1_offset(axisA))
c
c           irrep_tr symmetry
c
            call ccsd_lr_d1_hp(d_rho1(axisA),
     1           d_t1,d_t2,d_tr1(axisA),d_tr2(axisA),
     2           d_y1,d_y2,d_yr1(axisA),d_yr2(axisA),
     3           k_rho1_offset(axisA),
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_y1_offset,k_y2_offset,
     7           k_yr1_offset(axisA),k_yr2_offset(axisA))
c
c           irrep_tr symmetry
c
            call ccsd_lr_d1_pp(d_rho1(axisA),
     1           d_t1,d_t2,d_tr1(axisA),d_tr2(axisA),
     2           d_y1,d_y2,d_yr1(axisA),d_yr2(axisA),
     3           k_rho1_offset(axisA),
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_y1_offset,k_y2_offset,
     7           k_yr1_offset(axisA),k_yr2_offset(axisA))
c
c           irrep_yr symmetry
c
            call ccsd_lr_d1_ph(d_rho1(axisA),
     1           d_t1,d_yr1(axisA),
     2           k_rho1_offset(axisA),
     3           k_t1_offset,k_yr1_offset(axisA))
c
        endif ! respaxis(axisA)
      enddo ! axisA loop
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
c
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
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
c $Id$
