      subroutine ccsd_imaginary(d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     6           k_tr1_offset,k_tr2_offset,
     7           size_tr1,size_tr2)
      implicit none
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
      integer i,j,dummy,axis
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_lambda1,d_lambda2
      integer d_tr1(9),d_tr2(9)
      integer k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      integer zero,one
      external sym_abelian_axis
      double precision pi
      parameter(pi = 3.14159265358979323846264338327950288419D0)
      parameter(zero = 0)
      parameter(one  = 0)
      double precision omega
      double precision threeoverpi
      parameter (threeoverpi=3.0d0/pi)
      logical nodezero,guess
      character*4 irrepname
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
c
      nodezero=(ga_nodeid().eq.0)
c
      if(nodezero) then
        write(6,*)
        write(6,*) 'Casimir-Polder integration points are:'
        do i = 1, inumfreq
          write(6,999) i,ifreq(i)
        enddo
      endif
c
      do omegacount=1,inumfreq
        omega = ifreq(omegacount)
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
            if (nodezero) write(LuOut,9431) omega
c
            if (cc_ir_alg.eq.1) then
c
c             REAL COMPONENT
c
              dynaxis = 6
              if (guess_ir_real.and.(omegacount.eq.1)) then
                call ccsd_ir_guess(d_d1,d_t1,d_t2,d_tr1,d_tr2,
     1             k_d1_offset,k_t1_offset,k_t2_offset,
     2             k_tr1_offset,k_tr2_offset,
     4             size_tr1,size_tr2,zero,axis,omega)
              endif ! guess_ir_real
              call ccsd_ir_real_sq_it(axis,dynaxis,omega,
     1             d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2             k_f1_offset,k_v2_offset,k_d1_offset,
     3             k_t1_offset,k_t2_offset,
     4             k_tr1_offset,k_tr2_offset,
     5             size_tr1,size_tr2)
c
c             IMAGINARY COMPONENT
c
              if (guess_ir_imag.and.(omegacount.eq.1)) then
                call ccsd_ir_get_imag_from_real(axis,omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     4               size_tr1,size_tr2)
              endif
c
              dynaxis = 0 ! positive omega_I
              call ccsd_ir_imag_it(axis,dynaxis,omega,
     1             d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2             k_f1_offset,k_v2_offset,k_d1_offset,
     3             k_t1_offset,k_t2_offset,
     4             k_tr1_offset,k_tr2_offset,
     5             size_tr1,size_tr2)
              call daxfile(1,-1.0d0,d_tr1(axis+0),d_tr1(axis+3),
     1                      size_tr1(axis))
              call daxfile(1,-1.0d0,d_tr2(axis+0),d_tr2(axis+3),
     2                      size_tr2(axis))
#ifdef DEBUG_ONLY
              dynaxis = 3 ! is for negative omega_I
              call ccsd_ir_imag_it(axis,dynaxis,(-1.0d0*omega),
     1             d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2             k_f1_offset,k_v2_offset,k_d1_offset,
     3             k_t1_offset,k_t2_offset,
     4             k_tr1_offset,k_tr2_offset,
     5             size_tr1,size_tr2)
#endif
c
            elseif (cc_ir_alg.eq.2) then
              if (omega.gt.(1.0d-3)) then
c
c               IMAGINARY COMPONENT
c
                dynaxis = 0 ! positive omega_I
                call ccsd_ir_imag_sq_it(axis,dynaxis,1.0d0/omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,
     4               k_tr1_offset,k_tr2_offset,
     5               size_tr1,size_tr2)
                call daxfile(1,-1.0d0,d_tr1(axis+0),d_tr1(axis+3),
     1                       size_tr1(axis))
                call daxfile(1,-1.0d0,d_tr2(axis+0),d_tr2(axis+3),
     2                       size_tr2(axis))
#ifdef DEBUG_ONLY
                dynaxis = 3 ! negative omega_I
                call ccsd_ir_imag_sq_it(axis,dynaxis,(-1.0d0/omega),
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,
     4               k_tr1_offset,k_tr2_offset,
     5               size_tr1,size_tr2)
#endif
c
c               REAL COMPONENT
c
                dynaxis = 6 ! offset to the third set of response for real component
                call ccsd_ir_real_it(axis,dynaxis,1.0d0/omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,
     4               k_tr1_offset,k_tr2_offset,
     5               size_tr1,size_tr2)
              elseif (omega.eq.(0.0d0)) then
c
c               IMAGINARY COMPONENT
c
                call tce_zero(d_tr1(axis+0),size_tr1(axis))
                call tce_zero(d_tr2(axis+0),size_tr2(axis))
                call tce_zero(d_tr1(axis+3),size_tr1(axis))
                call tce_zero(d_tr2(axis+3),size_tr2(axis))
c
c               REAL COMPONENT
c
                dynaxis = 6 ! offset to the third set of response for real component
                call ccsd_lr_iter(axis,dynaxis,omega,
     1               d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2               k_f1_offset,k_v2_offset,k_d1_offset,
     3               k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     4               size_tr1,size_tr2)
              else
                if (nodezero) then
                  write(9121) 'Frequency is too close to zero.'
                endif
              endif
c
            elseif (cc_ir_alg.eq.3) then
c
              guess = .true. !(omegacount.le.2)
c
              call ccsd_ir_coupled_it(guess,axis,omega,
     1             d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2             k_f1_offset,k_v2_offset,k_d1_offset,
     3             k_t1_offset,k_t2_offset,
     4             k_tr1_offset,k_tr2_offset,
     5             size_tr1,size_tr2)
c
            endif ! cc_ir_alg
c
          endif ! respaxis(axis)
        enddo ! axis loop
c
c CCSD-IR evaluation step
c
        call  ccsd_ir_eval(omegacount,d_a0,d_f1,d_v2,d_d1,
     1                     d_t1,d_t2,d_lambda1,d_lambda2,
     2                     d_tr1,d_tr2,k_a0_offset,
     3                     k_f1_offset,k_v2_offset,k_d1_offset,
     4                     k_t1_offset,k_t2_offset,
     5                     k_l1_offset,k_l2_offset,
     6                     k_tr1_offset,k_tr2_offset)
c
      enddo ! omegacount loop
c
c     Casimir-Polder integration for determination of C6 coefficient
c
      if (ifreqauto) then
        do i = 1, 8
          integral(i) = 0.0d0
          do omegacount = 1,inumfreq-1
            num0 = ifreqval(omegacount,i)
            num1 = ifreqval(omegacount+1,i)
            den0 = ifreq(omegacount)
            den1 = ifreq(omegacount+1)
            num0sq = num0*num0
            num1sq = num1*num1
            numtot = num1sq+num0sq
            dentot = den1-den0
            integral(i) = integral(i) + 0.5d0*numtot*dentot
          enddo
          integral(i) = integral(i)*threeoverpi
        enddo
c
        if (nodezero) write(LuOut,9435) "CCSD Imaginary Response",
     1            integral(1),integral(2),integral(3),integral(4),
     2            integral(5),integral(6),integral(7),integral(8)
c
      endif
c
  999 format(3x,'ifreq(',i2,') = ',f14.8)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9420 format(1x,i4,f25.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
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
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end


      subroutine ccsd_ir_real_sq_it(axis,dynaxis,omega,
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
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer axisA
      integer axisB
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9),d_rr1(3),d_rr2(3),d_ir1(3),d_ir2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision rr1,rr2,residual
      double precision omega
      double precision cpu, wall
      double precision ddotfile
      external ddotfile
      logical nodezero
      character*4 irrepname
      character*255 filename
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1))
     &    write(LuOut,9400) "CCSD-IR (real component) - squared"
c
        call tce_filename(ir1filename(axis),filename)
        call createfile(filename,d_ir1(axis),size_tr1(axis))
        call tce_zero(d_ir1(axis),size_tr1(axis))
c
        call tce_filename(ir2filename(axis),filename)
        call createfile(filename,d_ir2(axis),size_tr2(axis))
        call tce_zero(d_ir2(axis),size_tr2(axis))
c
        call ccsd_o1(d_ir1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x1_old(d_f1,d_ir1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),k_f1_offset,
     2       k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_ir1(axis),size_tr1(axis))
c
        call ccsd_o2(d_ir2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x2_old(d_f1,d_ir2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_ir2(axis),size_tr2(axis))
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_ir1(axis),d_ir2(axis),k_f1_offset,
     2       k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
c
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_ir1(axis),d_ir2(axis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
c
        call deletefile(d_ir1(axis))
        call deletefile(d_ir2(axis))
c
        call daxpyfile(1,omega*omega,d_tr1(axis+dynaxis),
     1       d_rr1(axis),size_tr1(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call daxpyfile(1,omega*omega,d_tr2(axis+dynaxis),
     1       d_rr2(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
        residual = max(rr1,rr2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) write(LuOut,9410)
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
        call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1       d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2       size_tr1(axis),
     3       d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4       size_tr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,omega,2)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_ir_real_sq_it: maxiter exceeded',iter,CALC_ERR)
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


      subroutine ccsd_ir_imag_sq_it(axis,dynaxis,omega,
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
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis
      integer omegacount,omegasign
      integer dynfreq,dynaxis
      integer axisA,axisB
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9),d_rr1(3),d_rr2(3),d_ir1(3),d_ir2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision rr1,rr2,residual
      double precision omega
      double precision cpu, wall
      double precision ddotfile
      external ddotfile
      logical nodezero
      character*4 irrepname
      character*255 filename
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1))
     &    write(LuOut,9400) "CCSD-IR (imaginary component) - squared"
c
        call tce_filename(ir1filename(axis),filename)
        call createfile(filename,d_ir1(axis),size_tr1(axis))
        call tce_zero(d_ir1(axis),size_tr1(axis))
c
        call tce_filename(ir2filename(axis),filename)
        call createfile(filename,d_ir2(axis),size_tr2(axis))
        call tce_zero(d_ir2(axis),size_tr2(axis))
c
        call eomccsd_x1_old(d_f1,d_ir1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),k_f1_offset,
     2       k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_ir1(axis),size_tr1(axis))
c
        call eomccsd_x2_old(d_f1,d_ir2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_ir2(axis),size_tr2(axis))
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_ir1(axis),d_ir2(axis),k_f1_offset,
     2       k_tr1_offset(axis),k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
c
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_ir1(axis),d_ir2(axis),k_f1_offset,
     2       k_tr2_offset(axis),k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
c
        call tce_zero(d_ir1(axis),size_tr1(axis))
        call ccsd_o1(d_ir1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call daxpyfile(1,omega,d_ir1(axis),d_rr1(axis),size_tr1(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
c
        call tce_zero(d_ir2(axis),size_tr2(axis))
        call ccsd_o2(d_ir2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call daxpyfile(1,omega,d_ir2(axis),d_rr2(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
c
        call deletefile(d_ir1(axis))
        call deletefile(d_ir2(axis))
c
        call daxpyfile(1,omega*omega,d_tr1(axis+dynaxis),
     1       d_rr1(axis),size_tr1(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call daxpyfile(1,omega*omega,d_tr2(axis+dynaxis),
     1       d_rr2(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
        residual = max(rr1,rr2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) write(LuOut,9410)
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
        call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1       d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2       size_tr1(axis),
     3       d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4       size_tr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,omega,2)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_ir_imag_sq_it: maxiter exceeded',iter,CALC_ERR)
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



      subroutine ccsd_ir_real_it(axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,
     4           k_tr1_offset,k_tr2_offset,
     5           size_tr1,size_tr2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9),d_rr1(3),d_rr2(3),d_ir1(3),d_ir2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision rr1,rr2,residual
      double precision omega
      double precision cpu, wall
      double precision ddotfile
      external ddotfile
      logical nodezero
      character*4 irrepname
      character*255 filename
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1))
     &    write(LuOut,9400) "CCSD-IR (real component) - linear"
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call daxfile(1,(-1.0d0)*omega,d_tr1(axis),
     1       d_rr1(axis),size_tr1(axis))
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call daxfile(1,(-1.0d0)*omega,d_tr2(axis),
     1       d_rr2(axis),size_tr2(axis))
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
        residual = max(rr1,rr2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) write(LuOut,9410)
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
        call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1       d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2       size_tr1(axis),
     3       d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4       size_tr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,0.0d0,1)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_ir_real_it: maxiter exceeded',iter,CALC_ERR)
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


      subroutine ccsd_ir_imag_it(axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,
     4           k_tr1_offset,k_tr2_offset,
     5           size_tr1,size_tr2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2,d_tr1(9),d_tr2(9)
      integer d_rr1(3),d_rr2(3),d_ir1(3),d_ir2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision rr1,rr2,residual
      double precision omega
      double precision cpu, wall
      double precision ddotfile
      external ddotfile
      logical nodezero
      character*4 irrepname
      character*255 filename
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
        if (nodezero.and.(iter.eq.1))
     &    write(LuOut,9400) "CCSD-IR (imaginary component) - linear"
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call daxfile(1,(-1.0d0)*omega,d_tr1(axis+6),
     1       d_rr1(axis),size_tr1(axis))
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call daxfile(1,(-1.0d0)*omega,d_tr2(axis+6),
     1       d_rr2(axis),size_tr2(axis))
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
        residual = max(rr1,rr2)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) write(LuOut,9410)
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
        call tce_diis3(.false.,iter,.true.,.true.,.false.,.false.,
     1       d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2       size_tr1(axis),
     3       d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4       size_tr2(axis),
     5       dummy,dummy,dummy,dummy,
     6       dummy,dummy,dummy,dummy,0.0d0,1)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_ir_imag_it: maxiter exceeded',iter,CALC_ERR)
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



      subroutine ccsd_ir_coupled_it(guess,axis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,
     4           k_tr1_offset,k_tr2_offset,
     5           size_tr1,size_tr2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis,dynaxis
      integer omegacount
      integer omegasign
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9)
      integer d_rr1(3),d_rr2(3),d_ir1(3),d_ir2(3)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      integer zero,one
      parameter(zero = 0)
      parameter(one  = 1)
      external sym_abelian_axis
      double precision rr1,rr2,rr3,rr4,rres,ires,residual
      double precision omega
      double precision cpu, wall
      double precision ddotfile
      external ddotfile
      logical nodezero,guess
      character*4 irrepname
      character*255 filename
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
c
      nodezero=(ga_nodeid().eq.0)
c
      call tce_diis_init()
      do iter=1,maxiter
        cpu=-util_cpusec()
        wall=-util_wallsec()
c
c       REAL COMPONENT
c
        if (guess.and.(iter.eq.1).and.guess_ir_real) then
          call ccsd_lr_guess(d_f1,d_v2,d_d1,d_t1,d_t2,d_tr1,d_tr2,
     1         k_f1_offset,k_v2_offset,k_d1_offset,
     2         k_t1_offset,k_t2_offset,k_tr1_offset,k_tr2_offset,
     3         size_tr1,size_tr2,axis,6)
        endif
c
        if (guess.and.(iter.eq.1).and.guess_ir_imag) then
          call ccsd_ir_guess(d_d1,d_t1,d_t2,d_tr1,d_tr2,
     1         k_d1_offset,k_t1_offset,k_t2_offset,
     2         k_tr1_offset,k_tr2_offset,
     4         size_tr1,size_tr2,one,axis,omega)
        endif
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
        call daxfile(1,(-1.0d0)*omega,d_tr1(axis+0),
     1       d_rr1(axis),size_tr1(axis))
        call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x1_old(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+6),d_tr2(axis+6),
     2       k_f1_offset,k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_rr1(axis),size_tr1(axis))
        call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
c
        call daxfile(1,(-1.0d0)*omega,d_tr2(axis+0),
     1       d_rr2(axis),size_tr2(axis))
        call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x2_old(d_f1,d_rr2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+6),d_tr2(axis+6),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_rr2(axis),size_tr2(axis))
        call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
c
c       IMAGINARY COMPONENT
c
        call tce_filename(ir1filename(axis),filename)
        call createfile(filename,d_ir1(axis),size_tr1(axis))
        call tce_zero(d_ir1(axis),size_tr1(axis))
c
        call tce_filename(ir2filename(axis),filename)
        call createfile(filename,d_ir2(axis),size_tr2(axis))
        call tce_zero(d_ir2(axis),size_tr2(axis))
c
        call daxfile(1,(-1.0d0)*omega,d_tr1(axis+6),
     1       d_ir1(axis),size_tr1(axis))
        call eomccsd_x1_old(d_f1,d_ir1(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+0),d_tr2(axis+0),
     2       k_f1_offset,k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
        call reconcilefile(d_ir1(axis),size_tr1(axis))
        call tce_residual_tr1(d_ir1(axis),k_tr1_offset(axis),rr3)
c
        call daxfile(1,(-1.0d0)*omega,d_tr2(axis+6),
     1       d_ir2(axis),size_tr2(axis))
        call eomccsd_x2_old(d_f1,d_ir2(axis),d_t1,d_t2,d_v2,
     1       d_tr1(axis+0),d_tr2(axis+0),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
        call reconcilefile(d_ir2(axis),size_tr2(axis))
        call tce_residual_tr2(d_ir2(axis),k_tr2_offset(axis),rr4)
c
        if (nodezero.and.(iter.eq.1)) then
            write(LuOut,9050) "CCSD-IR (coupled iteration)"
        endif
c
        rres = max(rr1,rr2)
        ires = max(rr3,rr4)
        residual = max(rres,ires)
        cpu=cpu+util_cpusec()
        wall=wall+util_wallsec()
        if (nodezero) write(LuOut,9100) iter,rres,ires,cpu,wall
        if (residual .lt. thresh) then
          if (nodezero) write(LuOut,9060)
          call deletefile(d_rr2(axis))
          call deletefile(d_rr1(axis))
          call deletefile(d_ir2(axis))
          call deletefile(d_ir1(axis))
          call tce_diis_tidy()
          return
        endif
        call tce_diis3c(iter,.true.,.true.,.true.,.true.,
     1      d_rr1(axis),d_tr1(axis+6),k_tr1_offset(axis),size_tr1(axis),
     2      d_rr2(axis),d_tr2(axis+6),k_tr2_offset(axis),size_tr2(axis),
     3      d_ir1(axis),d_tr1(axis+0),k_tr1_offset(axis),size_tr1(axis),
     4      d_ir2(axis),d_tr2(axis+0),k_tr2_offset(axis),size_tr2(axis),
     5      0.0d0)
        call deletefile(d_rr2(axis))
        call deletefile(d_rr1(axis))
        call deletefile(d_ir2(axis))
        call deletefile(d_ir1(axis))
        if (nodezero) call util_flush(LuOut)
      enddo ! iter loop
      call errquit('ccsd_ir_coupled_it: maxiter exceeded',iter,CALC_ERR)
c
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9050 format(/,1x,A,' iterations',/,
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iter    Residuum (real)  Residuum (imag)     Cpu    Wall',/
     3  1x,'--------------------------------------------------------')
 9060 format(
     1  1x,'--------------------------------------------------------',/
     2  1x,'Iterations converged')
 9100 format(1x,i4,2f18.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end


      subroutine ccsd_ir_get_imag_from_real(axis,omega,
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
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9)
      integer k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision neg
      parameter (neg=-1.0d0)
      logical nodezero
      character*4 irrepname
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
      cpu=-util_cpusec()
      wall=-util_wallsec()
c
      if (omega.eq.(0.0d0)) then
c
c       x_I = 0 if w_I = 0
c
        call tce_zero(d_tr1(axis+0),size_tr1(axis))
        call tce_zero(d_tr2(axis+0),size_tr2(axis))
        call tce_zero(d_tr1(axis+3),size_tr1(axis))
        call tce_zero(d_tr2(axis+3),size_tr2(axis))
        return
      else
c
c       x_I(+) = ( A * x_R - b ) / omega_I
c       x_I(-) = -x_I(+)
c
        call tce_zero(d_tr1(axis+0),size_tr1(axis))
        call tce_zero(d_tr2(axis+0),size_tr2(axis))
c
        call ccsd_o1(d_tr1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x1_old(d_f1,d_tr1(axis+0),d_t1,d_t2,d_v2,
     1       d_tr1(axis+6),d_tr2(axis+6),k_f1_offset,
     2       k_tr1_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis))
c
        call ccsd_o2(d_tr2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call eomccsd_x2_old(d_f1,d_tr2(axis+0),d_t1,d_t2,d_v2,
     1       d_tr1(axis+6),d_tr2(axis+6),
     2       k_f1_offset,k_tr2_offset(axis),
     3       k_t1_offset,k_t2_offset,k_v2_offset,
     4       k_tr1_offset(axis),k_tr2_offset(axis),
     5       size_tr1(axis),size_tr2(axis))
c
        call reconcilefile(d_tr1(axis+0),size_tr1(axis))
        call reconcilefile(d_tr2(axis+0),size_tr2(axis))
c
        call dscalfile((1.0d0/omega),d_tr1(axis+0),size_tr1(axis))
        call dscalfile((1.0d0/omega),d_tr2(axis+0),size_tr2(axis))
c
c        call daxfile(1,neg,d_tr1(axis+0),d_tr1(axis+3),size_tr1(axis))
c        call daxfile(1,neg,d_tr2(axis+0),d_tr2(axis+3),size_tr2(axis))
c        call copyfile(d_tr1(axis+0),d_tr1(axis+3),size_tr1(axis))
c        call copyfile(d_tr2(axis+0),d_tr2(axis+3),size_tr2(axis))
        return
      endif
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
      end


      subroutine ccsd_ir_eval(omegacount,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_lambda1,d_lambda2,d_tr1,d_tr2,
     2           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset,
     6           k_tr1_offset,k_tr2_offset)
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
      integer axisA
      integer axisB
      integer omegacount
      integer omegasign
      integer dynfreq
      integer dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_a0,d_f1,d_v2,d_d1(3)
      integer d_t1,d_t2,d_lambda1,d_lambda2
      integer k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset
      integer d_tr1(9),d_tr2(9)
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
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
      character*4 ir1filename(3) ! File name stub
      data ir1filename/'ir1x','ir1y','ir1z'/
      character*4 ir2filename(3) ! File name stub
      data ir2filename/'ir2x','ir2y','ir1z'/
      character*5 rr1filename(3) ! File name stub
      data rr1filename/'rr1x ','rr1y ','rr1z '/
      character*5 rr2filename(3) ! File name stub
      data rr2filename/'rr2x ','rr2y ','rr2z '/
c
      nodezero=(ga_nodeid().eq.0)
c
c CCSD-IR evaluation step
c
      cpu=-util_cpusec()
      wall=-util_wallsec()

      omega = ifreq(omegacount)

      do axisA = 1, 3
        do axisB = 1, axisA
          alpha(axisA,axisB)=0.0d0
          if (respaxis(axisA).and.respaxis(axisB)) then
            irrep_a=sym_abelian_axis(geom,axisA)
            irrep_b=sym_abelian_axis(geom,axisB)
            irrep_y=irrep_g
c
            call tce_filename('a0',filename)
            call createfile(filename,d_a0,1) ! size_a0 = 1
c
c            if (nodezero) write(LuOut,*) "axisA = ",axisA
c            if (nodezero) write(LuOut,*) "axisB = ",axisB
c
            alpha1=0.0d0
            alpha2=0.0d0
c
            irrep_d=irrep_a
            irrep_tr=irrep_b
c
c            write(LuOut,*) "alpha_1 AB"
            call alpha_1(d_d1(axisA),d_a0,d_t1,d_t2,
     1           d_tr1(axisB+6),d_tr2(axisB+6),
     2           d_lambda1,d_lambda2,
     3           k_d1_offset(axisA),k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisB),k_tr2_offset(axisB),
     6           k_l1_offset,k_l2_offset)
c
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha1,1,0)
            call tce_zero(d_a0,1)
c            write(LuOut,*) "alpha_1 AB = ",alpha1
c
            if (axisA.eq.axisB) then
              alpha2=alpha1
            else
              irrep_d=irrep_b
              irrep_tr=irrep_a
c              write(LuOut,*) "alpha_1 BA"
              call alpha_1(d_d1(axisB),d_a0,d_t1,d_t2,
     1             d_tr1(axisA+6),d_tr2(axisA+6),
     2             d_lambda1,d_lambda2,
     3             k_d1_offset(axisB),k_a0_offset,
     4             k_t1_offset,k_t2_offset,
     5             k_tr1_offset(axisA),k_tr2_offset(axisA),
     6             k_l1_offset,k_l2_offset)
c
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,1)
            endif ! axisA.eq.axisB
c            write(LuOut,*) "alpha_1 BA = ",alpha2
c
            alpha(axisA,axisB)=alpha(axisA,axisB)+alpha1+alpha2
c
            irrep_tra=irrep_a
            irrep_trb=irrep_b
c
            if (omega.ne.(0.0d0)) then
              alpha1=0.0d0
              alpha2=0.0d0
c              write(LuOut,*) "alpha_2 A(+i) B(-i)"
              call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1             d_tr1(axisA),d_tr2(axisA),
     2             d_tr1(axisB+3),d_tr2(axisB+3),
     3             d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4             k_t1_offset,k_t2_offset,
     5             k_tr1_offset(axisA),k_tr2_offset(axisA),
     6             k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7             k_l1_offset,k_l2_offset)
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha1,1,0)
              call tce_zero(d_a0,1)
c              write(LuOut,*) "alpha_2 A(+i) B(-i) = ",alpha1
c
c              write(LuOut,*) "alpha_2 A(-i) B(+i)"
              call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1             d_tr1(axisA+3),d_tr2(axisA+3),
     2             d_tr1(axisB),d_tr2(axisB),
     3             d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4             k_t1_offset,k_t2_offset,
     5             k_tr1_offset(axisA),k_tr2_offset(axisA),
     6             k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7             k_l1_offset,k_l2_offset)
              call reconcilefile(d_a0,1)
              call get_block(d_a0,alpha2,1,0)
              call tce_zero(d_a0,1)
c              write(LuOut,*) "alpha_2 A(-i) B(+i) = ",alpha2
            elseif (omega.eq.(0.0d0)) then
              alpha1 = 0.0d0
              alpha2 = 0.0d0
            endif
c
            alpha3=0.0d0
c            write(LuOut,*) "alpha_2 A(r) B(r)"
            call alpha_2(d_f1,d_a0,d_t1,d_t2,
     1           d_tr1(axisA+6),d_tr2(axisA+6),
     2           d_tr1(axisB+6),d_tr2(axisB+6),
     3           d_v2,d_lambda1,d_lambda2,k_f1_offset,k_a0_offset,
     4           k_t1_offset,k_t2_offset,
     5           k_tr1_offset(axisA),k_tr2_offset(axisA),
     6           k_tr1_offset(axisB),k_tr2_offset(axisB),k_v2_offset,
     7           k_l1_offset,k_l2_offset)
            call reconcilefile(d_a0,1)
            call get_block(d_a0,alpha3,1,0)
c            write(LuOut,*) "alpha_2 A(r) B(r) = ",alpha3
            call deletefile(d_a0)
c
            alpha(axisA,axisB)=alpha(axisA,axisB)-alpha1-alpha2+alpha3
            alpha(axisA,axisB)=-1.0d0*alpha(axisA,axisB)
c
c            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
c            if (nodezero) write(LuOut,9020) cpu, wall
c
          endif ! respaxis(axis)
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
      if (nodezero) write(LuOut,9434) "CCSD Imaginary Response",
     1  ifreq(omegacount),auXnm/ifreq(omegacount),
     2  alpha(1,1),alpha(2,1),alpha(3,1),
     3  au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4  alpha(2,1),alpha(2,2),alpha(3,2),
     5  au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6  alpha(3,1),alpha(3,2),alpha(3,3),
     7  au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8  alphaevr(1),alphaevr(2),alphaevr(3),
     9  au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1  alphaiso,au2ang3*alphaiso,alphaani,au2ang3*alphaani
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
c
      ifreqval(omegacount,1) = alpha(1,1)
      ifreqval(omegacount,2) = alpha(2,2)
      ifreqval(omegacount,3) = alpha(3,3)
      ifreqval(omegacount,4) = alpha(2,1)
      ifreqval(omegacount,5) = alpha(3,1)
      ifreqval(omegacount,6) = alpha(3,2)
      ifreqval(omegacount,7) = alphaiso
      ifreqval(omegacount,8) = alphaani
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
      return
      end




      subroutine ccsd_ir_guess(d_d1,d_t1,d_t2,d_tr1,d_tr2,
     1           k_d1_offset,k_t1_offset,k_t2_offset,
     2           k_tr1_offset,k_tr2_offset,
     3           size_tr1,size_tr2,component,axis,omega)
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
      integer axis,component
      integer d_d1(3),d_t1,d_t2
      integer d_tr1(9),d_tr2(9),d_rr1(3),d_rr2(3)
      integer k_d1_offset(3),k_t1_offset,k_t2_offset
      integer k_tr1_offset(3),k_tr2_offset(3)
      integer size_tr1(3),size_tr2(3)
      double precision omega
      double precision cpu, wall
      logical nodezero
      character*4 irrepname
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
      if (component.eq.0) then
c
        if (nodezero) write(6,9121) 'Initial guess x = b/Adiag'
c
c       0. Create dummies
c
        call tce_filename(rr1filename(axis),filename)
        call createfile(filename,d_rr1(axis),size_tr1(axis))
        call tce_zero(d_rr1(axis),size_tr1(axis))
c
        call tce_filename(rr2filename(axis),filename)
        call createfile(filename,d_rr2(axis),size_tr2(axis))
        call tce_zero(d_rr2(axis),size_tr2(axis))
c
c       1. Form b
c
        call ccsd_o1(d_rr1(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr1_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
        call ccsd_o2(d_rr2(axis),d_d1(axis),d_t1,d_t2,
     1       k_tr2_offset(axis),k_d1_offset(axis),
     2       k_t1_offset,k_t2_offset)
c
c       2. Hit with preconditioner
c
        call tce_jacobi_ir1(d_rr1(axis),d_tr1(axis+6),
     1                       k_tr1_offset(axis),0.0d0,0.0d0,1.0d0)
        call tce_jacobi_ir2(d_rr2(axis),d_tr2(axis+6),
     1                       k_tr2_offset(axis),0.0d0,0.0d0,1.0d0)
c
c       3. Delete dummies
c
        call deletefile(d_rr1(axis))
        call deletefile(d_rr2(axis))
c
      elseif (component.eq.1) then
c
c     Initial guess x_I = omega*x_R/(Adiag-omega)
        if (nodezero) write(6,9121)
     1     'Initial guess x_I = omega*x_R/Adiag'
c
        if (omega.eq.(0.0d0)) then
          call tce_zero(d_tr1(axis+0),size_tr1(axis))
          call tce_zero(d_tr2(axis+0),size_tr2(axis))
          call tce_zero(d_tr1(axis+3),size_tr1(axis))
          call tce_zero(d_tr2(axis+3),size_tr2(axis))
          return
        endif
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
c     1. Form w_I = omega*x_R
c
        call daxfile(1,omega,d_tr1(axis+6),
     1       d_rr1(axis),size_tr1(axis))
        call daxfile(1,omega,d_tr2(axis+6),
     1       d_rr2(axis),size_tr2(axis))
c
c     3. Hit with preconditioner
c
        if (omega.gt.(0.0d0)) then
          call tce_jacobi_ir1(d_rr1(axis),d_tr1(axis+0),
     1         k_tr1_offset(axis),0.0d0,0.0d0,1.0d0)
          call tce_jacobi_ir2(d_rr2(axis),d_tr2(axis+0),
     1         k_tr2_offset(axis),0.0d0,0.0d0,1.0d0)
        elseif (omega.lt.(0.0d0)) then
          call tce_jacobi_ir1(d_rr1(axis),d_tr1(axis+3),
     1         k_tr1_offset(axis),0.0d0,0.0d0,1.0d0)
          call tce_jacobi_ir2(d_rr2(axis),d_tr2(axis+3),
     1         k_tr2_offset(axis),0.0d0,0.0d0,1.0d0)
        endif
c
c     4. Delete dummies
c
        call deletefile(d_rr1(axis))
        call deletefile(d_rr2(axis))
c
      endif ! component
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
