      subroutine eomccsd_energy(d_rx1,d_rx2,size_x1,size_x2,
     1           k_x1_offset,k_x2_offset,
     1           d_f1,d_v2,d_t1,d_t2,k_f1_offset,k_v2_offset,
     2           k_t1_offset,k_t2_offset,k_omegax)
      implicit none
#include "global.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
c
      integer size_x1,size_x2
      integer k_x1_offset,k_x2_offset
      integer d_f1,d_v2
      integer d_t1,d_t2
      integer k_f1_offset,k_v2_offset
      integer k_t1_offset,k_t2_offset
c
c      integer iter,ivec
      integer ivec
      logical nodezero
      double precision cpu, wall
      logical needt1,needt2,needt3,needt4
      integer size_x3,size_x4
      integer k_x3_offset,k_x4_offset
      integer d_rx1,d_rx2,d_rx3,d_rx4
      integer k_residual,l_residual
      integer k_omegax, l_omegax
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
c
      character*255 filename
      character*255 modelname
      logical converged
c
      nodezero=(ga_nodeid().eq.0)
      needt1=.true.
      needt2=.true.
      needt3=.false.
      needt4=.false.
c
      call tce_eom_xguess(needt1,needt2,needt3,needt4,
     1     size_x1,size_x2,size_x3,size_x4,
     2     k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset)
      if (nxtrials.eq.0) goto 200
      modelname = "EOM-CCSD right-hand side"
      if (nodezero) write(LuOut,9220)
     1   modelname(1:inp_strlen(modelname))
      do iter=1,maxiter
         if (nodezero.and.util_print('eom',print_default))
     1      write(LuOut,9210) iter,nxtrials
         do ivec = 1,nxtrials
            if (.not.xp1_exist(ivec)) then
               call tce_filenameindexed(ivec,'xp1',filename)
               call createfile(filename,xp1(ivec),size_x1)
               xp1_exist(ivec) = .true.
               call dratoga(x1(ivec))
               call dratoga(x2(ivec))
               call eomccsd_x1_grad(d_f1,xp1(ivec),d_t1,d_t2,d_v2,
     1            x1(ivec),x2(ivec),k_f1_offset,k_x1_offset,
     2            k_t1_offset,k_t2_offset,k_v2_offset,
     3            k_x1_offset,k_x2_offset)
               call reconcilefile(xp1(ivec),size_x1)
               call gatodra(x2(ivec))
               call gatodra(x1(ivec))
               call gatodra(xp1(ivec))
            endif
            if (.not.xp2_exist(ivec)) then
               call tce_filenameindexed(ivec,'xp2',filename)
               call createfile(filename,xp2(ivec),size_x2)
               xp2_exist(ivec) = .true.
               call dratoga(x1(ivec))
               call dratoga(x2(ivec))
               call eomccsd_x2_grad(d_f1,xp2(ivec),d_t1,d_t2,d_v2,
     1              x1(ivec),x2(ivec),k_f1_offset,k_x2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_x1_offset,k_x2_offset)
               call reconcilefile(xp2(ivec),size_x2)
               call gatodra(x2(ivec))
               call gatodra(x1(ivec))
               call gatodra(xp2(ivec))
            endif
         enddo
         if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1       l_residual,k_residual))
     2       call errquit('eomccsd_energy: ma problem',10,ma_err)
         call tce_eom_xdiagon_grad(needt1,needt2,needt3,needt4,
     1        size_x1,size_x2,size_x3,size_x4,
     2        k_x1_offset,k_x2_offset,k_x3_offset,k_x4_offset,
     3        d_rx1,d_rx2,d_rx3,d_rx4,
     4        dbl_mb(k_omegax),dbl_mb(k_residual))
         cpu=cpu+util_cpusec()
         wall=wall+util_wallsec()
         converged = .true.
         do ivec = 1,nroots_reduced
            if (nodezero.and.(ivec.ne.nroots_reduced))
     1         write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2         dbl_mb(k_omegax+ivec-1),
     3         dbl_mb(k_omegax+ivec-1)*au2ev
            if (nodezero.and.(ivec.eq.nroots_reduced))
     1         write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2         dbl_mb(k_omegax+ivec-1),
     3         dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
            if (nodezero) call util_flush(LuOut)
            if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1         converged = .false.
         enddo
         cpu=-util_cpusec()
         wall=-util_wallsec()
         if (.not.ma_pop_stack(l_residual))
     1      call errquit("eomccsd_energy: ma problem",20,ma_err)
         if (converged) then
            call tce_eom_xtidy
            if (nodezero) then
               write(LuOut,9240)
               call util_flush(LuOut)
            endif
            goto 200
         endif
      enddo
      call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
  200 continue
 9210 format(/,1x,'Iteration ',i3,' using ',i4,' trial vectors')
 9220 format(/,1x,A,' iterations',/,1x,
     1'--------------------------------------------------------------'
     2,/,1x,
     3'     Residuum       Omega / hartree  Omega / eV    Cpu    Wall'
     4,/,1x,
     5'--------------------------------------------------------------')
 9230 format(1x,f17.13,f18.13,f11.5,2f8.1)
 9240 format(1x,
     1'--------------------------------------------------------------'
     2,/,1x,'Iterations converged')
      end
c $Id$
