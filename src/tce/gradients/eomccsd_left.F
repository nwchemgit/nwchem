      subroutine eomccsd_left(d_ry1,d_ry2,size_y1,size_y2,
     1                        k_y1_offset,k_y2_offset,
     2                        d_f1,d_v2,d_t1,d_t2,k_f1_offset,
     3                        k_v2_offset,k_t1_offset,k_t2_offset,
     4                        k_omegay)
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
      integer d_ry1,d_ry2
      integer size_y1,size_y2
      integer k_y1_offset,k_y2_offset
      integer d_f1,d_v2,d_t1,d_t2
      integer k_f1_offset,k_v2_offset
      integer k_t1_offset,k_t2_offset
      integer k_omegay
c
      logical needt1,needt2,needt3,needt4
      integer size_y3,size_y4,k_y3_offset,k_y4_offset
      integer d_ry3,d_ry4
c
      character*255 modelname, filename
c      integer iter,ivec
      integer ivec
      logical nodezero
      integer k_residual,l_residual
      double precision au2ev   ! Conversion factor from a.u. to eV
      parameter (au2ev=27.2113961d0)
      logical converged
      double precision cpu, wall
c
      nodezero=(ga_nodeid().eq.0)
c
      needt1=.true.
      needt2=.true.
      needt3=.false.
      needt4=.false.
c
      call tce_eom_yguess(needt1,needt2,needt3,needt4,
     1       size_y1,size_y2,size_y3,size_y4,
     2       k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset)
      if (nytrials.eq.0) goto 400
      modelname = "EOM-CCSD left-hand side"
      if (nodezero) write(LuOut,9220)
     1   modelname(1:inp_strlen(modelname))
      do iter=1,maxiter
         if (nodezero.and.util_print('eom',print_default))
     1      write(LuOut,9210) iter,nytrials
         do ivec = 1,nytrials
            if (.not.yp1_exist(ivec)) then
               call tce_filenameindexed(ivec,'yp1',filename)
               call createfile(filename,yp1(ivec),size_y1)
               yp1_exist(ivec) = .true.
               call dratoga(y1(ivec))
               call dratoga(y2(ivec))
               call eomccsd_y1(d_f1,yp1(ivec),d_t1,d_t2,d_v2,
     1              y1(ivec),y2(ivec),k_f1_offset,k_y1_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_y1_offset,k_y2_offset)
               call reconcilefile(yp1(ivec),size_y1)
               call gatodra(y2(ivec))
               call gatodra(y1(ivec))
               call gatodra(yp1(ivec))
            endif
            if (.not.yp2_exist(ivec)) then
               call tce_filenameindexed(ivec,'yp2',filename)
               call createfile(filename,yp2(ivec),size_y2)
               yp2_exist(ivec) = .true.
               call dratoga(y1(ivec))
               call dratoga(y2(ivec))
               call eomccsd_y2(d_f1,yp2(ivec),d_t1,d_t2,d_v2,
     1              y1(ivec),y2(ivec),k_f1_offset,k_y2_offset,
     2              k_t1_offset,k_t2_offset,k_v2_offset,
     3              k_y1_offset,k_y2_offset)
               call reconcilefile(yp2(ivec),size_y2)
               call gatodra(y2(ivec))
               call gatodra(y1(ivec))
               call gatodra(yp2(ivec))
            endif
         enddo
         if (.not.ma_push_get(mt_dbl,nytrials,'residual',
     1       l_residual,k_residual))
     2      call errquit('eomccsd_left: ma problem',10,ma_err)
         call tce_eom_ydiagon(needt1,needt2,needt3,needt4,
     1          size_y1,size_y2,size_y3,size_y4,
     2          k_y1_offset,k_y2_offset,k_y3_offset,k_y4_offset,
     3          d_ry1,d_ry2,d_ry3,d_ry4,
     4          dbl_mb(k_omegay),dbl_mb(k_residual))
         cpu=cpu+util_cpusec()
         wall=wall+util_wallsec()
         converged = .true.
         do ivec = 1,nroots_reduced
            if (nodezero.and.(ivec.ne.nroots_reduced))
     1         write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2         dbl_mb(k_omegay+ivec-1),
     3         dbl_mb(k_omegay+ivec-1)*au2ev
            if (nodezero.and.(ivec.eq.nroots_reduced))
     1         write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2         dbl_mb(k_omegay+ivec-1),
     3         dbl_mb(k_omegay+ivec-1)*au2ev,cpu,wall
            if (nodezero) call util_flush(LuOut)
            if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1         converged = .false.
         enddo
         cpu=-util_cpusec()
         wall=-util_wallsec()
         if (.not.ma_pop_stack(l_residual))
     1      call errquit("eomccsd_left: ma problem",20,ma_err)
         if (converged) then
            if (nodezero) write(LuOut,9240)
            call tce_eom_ytidy
             modelname = "EOM-CCSD"
             goto 400
         endif
      enddo
      call errquit('eomccsd_left: maxiter exceeded',iter,CALC_ERR)
  400 continue
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
