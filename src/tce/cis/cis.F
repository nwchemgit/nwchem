      subroutine cis(LuOut,size_x1,k_x1_offset,d_f1,k_f1_offset,d_v2,
     1 k_v2_offset,k_omegax,l_omegax,d_rx1)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
c     
c
      integer LuOut
      integer size_x1
      integer k_x1_offset 
      integer d_f1
      integer k_f1_offset
      integer d_v2
      integer k_v2_offset
      integer k_omegax,l_omegax
      integer k_residual,l_residual
      integer d_rx1
      logical converged
      logical nodezero
      integer ivec
      integer maxiter_aa
      integer iopen_f
      double precision cpu, wall
      double precision au2ev
      parameter (au2ev=27.2113961d0)
      character*255 filename
c    
      maxiter_aa=200
      nodezero=(ga_nodeid().eq.0)
      iopen_f=1
c
c      cpu=-util_cpusec()
c      wall=-util_wallsec()
c
      call cis_xguess(size_x1,k_x1_offset,iopen_f)
c
c      cpu=cpu+util_cpusec()
c      wall=wall+util_wallsec()
c      if(nodezero) then
c       write(6,*)'CIS guess cpu vs. wall time: ',cpu,wall
c       call util_flush(6)
c      end if
c
      if (nxtrials.eq.0) goto 200
      do iter=1,maxiter_aa
         if (nodezero.and.util_print('eom',print_default))
     1        write(LuOut,9210) iter,nxtrials
      cpu=-util_cpusec()
      wall=-util_wallsec()
         do ivec = 1,nxtrials
            if (.not.xp1_exist(ivec)) then
               call tce_filenameindexed(ivec,'xp1',filename)
               call createfile(filename,xp1(ivec),size_x1)
               xp1_exist(ivec) = .true.
               call dratoga(x1(ivec))
               call cis_x1(d_f1,xp1(ivec),d_v2,x1(ivec),
     1              k_f1_offset,k_x1_offset,k_v2_offset,k_x1_offset)
               call reconcilefile(xp1(ivec),size_x1)
               call gatodra(x1(ivec))
               call gatodra(xp1(ivec))
            endif
         enddo
c      cpu=cpu+util_cpusec()
c      wall=wall+util_wallsec()
c      if(nodezero) then
c       write(6,*)'CIS HC1 time: ',cpu,wall
c       call util_flush(6)
c      end if

         if (.not.ma_push_get(mt_dbl,nxtrials,'residual',
     1        l_residual,k_residual))
     2        call errquit('cis.F: MA problem',101,MA_ERR)
c      cpu=-util_cpusec()
c      wall=-util_wallsec()
         call cis_xdiagon(size_x1, k_x1_offset,d_rx1,
     1        dbl_mb(k_omegax),dbl_mb(k_residual))
c      cpu=cpu+util_cpusec()
c      wall=wall+util_wallsec()
c      if(nodezero) then
c       write(6,*)'CIS xdiagon time: ',cpu,wall
c       call util_flush(6)
c      end if
c
         converged = .true.
         do ivec = 1,nroots_reduced
            if (nodezero.and.(ivec.ne.nroots_reduced))
     1           write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2           dbl_mb(k_omegax+ivec-1),
     3           dbl_mb(k_omegax+ivec-1)*au2ev
            if (nodezero.and.(ivec.eq.nroots_reduced)) then
              cpu=cpu+util_cpusec()
              wall=wall+util_wallsec()
                 write(LuOut,9230) dbl_mb(k_residual+ivec-1),
     2           dbl_mb(k_omegax+ivec-1),
     3           dbl_mb(k_omegax+ivec-1)*au2ev,cpu,wall
            end if
            if (nodezero) call util_flush(LuOut)
            if (dbl_mb(k_residual+ivec-1).gt.thresh)
     1           converged = .false.
         enddo
c         cpu=-util_cpusec()
c         wall=-util_wallsec()
         if (.not.ma_pop_stack(l_residual))
     1        call errquit("cis.F: MA problem",102,MA_ERR)
         if (converged) then
ccx            call cis_xtidy
            if (nodezero) then
               write(LuOut,9240)
               call util_flush(LuOut)
            endif
            goto 200
         endif
      enddo
ccx      call errquit('cis.F: maxiter exceeded',iter,CALC_ERR)
      write(6,*)'CIS maxiter exceeded'
  200 continue
ccx
            call cis_xtidy
ccx
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
