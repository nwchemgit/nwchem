      subroutine ccsd_left(d_f1,d_lr1,d_lr2,d_t1,d_t2,d_v2,d_lambda1,
     1       d_lambda2,k_f1_offset,k_l1_offset,k_t1_offset,
     2       k_t2_offset,k_v2_offset,k_l2_offset,size_l1,size_l2,
     3       iter,nodezero,residual)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
      integer d_f1,d_lr1,d_lr2,d_t1,d_t2,d_v2,d_lambda1,d_lambda2
      integer k_f1_offset,k_l1_offset,k_t1_offset,k_t2_offset
      integer k_v2_offset,k_l2_offset,size_l1,size_l2
      character*255 filename
      double precision cpu, wall
      logical nodezero
c      integer iter,dummy
      integer dummy
      double precision r1, r2, residual
      double precision ddotfile
      external ddotfile
c
      call tce_diis_init()
      do iter=1,maxiter
         cpu=-util_cpusec()
         wall=-util_wallsec()
         if (nodezero.and.(iter.eq.1)) write(LuOut,9400) "CCSD Lambda"
         call tce_filename('lr1',filename)
         call createfile(filename,d_lr1,size_l1)
         call ccsd_lambda1(d_f1,d_lr1,d_t1,d_t2,d_v2,d_lambda1,
     1         d_lambda2,k_f1_offset,k_l1_offset,k_t1_offset,
     2         k_t2_offset,k_v2_offset,k_l1_offset,k_l2_offset)
c         print *,iter,'ccsd_lambda1'
         call reconcilefile(d_lr1,size_l1)
         call tce_filename('lr2',filename)
         call createfile(filename,d_lr2,size_l2)
         call ccsd_lambda2(d_f1,d_lr2,d_t1,d_t2,d_v2,d_lambda1,
     1         d_lambda2,k_f1_offset,k_l2_offset,k_t1_offset,
     2         k_t2_offset,k_v2_offset,k_l1_offset,k_l2_offset)
c         print *,iter,'ccsd_lambda2'
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
            call deletefile(d_lr2)
            call deletefile(d_lr1)
            call tce_diis_tidy()
            goto 2000
         endif
         call tce_diis(.true.,iter,.true.,.true.,.false.,.false.,
     1                  d_lr1,d_lambda1,k_l1_offset,size_l1,
     2                  d_lr2,d_lambda2,k_l2_offset,size_l2,
     3                  dummy,dummy,dummy,dummy,
     4                  dummy,dummy,dummy,dummy)
         call deletefile(d_lr2)
         call deletefile(d_lr1)
         if (nodezero) call util_flush(LuOut)
      enddo
      call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)

 2000 continue
 9400 format(/,1x,A,' iterations',/,
     1  1x,'--------------------------------------',/
     2  1x,'Iter          Residuum     Cpu    Wall',/
     3  1x,'--------------------------------------')
 9410 format(
     1  1x,'--------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f18.13,2f8.1)
      end
c $Id$
