      Subroutine dft_get2eri(CD_bas_han, g_2ceri,oskel)
c
C$Id: dft_get2eri.F,v 1.12 2001-12-28 19:38:53 edo Exp $
c
      implicit none
      integer CD_bas_han
      integer g_2ceri
c
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "sym.fh"
c      
      integer nshells_cd, MXnshbf_cd, me, nproc, icount,
     &        maxg, mscratch_2e2c, n1, n2, Nintegrals,
     &        ishp, ifirstp, ilastp, nshp,
     &        ishq, ifirstq, ilastq, nshq
      integer lbuf, ibuf, lscr, iscr, lvec, ivec
      integer nxtask
      logical oskel,odoit
      external nxtask
      double precision q2
c     
c     Compute the matrix of 2-ctr ERIs.
c     
      if ( .not. bas_numcont(CD_bas_han,nshells_cd) )then
         call errquit('Exiting in gaget2eri.',1)
      endif
      if ( .not. bas_nbf_cn_max(CD_bas_han,MXnshbf_cd) )then
         call errquit('Exiting in gaget2eri.',2)
      endif
c
c     Initialize ERI.
c
      call ga_sync
      call ga_zero(g_2ceri)
c
c     Allocate scratch space.
c
      me  = ga_nodeid()
      nproc  = ga_nnodes()
c     niter  = (nshells_cd*(nshells_cd+1))/2
      icount = (nshells_cd + mod(nshells_cd,2))/(2*nproc)
      icount = max(icount,1)
      call int_mem_2e2c(maxg, mscratch_2e2c)
c
c      write(6,*)' nproc, me, icount, maxg, mscratch_2e2c: ',
c     &            nproc, me, icount, maxg, mscratch_2e2c
c
      if (.not.MA_Push_Get(MT_Dbl,maxg,'buf',lbuf,ibuf))
     &   call errquit('dft_get2eri: cannot allocate buf',0)
      if (.not.MA_Push_Get(MT_Dbl,maxg,'vec',lvec,ivec))
     &   call errquit('dft_get2eri: cannot allocate buf',0)
      if (.not.MA_Push_Get(MT_Dbl,mscratch_2e2c,'scr',lscr,iscr))
     &   call errquit('dft_get2eri: cannot allocate scr',0)
c
      n1 = nxtask(nproc,icount)
      n2 = 0
c
      do 20 ishp = 1, nshells_cd
c
         if (.not. bas_cn2bfr(CD_bas_han, ishp, ifirstp, ilastp))
     &      call errquit('dft_get2eri: problem with bas_cn2bfr',0)
c
         nshp = ilastp - ifirstp + 1
c
         do 10 ishq = 1, ishp
c     
            odoit=.true.
            q2=1d0
            if (oskel.and.(.not.
     .           sym_shell_pair(cd_bas_han,ishp,ishq,q2))) goto 10
            if ( n2.eq.n1 )then
               if (.not. bas_cn2bfr(CD_bas_han, ishq, ifirstq, ilastq))
     &            call errquit('dft_get2eri: problem with bas_cn2bfr',0)
c
               nshq = ilastq - ifirstq + 1
c
               Nintegrals=maxg
c
               call int_2e2c(CD_bas_han, ishp, CD_bas_han, ishq,
     &                       mscratch_2e2c, Dbl_MB(iscr), Nintegrals,
     &                       Dbl_MB(ibuf))
               if(oskel) then
                  if(ishq.ne.ishp) q2=q2*.5d0
                  if(q2.gt.1d0) call dscal(
     .              nshp*nshq,q2,dbl_mb(ibuf),1)
               endif
c
c              Re-order.  
c
               call sca_lab_mat(CD_bas_han, ishp, ishq, nshp, nshq,
     &                          g_2ceri, DBL_MB(ibuf), DBL_MB(ivec),
     .              'put')
c
               n1 = nxtask(nproc,icount)
c
            endif
c
            n2 = n2 + 1
c
   10    continue
c
   20 continue
c
      if (.not.ma_chop_stack(lbuf))
     &   call errquit('dft_get2eri: cannot pop stack',0)
c
      n1 = nxtask(-nproc,icount)
c
      return
      end

