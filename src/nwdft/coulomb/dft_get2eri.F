      Subroutine dft_get2eri(CD_bas_han, g_2ceri)
c
C$Id: dft_get2eri.F,v 1.9 2000-04-22 01:05:43 edo Exp $
c
      implicit none
      integer CD_bas_han
      integer g_2ceri
c
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c      
      integer nshells_cd, MXnshbf_cd, me, nproc, icount, MXint,
     &        maxg, mscratch_2e2c, n1, n2, Nintegrals,
     &        ishp, ifirstp, ilastp, nshp,
     &        ishq, ifirstq, ilastq, nshq
      integer lbuf, ibuf, lscr, iscr, lvec, ivec
      integer nxtask
      external nxtask
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
      MXint = MXnshbf_cd*MXnshbf_cd
c
c      write(6,*)' nproc, me, icount, maxg, mscratch_2e2c, MXint: ',
c     &            nproc, me, icount, maxg, mscratch_2e2c, MXint
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
            if ( n2.eq.n1 )then
               if (.not. bas_cn2bfr(CD_bas_han, ishq, ifirstq, ilastq))
     &            call errquit('dft_get2eri: problem with bas_cn2bfr',0)
c
               nshq = ilastq - ifirstq + 1
c
               Nintegrals=MXint
c
               call int_2e2c(CD_bas_han, ishp, CD_bas_han, ishq,
     &                       mscratch_2e2c, Dbl_MB(iscr), Nintegrals,
     &                       Dbl_MB(ibuf))
c
c              Re-order.  
c
               call sca_lab_mat(CD_bas_han, ishp, ishq, nshp, nshq,
     &                          g_2ceri, DBL_MB(ibuf), DBL_MB(ivec))
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
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dft_get2eri: cannot pop stack',0)
      if (.not.ma_pop_stack(lvec))
     &   call errquit('dft_get2eri: cannot pop stack',0)
      if (.not.ma_pop_stack(lbuf))
     &   call errquit('dft_get2eri: cannot pop stack',0)
c
      n1 = nxtask(-nproc,icount)
c
      call ga_sync
c      call ga_print(g_2ceri)
c      call ga_sync
c
      return
      end

