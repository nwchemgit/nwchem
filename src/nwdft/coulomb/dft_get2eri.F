      Subroutine dft_get2eri(CD_bas_han,g_2ceri,geom)

C$Id: dft_get2eri.F,v 1.5 1996-10-18 21:49:59 d3h449 Exp $
      implicit none
      integer CD_bas_han, geom
      integer g_2ceri
c**** 
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c      
      integer nshells_cd,MXnshbf_cd,me,nproc,icount,MXint,
     ,     maxg,mscratch_2e2c,n1,n2,Nint,
     ,     ishp,ifirstp,ilastp,nshp,
     ,     ishq,ifirstq,ilastq,nshq
      integer lbuf,ibuf,lscr,iscr,lvec,ivec
      integer nxtask
      external nxtask


c--   > Electron Repulsion Integrals

c     
c     Compute the matrix of 2-ctr ERIs.
c     
c******************************************************************************
      if( .not. bas_numcont(CD_bas_han,nshells_cd) )then
        call errquit('Exiting in gaget2eri.',1)
      end if
      if( .not. bas_nbf_cn_max(CD_bas_han,MXnshbf_cd) )then
        call errquit('Exiting in gaget2eri.',2)
      end if
c     Initialize ERI.
      call ga_sync
      call ga_zero(g_2ceri)
c     Allocate scratch space.

      me  = ga_nodeid()
      nproc  = ga_nnodes()
c     niter  = (nshells_cd*(nshells_cd+1))/2
      icount = (nshells_cd + mod(nshells_cd,2))/2
      icount = max(icount,1)

      call int_mem_2e2c(maxg, mscratch_2e2c)

      MXint=MXnshbf_cd*MXnshbf_cd

      if(.not.MA_Push_Get(MT_Dbl,MXint,'buf',lbuf,ibuf))
     &  call errquit('dft_get2eri: cannot allocate buf',0)
      if(.not.MA_Push_Get(MT_Dbl,MXint,'vec',lvec,ivec))
     &  call errquit('dft_get2eri: cannot allocate buf',0)
      if(.not.MA_Push_Get(MT_Dbl,mscratch_2e2c,'scr',lscr,iscr))
     &  call errquit('dft_get2eri: cannot allocate scr',0)

      n1 = nxtask(nproc,icount)
      n2 = 0

      do 20 ishp = 1,nshells_cd

        if (.not. bas_cn2bfr( CD_bas_han,ishp,
     &       ifirstp,ilastp))
     &       call errquit('Exiting in gaget2eri.',1)
            nshp=ilastp-ifirstp+1

        do 10 ishq = 1,ishp

          if( n2.eq.n1 )then
            if (.not. bas_cn2bfr( CD_bas_han,ishq,
     &           ifirstq,ilastq))
     &           call errquit('Exiting in gaget2eri.',2)
            nshq=ilastq-ifirstq+1

            Nint=MXint

            call int_2e2c(CD_bas_han, ishp, CD_bas_han, ishq,
     &                    mscratch_2e2c, Dbl_MB(iscr), Nint,
     &                    Dbl_MB(ibuf))


c           Re-order.  
            call sca_lab_mat( CD_bas_han,ishp,ishq,nshp,nshq,
     &                        g_2ceri,DBL_MB(ibuf),DBL_MB(ivec))

            n1 = nxtask(nproc,icount)

          end if

          n2 = n2 + 1

   10   continue
   20 continue

      if(.not.ma_pop_stack(lscr))
     &  call errquit('dft_get2eri: cannot pop stack',0)
      if(.not.ma_pop_stack(lvec))
     &  call errquit('dft_get2eri: cannot pop stack',0)
      if(.not.ma_pop_stack(lbuf))
     &  call errquit('dft_get2eri: cannot pop stack',0)

      n1 = nxtask(-nproc,icount)

      return
      end

