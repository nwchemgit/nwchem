      Subroutine dft_get2eri(CD_bas_han,iga_2ceri,geom)

C$Id: dft_get2eri.F,v 1.1 1995-04-26 03:43:36 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical LResult

c****
c**** nwchem handles
c****
      integer CD_bas_han, geom
c**** 
#include <bas.fh>

#include <mafdecls.h>
#include "global.fh"
#include "tcgmsg.fh"

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
      call ga_zero(iga_2ceri)
c     Allocate scratch space.

      me  = ga_nodeid()
      nproc  = ga_nnodes()
c     niter  = (nshells_cd*(nshells_cd+1))/2
      icount = (nshells_cd + mod(nshells_cd,2))/2
      icount = max(icount,1)

      call int_mem_2e2c(maxg, mscratch_2e2c)

      MXint=MXnshbf_cd*MXnshbf_cd

      LResult = MA_Push_Get(MT_Dbl,MXint,'buf',lbuf,ibuf)
      LResult = MA_Push_Get(MT_Dbl,MXint,'buf',lvec,ivec)
      LResult = MA_Push_Get(MT_Dbl,mscratch_2e2c,'scr',lscr,iscr)

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
     &                        iga_2ceri,DBL_MB(ibuf),DBL_MB(ivec))

            n1 = nxtask(nproc,icount)

          end if

          n2 = n2 + 1

   10   continue
   20 continue

      LResult = MA_Pop_Stack(lscr)
      LResult = MA_Pop_Stack(lvec)
      LResult = MA_Pop_Stack(lbuf)

      n1 = nxtask(-nproc,icount)

      return
      end

