      Subroutine dft_mxspin_ovlp(basis, g_alpha, g_beta, g_tmp)

C$Id: dft_mxspin_ovlp.F,v 1.1 1997-01-30 01:43:00 d3h449 Exp $
      Implicit none
      integer basis
      integer g_s              ! overlap
      integer g_alpha          ! alpha eigenvecs [input]
      integer g_beta           ! beta eigenvecs [input]
      integer g_tmp            ! scratch space

#include "bas.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "stdio.fh"
c
c     local
c
      integer me,nproc
c
      integer i,j,n,m,mend
      integer k_ord,l_ord,k_tmp,l_tmp,k_scr,l_scr
c      
      integer at(10), bt(10)
      double precision prodmin, prodtmp, pt(10)
c
c*****************************************************************************
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in max_ovl',1)
c
      call gacrea(geom,basis,nbf,nbf,'ovl ',g_s,'atom')
c
      if(.not.MA_Push_Get(MT_Dbl,nbf,'tmpm',l_tmp, k_tmp))
     &     call errquit('dft_mxovl: cannot allocate tmpm',0)
c
      call ga_zero(g_s)
      call int_1e_ga(basis,basis,g_s,'overlap',.false.)
c
c     Compute matrix mult (C_alpha)T * S * C_beta = S`
c     
      call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_alpha,g_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_beta,0.d0,g_s)
c
      if(me.eq.0) then
        write(LuOut, 9998)
c        call ga_print(g_s)
      endif
      prodmin=0.1d0
      if(me.eq.0) then
        n = 0
        do i=1,nbf
c
c         get row of g_s 
c
          call ga_get(g_s,i,i,1,nbf,DBL_MB(k_tmp),1)
          do j=1,nbf
            prodtmp = abs(dbl_mb(k_tmp+j-1))
            if(prodtmp.gt.prodmin) then
              n = n + 1
              at(n) = i
              bt(n) = j
              pt(n) = prodtmp
              if(n.ge.10.or.(i.eq.nbf.and.j.eq.nbf))then
                mend = 10
c                if(n.lt.10)mend = n
                write(LuOut,9999)(at(m),m=1,mend),
     &                           (bt(m),m=1,mend),
     &                           (pt(m),m=1,mend)
                n = 0
                do m = 1, 10
                  at(m) = 0
                  bt(m) = 0
                  pt(m) = 0.0d0
                enddo
              endif
            endif
          enddo
        enddo
      endif
 9998 format(/,1x,'  alpha - beta orbital overlaps ',/,
     &         1x,'  ----------------------------- ',/)
 9999 format(1x,'  alpha ',10(1x,i5,1x),/
     &       1x,'   beta ',10(1x,i5,1x),/
     &       1x,'overlap ',10(f7.3),/)
c
      call gadest(g_s)
c
      if(.not.MA_Pop_Stack(l_tmp))
     & call errquit('dft_mxovl: cannot pop stack',0)
c
      return
      end
