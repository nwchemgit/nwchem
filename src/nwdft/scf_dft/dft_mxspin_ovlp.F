      Subroutine dft_mxspin_ovlp(basis, g_alpha, g_beta, g_tmp)

C$Id: dft_mxspin_ovlp.F,v 1.6 2000-10-25 22:50:14 windus Exp $
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
      integer i,j,jbig,n,ichunks,nbe,nend
      integer k_tmpr1, l_tmpr1, k_tmpr2, l_tmpr2, k_tmpi1, l_tmpi1
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision prodbig, prodtmp
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in max_ovl',1)
c
c      call gacrea(geom,basis,nbf,nbf,'ovl ',g_s,'atom')
      g_s = ga_create_atom_blocked(geom, basis, 'AO ovl')
c
      if(.not.MA_Push_Get(MT_Dbl,nbf,'real_tmp1',l_tmpr1, k_tmpr1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp1',0)
      if(.not.MA_Push_Get(MT_Dbl,nbf,'real_tmp2',l_tmpr2, k_tmpr2))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp2',0)
      if(.not.MA_Push_Get(MT_Int,nbf,'int_tmp1',l_tmpi1, k_tmpi1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate int_tmp1',0)
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
        write(LuOut, 9996)
c        call ga_print(g_s)
      endif
      if(me.eq.0) then
        jbig = 1         ! take care of compiler warnings
        do i = 1, nbf
c
c         get row of g_s 
c
          call ga_get(g_s,i,i,1,nbf,DBL_MB(k_tmpr1),1)
          prodbig=0.0d0
          do j = 1, nbf
            prodtmp = abs(dbl_mb(k_tmpr1+j-1))
            if(prodtmp.gt.prodbig) then
              prodbig = prodtmp
              jbig = j
            endif
          enddo
          dbl_mb(k_tmpr2+i-1) = prodbig
          int_mb(k_tmpi1+i-1) = jbig
        enddo
        ichunks = nbf/10
        if(ichunks*10.lt.nbf)ichunks = ichunks +1
        do i = 1, ichunks
          nbe = 10*(i-1) + 1
          nend = nbe + 9
          if(nend.gt.nbf)nend = nbf
          write(LuOut,9997)(n,n=nbe,nend)
          write(LuOut,9998)(int_mb(k_tmpi1+n-1),n=nbe,nend)
          write(LuOut,9999)(dbl_mb(k_tmpr2+n-1),n=nbe,nend)
        enddo
      endif
 9996 format(/,1x,'  alpha - beta orbital overlaps ',/,
     &         1x,'  ----------------------------- ',/)
 9997 format(/,1x,'  alpha ',10(1x,i5,1x))
 9998 format(  1x,'   beta ',10(1x,i5,1x))
 9999 format(  1x,'overlap ',10(f7.3),/)
c
c      call gadest(g_s)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_s', 0)
c
      if(.not.MA_Pop_Stack(l_tmpi1))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
      if(.not.MA_Pop_Stack(l_tmpr2))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
      if(.not.MA_Pop_Stack(l_tmpr1))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
c
      return
      end
