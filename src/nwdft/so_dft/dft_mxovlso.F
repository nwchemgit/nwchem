      Subroutine dft_mxovlso(basis, evals, g_evecs, g_oevecs,g_tmp)

C$Id$
      Implicit none
#include "errquit.fh"
      integer basis
      integer g_s              ! overlap       [input] overwritten!
      integer g_evecs          ! new eigenvecs [input/output]
      integer g_oevecs         ! old eigenvecs [input]
      integer g_tmp            ! scratch space
      double precision evals(*)! eigenvalues  [input/output]

#include "bas.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
c
c     local
c
      integer me,nproc
c      integer nbf,i,j,ix
      integer i,j,ix
      integer k_ord,l_ord,k_tmp,l_tmp,k_scr,l_scr
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision prodmax
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in max_ovl',1, BASIS_ERR)
c
c      call gacrea(geom,basis,nbf,nbf,'ovl ',g_s,'atom')
      g_s = ga_create_atom_blocked(geom, basis, 'AO ovl')
c
      if(.not.MA_Push_Get(MT_Dbl,nbf,'tmpm',l_tmp, k_tmp))
     &     call errquit('dft_mxovl: cannot allocate tmpm',0, MA_ERR)
      if(.not.MA_Push_Get(MT_Dbl,nbf,'scr',l_scr, k_scr))
     &     call errquit('dft_mxovl: cannot allocate scr',0, MA_ERR)
      if(.not.MA_Push_Get(MT_int,nbf,'ord',l_ord, k_ord))
     &     call errquit('dft_mxovl: cannot allocate ord',0, MA_ERR)

      call ga_zero(g_s)
      call int_1e_ga(basis,basis,g_s,'overlap',.false.)
      call ifill(nbf, 0, int_mb(k_ord), 1)

C
C     compute matrix mult (V_old)T * S * V_new = S`
C     
      call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_oevecs,g_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_evecs,0.d0,g_s)
c
c      call ga_print(g_s)
c
      if(me.eq.0) then
        do i=1,nbf
          prodmax=-30.d0
C
C         get row of g_s 
C
          call ga_get(g_s,i,i,1,nbf,DBL_MB(k_tmp),1)
          do j=1,nbf
            if(abs(dbl_mb(k_tmp+j-1)).gt.prodmax) then
              prodmax=abs(dbl_mb(k_tmp+j-1))
              ix=j
            endif
          enddo

          int_mb(k_ord+i-1) = ix
C
C         invalidate column that you have taken
C
          call dfill(nbf, 0.0d0, dbl_mb(k_tmp), 1)
          call put_col(g_s,nbf,ix,DBL_MB(k_tmp))
        enddo
      endif
      call ga_igop(Msg_mxovl_1,int_mb(k_ord),nbf,'+')
c      if(me.eq.0) then
c        write(*,*)' orbital order: ',(int_mb(k_ord+i-1),i=1,nbf)
c      endif
c
      call ga_copy(g_evecs,g_oevecs)
      call dcopy(nbf,evals,1,dbl_mb(k_tmp),1)
      call dfill(nbf,0.d0,evals,1)
      do i=me+1,nbf,nproc
        ix=int_mb(k_ord+i-1)
        call get_col(g_oevecs,nbf,ix,DBL_MB(k_scr))
        evals(i)=dbl_mb(k_tmp+ix-1)
        call put_col(g_evecs,nbf,i,DBL_MB(k_scr))
      enddo
c
c      call gadest(g_s)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('dft_mxovl: could not destroy g_s', 0, GA_ERR)
c
      if(.not.MA_Pop_Stack(l_ord))
     & call errquit('dft_mxovl: cannot pop stack',0, MA_ERR)
      if(.not.MA_Pop_Stack(l_scr))
     & call errquit('dft_mxovl: cannot pop stack',0, MA_ERR)
      if(.not.MA_Pop_Stack(l_tmp))
     & call errquit('dft_mxovl: cannot pop stack',0, MA_ERR)
c
      call ga_dgop(Msg_mxovl_2,evals,nbf,'+')
c
      return
      end
