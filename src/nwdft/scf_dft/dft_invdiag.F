      Subroutine dft_invdiag(geom,basis,g_A,g_cdinv,n)

C$Id: dft_invdiag.F,v 1.5 1996-10-18 21:50:06 d3h449 Exp $
      implicit none
c      
      integer geom,basis
      integer g_a   ! [input]
      integer g_cdinv ! [output]
      integer n
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      integer me,nproc,i,j,g_tmp2
      integer lev,iev,ltmp,itmp
      double precision toll,THRESHOLD
      parameter (toll=1.d-6,THRESHOLD=1.D-12)
c
      me=ga_nodeid()
      nproc=ga_nnodes()
      call ga_sync
      call gacrea(geom,basis,n,n,'ga temp2',g_tmp2,'row')

      if(.not.MA_Push_Get(MT_Dbl,n,'evals',lev,iev))
     &  call errquit('dft_invdiag: cannot allocate evals',0)
      if(.not.MA_Push_Get(MT_Dbl,n,'itmp',ltmp,itmp))
     &  call errquit('dft_invdiag: cannot allocate itmp',0)
      
C     diag
      call ga_sync
#if defined(PARALLEL_DIAG)         
      call ga_diag_std(g_A,g_tmp2,DBl_MB(iev))
#else
      call ga_diag_std_seq(g_A,g_tmp2,DBl_MB(iev))
#endif

C     check on eigenvalues

      do i=0,n-1
        if(abs(DBl_MB(iev+i)).lt.toll) then
          write(0,*) ' GAFACT - singular eigenvalue',i
          call flush(0)
          DBl_MB(iev+i)=0.d0
        else
          DBl_MB(iev+i)=1.d0/DBl_MB(iev+i)
        endif 
      enddo

C     (U * sigma^-1)

      do i=me+1,n,nproc
        call ga_get(g_tmp2,1,n,i,i,DBL_MB(itmp),1)
        do j=0,n-1
          DBL_MB(itmp+j)=DBL_MB(itmp+j)*DBl_MB(iev+i-1)
        enddo
        call ga_put(g_A,1,n,i,i,DBL_MB(itmp),1)
      enddo

      call ga_sync

C     (U * sigma^-1) * U(transp) 
      
      call ga_dgemm('N','T',n,n,n,1.d0,g_A,g_tmp2,0.d0,g_cdinv)

      if(.not.ma_pop_stack(ltmp))
     &  call errquit('dft_invdiag: cannot pop stack',0)
      if(.not.ma_pop_stack(lev))
     &  call errquit('dft_invdiag: cannot pop stack',0)
      call ga_SYNC
      call gadest(g_tmp2)
      return
      end
