      Subroutine dft_invdiag(geom,basis,g_A,g_cdinv,n)

C$Id: dft_invdiag.F,v 1.3 1995-12-20 03:27:05 og845 Exp $
      implicit none
      
      integer geom,basis
      integer g_a   ! [input]
      integer g_cdinv ! [output]
      integer n
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
C
      logical LResult
      integer me,nproc,i,j,g_tmp2
      integer lev,iev,ltmp,itmp
      double precision toll,THRESHOLD
      parameter (toll=1.d-6,THRESHOLD=1.D-12)
C
      me=ga_nodeid()
      nproc=ga_nnodes()
      call ga_sync
      call gacrea(geom,basis,
     &     n,n,'ga temp2',g_tmp2,'row')
      LResult = MA_Push_Get(MT_Dbl,n,'evals',lev,iev)
      LResult = MA_Push_Get(MT_Dbl,n,'itmp',ltmp,itmp)
      
      
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
        call ga_get(g_tmp2,
     *       1,n,i,i,DBL_MB(itmp),1)
        do j=0,n-1
          DBL_MB(itmp+j)=DBL_MB(itmp+j)*DBl_MB(iev+i-1)
        enddo
        call ga_put(g_A,
     *       1,n,i,i,DBL_MB(itmp),1)
      enddo
      call ga_sync

C     (U * sigma^-1) * U(transp) 
      
      call ga_dgemm('N','T',n,n,n,1.d0,
     &     g_A,g_tmp2,0.d0,g_cdinv)


      LResult = MA_Pop_Stack(ltmp)
      LResult = MA_Pop_Stack(lev)
      call ga_SYNC
      call gadest(g_tmp2)
      return
      end
