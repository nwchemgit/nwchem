      Subroutine dft_invdiag(geom,basis,iga_A,iga_cdinv,n)

C$Id: dft_invdiag.F,v 1.1 1995-04-27 03:16:19 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      parameter (toll=1.d-6)
      Parameter (THRESHOLD=1.D-12)
      logical LResult
      integer geom,basis
#include <mafdecls.h>
#include <global.fh>
#include <tcgmsg.fh>


      me=ga_nodeid()
      nproc=ga_nnodes()
      call ga_sync
      call gacrea(geom,basis,
     &     n,n,'ga temp2',iga_tmp2,'row')
      LResult = MA_Push_Get(MT_Dbl,n,'evals',lev,iev)
      LResult = MA_Push_Get(MT_Dbl,n,'itmp',ltmp,itmp)
      
      
C     diag
      call ga_sync
#if defined(PARALLEL_DIAG)         
      call ga_diag_std(iga_A,iga_tmp2,DBl_MB(iev))
#else
      call ga_diag_std_seq(iga_A,iga_tmp2,DBl_MB(iev))
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
        call ga_get(iga_tmp2,
     *       1,n,i,i,DBL_MB(itmp),1)
        do j=0,n-1
          DBL_MB(itmp+j)=DBL_MB(itmp+j)*DBl_MB(iev+i-1)
        enddo
        call ga_put(iga_A,
     *       1,n,i,i,DBL_MB(itmp),1)
      enddo
      call ga_sync

C     (U * sigma^-1) * U(transp) 
      
      call ga_dgemm('N','T',n,n,n,1.d0,
     &     iga_A,iga_tmp2,0.d0,iga_cdinv)


      LResult = MA_Pop_Stack(ltmp)
      LResult = MA_Pop_Stack(lev)
      call ga_SYNC
      call gadest(iga_tmp2)
      return
      end
