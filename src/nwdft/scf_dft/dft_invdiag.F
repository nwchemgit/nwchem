      Subroutine dft_invdiag(g_A, g_cdinv, n)

C$Id$
      implicit none
#include "errquit.fh"
c      
      integer g_a   ! [input]
      integer g_cdinv ! [output]
      integer n
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c
      integer me,nproc,i,j,g_tmp2
      integer lev,iev,ltmp,itmp
      double precision toll,THRESHOLD
      parameter (toll=1.d-6,THRESHOLD=1.D-12)
c
      me=ga_nodeid()
      nproc=ga_nnodes()
      call ga_sync
*ga:1:0
      if (.not. ga_create(mt_dbl, n, n, 'ga_temp2', n, 0, g_tmp2)) 
     &   call errquit('dft_main0d: error creating g_xcinv',0, GA_ERR)
c
      if(.not.MA_Push_Get(MT_Dbl,n,'evals',lev,iev))
     &  call errquit('dft_invdiag: cannot allocate evals',0, MA_ERR)
      if(.not.MA_Push_Get(MT_Dbl,n,'itmp',ltmp,itmp))
     &  call errquit('dft_invdiag: cannot allocate itmp',0, MA_ERR)
      
C     diag
      call ga_sync
#if defined(PARALLEL_DIAG)         
#ifdef SCALAPACK
      call dft_scaleig(g_a, g_tmp2, DBl_MB(iev),
     P     0) 
#else
      call ga_diag_std(g_A,g_tmp2,DBl_MB(iev))
#endif
#else
      call ga_diag_std_seq(g_A,g_tmp2,DBl_MB(iev))
#endif

C     check on eigenvalues

      do i=0,n-1
        if(abs(DBl_MB(iev+i)).lt.toll) then
          if(me.eq.0) write(LuOut,*) ' GAFACT - singular eigenvalue',i
          call util_flush(LuOut)
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
     &  call errquit('dft_invdiag: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lev))
     &  call errquit('dft_invdiag: cannot pop stack',0, MA_ERR)
      call ga_SYNC
c      call gadest(g_tmp2)
      if (.not. ga_destroy(g_tmp2)) call errquit
     &   ('dft_invdiag: could not destroy g_tmp2', 0, GA_ERR)
      return
      end
