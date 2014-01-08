C> \ingroup util_ga
C> @{
C>
C> \brief Anti-symmetrizes the matrix \f$A\f$
C>
C> Anti-symmetrizes the matrix \f$A\f$ by computing
C> \f{eqnarray*}{
C>    A = \left(A-A^T\right)/2
C> \f}
C> 
      subroutine ga_antisymmetrize(g_a) 
c
c     Antisymmetrizes matrix A:  A := .5 * (A-A`)
c     copied from ga_symmetrize, A+A' replaced by A-A'
c     J. Nieplocha 08.22.93
c
      implicit none
#include "mafdecls.fh" 
#include "global.fh"
      integer g_a               ! Matrix to symmetrize 
c
      Double Precision alpha   
      Integer myproc, proc 
      Integer ilo, ihi, jlo, jhi, myelem, nrow,ncol, lda
      Integer dim1, dim2, type
      Logical have_data

      Integer     adrA,handa          !A
      Integer hB, adrB          !B
c
c***  check environment
c
      myproc = ga_nodeid() 
      proc   = ga_nnodes()
      alpha = 0.5d0
c  
      call ga_inquire(g_a,  type, dim1, dim2)
c
      if (dim1 .ne. dim2) then
         call ga_error('ga_sym: can only sym square matrix', 0)
      endif
c
c***  Find the local distribution
      call ga_distribution(g_a, myproc, ilo, ihi, jlo, jhi)
#ifdef  BAD_GACCESS
      if(.not.ma_push_get(MT_DBL,(ihi-ilo+1)*(jhi-jlo+1),
     $  'scratch buff', handa, adra)) call 
     $  errquit('gaantisymm: pushget failed',0)
#endif
c
      have_data = ihi.gt.0 .and. jhi.gt.0
c
c     if(myproc.eq.0)call ga_print(g_a)
      call ga_sync()
      if (have_data) then
#ifdef  BAD_GACCESS
         lda=ihi-ilo+1
         call ga_get(g_a, ilo, ihi, jlo, jhi, dbl_mb(adrA), lda)
#else
         call ga_access(g_a, ilo, ihi, jlo, jhi, adrA, lda)
#endif
         nrow = ihi-ilo+1
         ncol = jhi-jlo+1
         myelem = nrow * ncol
         if (.not. ma_push_get(MT_DBL,myelem,'ga_symtr:B',hB,adrB))
     $      call ga_error('ga_symmetrize: insufficient memory', myelem)
c        
         call GA_GET(g_a, jlo, jhi, ilo, ihi, dbl_mb(adrB),ncol) 
c        print *, myproc,'a, aT', dbl_mb(adrA), dbl_mb(adrB),
c    $            ilo, ihi, jlo, jhi
      endif

      call ga_sync()
      
      if(have_data) then
         call gai_subtr(nrow, ncol, dbl_mb(adrA), lda, dbl_mb(adrB),
     &                ncol, alpha)
#ifdef  BAD_GACCESS
         call ga_put(g_a, ilo, ihi, jlo, jhi, dbl_mb(adrA), lda)
#else
         call ga_release_update(g_a, ilo, ihi, jlo, jhi)
#endif
         if (.not. ma_pop_stack(hB))
     $        call ga_error('ga_symmetrize: ma_pop?', 0)
      endif
#ifdef  BAD_GACCESS
      if(.not.ma_pop_stack(handa)) call 
     $  errquit('gaantisymm: popstack failed',0)
#endif
c
      call ga_sync()
c
      end


      subroutine nga_antisymmetrize(g_a) 
c
c     Antisymmetrizes matrix A:  A := .5 * (A-A`)
c     This is silly hack at the moment that loops over each of the
c     first indices, gets a matrix, symmetrizes, and puts it back.
c     This is limited to 3 dimensions.  Eventually, we should be able
c     to grab much larger chunks at a time.
c 
      implicit none
#include "mafdecls.fh" 
#include "global.fh"
      integer g_a               ! Matrix to symmetrize 
c
      Double Precision alpha   
      Integer myproc, proc 
      Integer flo, fhi,ilo, ihi, jlo, jhi, myelem, nrow,ncol, lda(2)
      Integer lo(3), hi(3), lob(3), hib(3), ldb(2)
      Integer dim1, dim2, type, nmat
      Integer ndim, dims(3)
      Logical have_data

      Integer     adrA          !A
      Integer hB, adrB          !B
c
c***  check environment
c
      myproc = ga_nodeid() 
      proc   = ga_nnodes()
      alpha = 0.5d0
c  
      ndim = ga_ndim(g_a)
      if (ndim.ne.3) 
     $   call ga_error('nga_sym: needs 3 dims',ndim)
      call nga_inquire(g_a,  type, ndim, dims)
      dim1 = dims(2)
      dim2 = dims(3)
c
      if (dim1 .ne. dim2) then
         call ga_error('nga_sym: can only sym square matrix', 0)
      endif
c
c***  Find the local distribution
      call nga_distribution(g_a, myproc, lo, hi)
c
c Take care of some indexing
c
      flo = lo(1)
      ilo = lo(2)
      jlo = lo(3)
      fhi = hi(1)
      ihi = hi(2)
      jhi = hi(3)
      lob(1) = flo
      hib(1) = fhi
      lob(2) = jlo  ! switch i and j
      hib(2) = jhi
      lob(3) = ilo
      hib(3) = ihi
      nmat = fhi-flo+1
      nrow = ihi-ilo+1
      ncol = jhi-jlo+1
      myelem = nrow * ncol * nmat
      ldb(1) = nmat
      ldb(2) = ncol
      have_data = fhi.gt.0 .and. ihi.gt.0 .and. jhi.gt.0
c
c Get memory (if we have data)
c
      if (have_data) then
         if (.not. ma_push_get(MT_DBL,myelem,'ga_symtr:B',hB,adrB))
     $      call ga_error('nga_symmetrize: insufficient memory', 
     $                     myelem)
      endif
c
c Get the differenc pieces
c
      call ga_sync()
      if (have_data) then
         call nga_access(g_a, lo, hi, adrA, lda)
c        
         call nga_get(g_a, lob, hib, dbl_mb(adrB),ldb) 
      endif

      call ga_sync()
      
      if(have_data) then
         call ngai_subtr(nmat, nrow, ncol, dbl_mb(adrA), lda(1), lda(2), 
     &                 dbl_mb(adrB), ldb(1), ldb(2), alpha)
         call nga_release_update(g_a, lo, hi)
      endif
c
      call ga_sync()
c
      if (have_data) then
         if (.not. ma_pop_stack(hB))
     $        call ga_error('ga_symmetrize: ma_pop?', 0)
      endif
c
      end


      subroutine gai_subtr(nrow, ncol, A, lda, B, ldb, alpha)
c
c     A := alpha * (A - B^)
c
      implicit none
      Integer nrow, ncol, lda, ldb
      Double Precision A(lda,ncol), B(ldb,nrow), alpha
      Integer i,j
c
      do j = 1, ncol
         do i = 1, nrow
            a(i,j) = alpha*(a(i,j)-b(j,i))
         enddo
      enddo
      END


      subroutine ngai_subtr(nmat, nrow, ncol, A, lda1, lda2, 
     $                    B, ldb1, ldb2, alpha)
c
c     A := alpha * (A - B^)
c     Again, assumes 3 dimensions
c
      implicit none
      Integer nmat, nrow, ncol, lda1, lda2, ldb1, ldb2
      Double Precision A(lda1,lda2,ncol), B(ldb1,ldb2,nrow), alpha
      Integer i,j,f
c
      do j = 1, ncol
         do i = 1, nrow
            do f = 1, nmat
               a(f,i,j) = alpha*(a(f,i,j)-b(f,j,i))
            enddo
         enddo
      enddo
      END
C>
C> @}
c $Id$
