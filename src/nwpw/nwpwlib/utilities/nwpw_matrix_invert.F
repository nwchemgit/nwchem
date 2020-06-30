
      subroutine nwpw_matrix_invert(N,A,LDA)
      implicit none
      integer N,LDA
      real*8 A(LDA,*)

#include "bafdecls.fh"
#include "errquit.fh"

      integer B(2),ipiv(2)

      if (.not.BA_push_get(mt_dbl,LDA*N,'B',B(2),B(1)))
     > call errquit('nwpw_matrix_invert:out of stack',0,MA_ERR)
      if (.not.BA_push_get(mt_int,N,'ipiv',ipiv(2),ipiv(1)))
     > call errquit('nwpw_matrix_invert:out of stack',0,MA_ERR)

      call nwpw_matrix_invert_sub(N,A,LDA,dbl_mb(B(1)),LDA,
     >                            int_mb(ipiv(1)))

      if (.not.BA_pop_stack(ipiv(2)))
     > call errquit('nwpw_matrix_invert:pop stack',0,MA_ERR)
      if (.not.BA_pop_stack(B(2)))
     > call errquit('nwpw_matrix_invert:pop stack',0,MA_ERR)

      return
      end 

      subroutine nwpw_matrix_invert_sub(n,a,lda,b,ldb,ipiv)
      implicit none
      integer n,lda,ldb
      real*8 a(lda,*),b(ldb,*)
      integer ipiv(*)

      integer i,j,info

!$OMP MASTER
      do j=1,n
         b(j,j) = 1.0d0
         do i=j+1,n
            b(i,j) = 0.0d0
            b(j,i) = 0.0d0
         end do
      end do
      call dgesv(n,n,a,lda,ipiv,b,ldb,info)
      !write(*,*) "nwpw_matrix_invert_sub: info=",info
      do j=1,n
         a(j,j) = b(j,j)
         do i=j+1,n
            a(i,j) = b(i,j)
            a(j,i) = b(i,j)
         end do
      end do
!$OMP END MASTER


      return
      end
  
c $Id$
