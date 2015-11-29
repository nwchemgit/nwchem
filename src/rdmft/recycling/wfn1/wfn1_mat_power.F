C> \ingroup wfn1
C> @{
C>
C> \brief Evaluate the power of a matrix
C>
C> This subroutine calculates \f$A' = A^{p}\f$ where \f$A\f$ is a
C> square matrix and \f$p\f$ is a floating point number. The function
C> evaluates the expression in place. I.e. the input matrix is
C> overwritten with the result.
C>
      subroutine wfn1_mat_power(n,a,lda,s,lds,p)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer n !< [Input] The dimension of the matrix
      integer lda !< [Input] The leading dimension of array \f$A\f$.
      integer lds !< [Input] The leading dimension of array \f$S\f$.
c
      double precision a(lda,n) !< [In/Output] The matrix
      !< - On input the matrix to calculate the power off
      !< - On output the input matrix to the power \f$p\f$
      double precision s(lds,n) !< [Input] The metric matrix
      double precision p        !< [Input] The power
c
      integer l_matvec, k_matvec
      integer l_diag, k_diag
      integer l_metric, k_metric
      integer l_work, k_work
      integer lwork, info
      integer itype
      integer i, j, k
      double precision size
c
      return ! LAPACK dsygv does not work!
c
      itype = 1
      if (.not.ma_push_get(MT_DBL,n*n,"matvec",l_matvec,k_matvec))
     &  call errquit("wfn1_mat_power: could not allocate matvec",0,
     &               MA_ERR)
      if (.not.ma_push_get(MT_DBL,n*n,"metric",l_metric,k_metric))
     &  call errquit("wfn1_mat_power: could not allocate metric",0,
     &               MA_ERR)
      if (.not.ma_push_get(MT_DBL,n,"diag",l_diag,k_diag))
     &  call errquit("wfn1_mat_power: could not allocate diag",0,
     &               MA_ERR)
c
      do j = 1, n
        do i = 1, n
          dbl_mb(k_matvec+(j-1)*n+i-1) = 0.0d0
          dbl_mb(k_metric+(j-1)*n+i-1) = 0.0d0
        enddo
        dbl_mb(k_diag+j-1) = 0.0d0
      enddo
      do j = 1, n
        do i = 1, j
          dbl_mb(k_matvec+(j-1)*n+i-1) = a(i,j)
          dbl_mb(k_metric+(j-1)*n+i-1) = s(i,j)
        enddo
      enddo
c
      info = 0
      lwork = -1
      call dsygv(itype,'V','U',n,dbl_mb(k_matvec),n,dbl_mb(k_metric),n,
     &           dbl_mb(k_diag),size,lwork,info)
      if (info.ne.0) call errquit("wfn1_mat_power w: info non-zero",
     &                            info,UERR)
      lwork = size+0.5d0
cDEBUG
c     lwork = 3*n-1
      write(*,*)'wfn1_mat_power: lwork, n: ',lwork,n,size
cDEBUG
      if (.not.ma_push_get(MT_DBL,lwork,"work",l_work,k_work))
     &  call errquit("wfn1_mat_power: could not allocate work",0,
     &               MA_ERR)
      info = 0
      itype = 1
      call dsygv(itype,'V','U',n,dbl_mb(k_matvec),n,dbl_mb(k_metric),n,
     &           dbl_mb(k_diag),dbl_mb(k_work),lwork,info)
      if (info.ne.0) call errquit("wfn1_mat_power: info non-zero",
     &                            info,UERR)
cDEBUG
      do i = 0, n-1
        write(*,*)'wfn1_mat_power: diag: ',i+1,dbl_mb(k_diag+i)
      enddo
      write(*,*)'wfn1_mat_power: vectors'
      call hess_hssout(dbl_mb(k_matvec),n,n,n)
cDEBUG
      do i = 0, n-1
c       dbl_mb(k_diag+i) = max(dbl_mb(k_diag+i),0.0d0)**p
        dbl_mb(k_diag+i) = max(dbl_mb(k_diag+i),0.0d0)
      enddo
      do j = 1, n
        do i = 1, n
          a(i,j) = 0.0d0
          do k = 1, n
            a(i,j) = a(i,j) + dbl_mb(k_diag+k-1)
     &                      * dbl_mb(k_matvec+(k-1)*n+i-1)
     &                      * dbl_mb(k_matvec+(k-1)*n+j-1)
          enddo ! k
        enddo ! i
      enddo ! j
      if (.not.ma_pop_stack(l_work))
     &  call errquit("wfn1_mat_power: could not deallocate work",
     &               0,MA_ERR)
      if (.not.ma_pop_stack(l_diag))
     &  call errquit("wfn1_mat_power: could not deallocate diag",
     &               0,MA_ERR)
      if (.not.ma_pop_stack(l_metric))
     &  call errquit("wfn1_mat_power: could not deallocate metric",
     &               0,MA_ERR)
      if (.not.ma_pop_stack(l_matvec))
     &  call errquit("wfn1_mat_power: could not deallocate matvec",
     &               0,MA_ERR)
      end
C>
C> @}
