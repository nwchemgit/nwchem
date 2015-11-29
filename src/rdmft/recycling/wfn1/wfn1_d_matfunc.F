#include "wfn1_ad.F"
#include "blas/double/wfn1_ad_daxpy.f"
#include "blas/double/wfn1_ad_dcopy.f"
#include "blas/double/wfn1_ad_ddot.f"
#include "blas/double/wfn1_ad_dgemm.F"
#include "blas/double/wfn1_ad_dger.f"
#include "blas/double/wfn1_ad_dgemv.f"
#include "blas/double/wfn1_ad_dnrm2.f"
#include "blas/double/wfn1_ad_dscal.f"
#include "blas/double/wfn1_ad_dswap.f"
#include "blas/double/wfn1_ad_dsymm.f"
#include "blas/double/wfn1_ad_dsymv.f"
#include "blas/double/wfn1_ad_dsyr2.f"
#include "blas/double/wfn1_ad_dsyr2k.f"
#include "blas/double/wfn1_ad_dsyrk.f"
#include "blas/double/wfn1_ad_dtrmm.f"
#include "blas/double/wfn1_ad_dtrmv.f"
#include "blas/double/wfn1_ad_dtrsm.f"
#include "blas/double/wfn1_ad_dtrsv.f"
#include "lapack/lapack_routine/wfn1_ad_iladlc.f"
#include "lapack/lapack_routine/wfn1_ad_iladlr.f"
#include "lapack/double/wfn1_ad_disnan.f"
#include "lapack/double/wfn1_ad_dlae2.f"
#include "lapack/double/wfn1_ad_dlaev2.f"
cinclude "lapack/double/wfn1_ad_dlaisnan.f"
#include "lapack/double/wfn1_ad_dlanst.f"
#include "lapack/double/wfn1_ad_dlansy.f"
#include "lapack/double/wfn1_ad_dlapy2.f"
#include "lapack/double/wfn1_ad_dlarf.f"
#include "lapack/double/wfn1_ad_dlarfb.f"
#include "lapack/double/wfn1_ad_dlarfg.f"
#include "lapack/double/wfn1_ad_dlarft.f"
#include "lapack/double/wfn1_ad_dlartg.f"
#include "lapack/double/wfn1_ad_dlascl.f"
#include "lapack/double/wfn1_ad_dlaset.f"
#include "lapack/double/wfn1_ad_dlasr.f"
#include "lapack/double/wfn1_ad_dlasrt.f"
#include "lapack/double/wfn1_ad_dlascale.f"
#include "lapack/double/wfn1_ad_dlassq.f"
#include "lapack/double/wfn1_ad_dlatrd.f"
#include "lapack/double/wfn1_ad_dorg2l.f"
#include "lapack/double/wfn1_ad_dorg2r.f"
#include "lapack/double/wfn1_ad_dorgql.f"
#include "lapack/double/wfn1_ad_dorgqr.f"
#include "lapack/double/wfn1_ad_dorgtr.f"
#include "lapack/double/wfn1_ad_dpotf2.f"
#include "lapack/double/wfn1_ad_dpotrf.f"
#include "lapack/double/wfn1_ad_dsteqr.f"
#include "lapack/double/wfn1_ad_dsterf.f"
#include "lapack/double/wfn1_ad_dsyev.f"
#include "lapack/double/wfn1_ad_dsygs2.f"
#include "lapack/double/wfn1_ad_dsygst.f"
#include "lapack/double/wfn1_ad_dsygv.f"
#include "lapack/double/wfn1_ad_dsytd2.f"
#include "lapack/double/wfn1_ad_dsytrd.f"
C>
C> \brief Power function
C>
      function wfn1_ad_pow(x)
      use wfn1_ad1
      implicit none
#include "wfn1_nxt.fh"
      type(wfn1_ad_dble) :: wfn1_ad_pow
      type(wfn1_ad_dble) :: x
c
      wfn1_ad_pow = x**wfn1_nxt_x_pow
c
      return
      end
C>
C> \brief Compute the gradient of a matrix function
C>
C> Compute the gradient of a matrix function using automatic
C> differentiation. The general procedure is to diagonalize the
C> matrix, apply the function to the eigenvalues, transform the result
C> back. Using automatic differentiation the elements of the result
C> will contain the matrix function values as well as their gradients.
C>
C> The result is a 4th order tensor. We arrange the elements of the
C> derivative as
C> \f{eqnarray*}{
C>   D_{klij} &=& \frac{\partial F_{kl}(X)}{\partial X_{ij}}
C> \}
C> i.e. for every pair \f$(i,j)\f$ a whole matrix is stored.
C>
      subroutine wfn1_d_matfunc(nbf,x,s,s12,sm12,ifunc,d)
      use wfn1_ad1
      implicit none
#include "lapack/double/intf_wfn1_ad_dsyev.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer nbf !< [Input] Then number of basis functions
      double precision x(nbf,nbf) !< [Input] The matrix \f$X\f$
      double precision s(nbf,nbf) !< [Input] The overlap matrix
      double precision s12(nbf,nbf) !< [Input] The sqrt(overlap) matrix
      double precision sm12(nbf,nbf) !< [Input] The 1/sqrt(overlap) matrix
      integer ifunc  !< [Input] The function identifier
      double precision d(nbf,nbf,nbf,nbf) !< [Output] The matrix
      !< function derivative
c
c     local
c
      type(wfn1_ad_dble) :: typemem(4)
      double precision   :: dblemem(2)
      integer nwords, nbytedble, nbytetype
      integer l_x,  k_x  !< Offset for matrix X
      integer l_s,  k_s  !< Offset for matrix S
      integer l_w,  k_w  !< Offset for vector W
      integer l_wk, k_wk !< Offset for vector Work
c
      integer info, lwork, i, j, ii
c
      type(wfn1_ad_dble), external :: wfn1_ad_pow
c
      nbytedble = loc(dblemem(2))-loc(dblemem(1))
      nbytetype = loc(typemem(2))-loc(typemem(1))
      nwords = (nbytetype+nbytedble-1)/nbytedble
c
      info  = 0
      lwork = -1
c     call wfn1_ad_dsygv(1,'V','U',nbf,typemem(4),nbf,typemem(3),nbf,
c    &                   typemem(2),typemem,lwork,info)
      call wfn1_ad_dsyev('V','U',nbf,typemem(3),nbf,
     &                   typemem(2),typemem,lwork,info)
      lwork = value(typemem(1))
      lwork = max(nbf*nbf,lwork)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nwords,"mat X",l_x,k_x))
     &  call errquit("wfn1_d_matfunc: alloc X failed",nbf*nbf*nwords,
     &               MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf*nwords,"mat S",l_s,k_s))
     &  call errquit("wfn1_d_matfunc: alloc S failed",nbf*nbf*nwords,
     &               MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nwords,"vec W",l_w,k_w))
     &  call errquit("wfn1_d_matfunc: alloc W failed",nbf*nwords,
     &               MA_ERR)
      if (.not.ma_push_get(MT_DBL,lwork*nwords,"vec WK",l_wk,k_wk))
     &  call errquit("wfn1_d_matfunc: alloc WK failed",lwork*nwords,
     &               MA_ERR)
c
      if (ifunc.eq.1) then
        call wfn1_d_matfunc_doit(nbf,lwork,x,s12,sm12,wfn1_ad_pow,d,
     &       dbl_mb(k_x),dbl_mb(k_s),dbl_mb(k_w),dbl_mb(k_wk))
      else
        call errquit("wfn1_d_matfunc: illegal func id: ",ifunc,UERR)
      endif
c
      if (.not.ma_pop_stack(l_wk))
     &  call errquit("wfn1_d_matfunc: dealloc WK failed",lwork*nwords,
     &               MA_ERR)
      if (.not.ma_pop_stack(l_w))
     &  call errquit("wfn1_d_matfunc: dealloc W failed",nbf*nwords,
     &               MA_ERR)
      if (.not.ma_pop_stack(l_s))
     &  call errquit("wfn1_d_matfunc: dealloc S failed",nbf*nbf*nwords,
     &               MA_ERR)
      if (.not.ma_pop_stack(l_x))
     &  call errquit("wfn1_d_matfunc: dealloc X failed",nbf*nbf*nwords,
     &               MA_ERR)
c
      end
C>
C> \brief actually do it
C>
      subroutine wfn1_d_matfunc_doit(nbf,lwork,x,s12,sm12,func,d,ad_x,
     &                               ad_s,ad_w,ad_wk)
      use wfn1_ad1
      implicit none
#include "blas/double/intf_wfn1_ad_dcopy.fh"
#include "blas/double/intf_wfn1_ad_dgemm.fh"
#include "lapack/double/intf_wfn1_ad_dsyev.fh"
#include "errquit.fh"
      integer nbf !< [Input] The number of basis functions
      integer lwork !< [Input] The length of the LAPACK work space
      double precision x(nbf,nbf) !< [Input] The matrix X
      double precision s12(nbf,nbf) !< [Input] The matrix sqrt(S)
      double precision sm12(nbf,nbf) !< [Input] The matrix 1/sqrt(S)
      type(wfn1_ad_dble), external :: func !< [Input] The function
      double precision d(nbf,nbf,nbf,nbf) !< [Output] The derivative
      !< tensor
      type(wfn1_ad_dble) :: ad_x(nbf,nbf) !< [Scratch]
      type(wfn1_ad_dble) :: ad_s(nbf,nbf) !< [Scratch]
      type(wfn1_ad_dble) :: ad_w(nbf)     !< [Scratch]
      type(wfn1_ad_dble) :: ad_wk(lwork)  !< [Scratch]
c
c     Constants
c
      type(wfn1_ad_dble) :: one
      type(wfn1_ad_dble) :: zero
c
c     local
c
      integer i, j, k, l, it, iv, nv
      double precision t
      integer info
c
c     Initialize the matrices
c
      one  = 1.0d0
      zero = 0.0d0
      do j = 1, nbf
        do i = 1, nbf
          ad_s(i,j) = s12(i,j)
          ad_x(i,j) = set_dij(x(i,j),nbf,i,j)
        enddo
      enddo
c
c     Transform the input matrix to S^{1/2}XS^{1/2}
c     I.e. into an orthogonal basis
c
      call wfn1_ad_dgemm('n','n',nbf,nbf,nbf,one,ad_x,nbf,ad_s,nbf,
     &                   zero,ad_wk,nbf)
      call wfn1_ad_dgemm('n','n',nbf,nbf,nbf,one,ad_s,nbf,ad_wk,nbf,
     &                   zero,ad_x,nbf)
c
c     Diagonalize the matrix X (we now use dsyev as in an orthogonal
c     basis the metric is the identity matrix).
c
      info = 0
      call wfn1_ad_dsyev('V','U',nbf,ad_x,nbf,ad_w,ad_wk,
     &                   lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_ad_dsygv: info.ne.0: ",info,UERR)
      endif
c
c     Transform the eigenvectors back into the original basis
c
      do j = 1, nbf
        do i = 1, nbf
          ad_s(i,j) = sm12(i,j)
        enddo
      enddo
      call wfn1_ad_dcopy(nbf*nbf,ad_x,1,ad_wk,1)
      call wfn1_ad_dgemm('n','n',nbf,nbf,nbf,one,ad_s,nbf,ad_wk,nbf,
     &                   zero,ad_x,nbf)
c
c     Apply the function to the eigenvalues
c
      do i = 1, nbf
        ad_w(i) = func(max(zero,ad_w(i)))
      enddo
c
c     Back transform the resulting diagonal matrix
c
      do j = 1, nbf
        do i = 1, nbf
          ad_s(i,j) = 0.0d0
        enddo
      enddo
      do j = 1, nbf
        do i = 1, nbf
          do k = 1, nbf
            ad_s(i,j) = ad_s(i,j) + ad_x(i,k)*ad_w(k)*ad_x(j,k)
          enddo
        enddo
      enddo
c
c     Extract tensor D
c
      do j = 1, nbf
        do i = 1, nbf
          do l = 1, nbf
            do k = 1, nbf
              d(k,l,i,j) = 0.0d0
            enddo
          enddo
        enddo
      enddo
      do l = 1, nbf
        do k = 1, nbf
          nv = get_nvar(ad_s(k,l))
          do iv = 1, nv
            call get_d1(ad_s(k,l),iv,t,it)
            j = (it-1)/nbf+1
            i = it-(j-1)*nbf
            d(k,l,i,j) = d(k,l,i,j) + 0.5d0*t
            d(k,l,j,i) = d(k,l,j,i) + 0.5d0*t
          enddo
        enddo
      enddo
c
c     Done
c
      end
      
