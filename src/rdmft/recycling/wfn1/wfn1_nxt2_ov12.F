C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Compute \f$S^{-1/2}\f$
C>
C> We need \$S^{-1/2}\f$ to compute derivatives of functions of the
C> density matrix. Here we compute and store this quantity for later
C> use.
C>
      subroutine wfn1_nxt2_ov12(nbf,ov,ov12,ovm12)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer          nbf           !< [Input] The number of basis
                                     !< functions
      double precision ov(nbf,nbf)   !< [Input] The overlap matrix
                                     !< \f$S\f$
      double precision ov12(nbf,nbf) !< [Output] The matrix
                                     !< \f$S^{1/2}\f$
      double precision ovm12(nbf,nbf) !< [Output] The matrix
                                      !< \f$S^{-1/2}\f$
c
      integer l_a, k_a   !< Copy of the matrix \f$S\f$ or the
                         !< eigenvectors
      integer l_w, k_w   !< The eigenvalues
      integer l_wk, k_wk !< Workspace
      integer lwork      !< The size of the workspace
      integer info       !< Error flag
c
      lwork = -1 !< Ask for memory requirements
      info  = 0
      call dsyev('V','U',nbf,ov,nbf,ov12,ov12,lwork,info)
      lwork = ov12(1,1)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,'mat A',l_a,k_a))
     &   call errquit("wfn1_nxt2_ov12: could not allocate k_a",
     &                nbf*nbf,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,'vec W',l_w,k_w))
     &   call errquit("wfn1_nxt2_ov12: could not allocate k_w",
     &                nbf,MA_ERR)
      if (.not.ma_push_get(MT_DBL,lwork,'Work',l_wk,k_wk))
     &   call errquit("wfn1_nxt2_ov12: could not allocate k_wk",
     &                nbf,MA_ERR)
      call wfn1_nxt2_ov12_doit(nbf,lwork,ov,ov12,ovm12,dbl_mb(k_a),
     &                         dbl_mb(k_w),dbl_mb(k_wk))
      if (.not.ma_chop_stack(l_a))
     &  call errquit("wfn1_nxt2_ov12: could not dealloc memory",
     &               0,MA_ERR)
c
      end
C>
C> \brief Do the real work of wfn1_nxt2_ov12
C>
      subroutine wfn1_nxt2_ov12_doit(nbf,lwork,ov,ov12,ovm12,a,w,work)
      implicit none
c
#include "errquit.fh"
c
      integer nbf   !< [Input] The number of basis functions
      integer lwork !< [Input] The size of the work space
c
      double precision ov(nbf,nbf)   !< [Input] The overlap
                                     !< matrix \f$S\f$
      double precision ov12(nbf,nbf) !< [Output] The matrix
                                     !< \f$S^{1/2}\f$
      double precision ovm12(nbf,nbf) !< [Output] The matrix
                                      !< \f$S^{-1/2}\f$
      double precision a(nbf,nbf)    !< [Scratch]
      double precision w(nbf)        !< [Scratch] Eigenvalues
      double precision work(lwork)   !< [Scratch]
c
      integer i, j, k                !< Counters
      integer info                   !< Error flag
c
c     Create a copy of the matrix (will be destroyed)
c
      do j = 1, nbf
        do i = 1, nbf
          a(i,j) = ov(i,j)
        enddo
      enddo
c
c     Diagonalize the matrix (overwriting the input matrix with the
c     eigenvectors)
c
      call dsyev('V','U',nbf,a,nbf,w,work,lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_ov12_doit: dsyev failed",info,UERR)
      endif
c
c     Apply the function sqrt(lambda)
c
      do i = 1, nbf
        w(i) = sqrt(w(i))
      enddo
c
c     Construct the matrix \f$S^{1/2}\f$ in the original basis
c
      do j = 1, nbf
        do i = 1, nbf
          ov12(i,j) = 0.0d0
        enddo
      enddo
      do j = 1, nbf
        do i = 1, nbf
          do k = 1, nbf
            ov12(i,j) = ov12(i,j) + a(i,k)*w(k)*a(j,k)
          enddo
        enddo
      enddo
c
c     Apply the function 1/sqrt(lambda)
c
      do i = 1, nbf
        w(i) = 1.0d0/w(i)
      enddo
c
c     Construct the matrix \f$S^{-1/2}\f$ in the original basis
c
      do j = 1, nbf
        do i = 1, nbf
          ovm12(i,j) = 0.0d0
        enddo
      enddo
      do j = 1, nbf
        do i = 1, nbf
          do k = 1, nbf
            ovm12(i,j) = ovm12(i,j) + a(i,k)*w(k)*a(j,k)
          enddo
        enddo
      enddo
c
c     Done
c
      end
C> 
C> @}
