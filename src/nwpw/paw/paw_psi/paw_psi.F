


*     ***********************************
*     *                                 *
*     *          paw_psi_MakeOrtho      *
*     *                                 *
*     ***********************************

*   This routine orthonormalizes the orbitals using a modified
* Gram-Schmidt algorithm.
*
      subroutine paw_psi_MakeOrtho(npack,ne,psi)
      implicit none
      integer npack,ne
      double complex psi(npack,ne)

*     **** local variables ****
      integer j,k
      real*8  w

      do k=1,ne
         call paw_overlap_matrix_gen(1,psi(1,k),psi(1,k),w)
         w = 1.0d0/dsqrt(w)
         call Pack_c_SMul(1,w,psi(1,k),psi(1,k))

         do j=k+1,ne
            call paw_overlap_matrix_gen(1,psi(1,k),psi(1,j),w)
            w = -w
            call Pack_cc_daxpy(1,w,psi(1,k),psi(1,j))
         end do
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          paw_psi_CheckOrtho     *
*     *                                 *
*     ***********************************

*   This routine return true if the orbitals are
* orthonormal.

      function paw_psi_CheckOrtho(npack,ne,psi)
      implicit none
      integer npack,ne
      double complex psi(npack,ne)

      real*8 paw_psi_CheckOrtho  !** RESULT **

*     **** local variables ****
      integer j,k
      real*8  w

      error = 0.0d0
      do k=1,ne
        call paw_overlap_matrix_gen(1,psi(1,k),psi(1,k),w)
        error = error + dabs(1.0d0-w)

        do j=k+1,ne
           call paw_overlap_matrix_gen(1,psi(1,j),psi(1,k),w)
           error = error + dabs(w)
        end do
      end do

      paw_psi_CheckOrtho = error
      return
      end
