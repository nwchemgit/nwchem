c-----------------------------------------------------------------------
c
      subroutine wfn1_exp(alpha,g_x,g_t,g_exp)
      implicit none
C>
C> \brief Calculate exp(alpha*X) where X is a skew-symmetric matrix
C>
C> NOTE: This routine is serial and does not synchronize. The rational
C>       for the latter is that in the application it is always followed
C>       by a ga_dgemm (which does synchronize). In addition where we
C>       invoke this is invoked 4 times where we use it, so where a
C>       single synchronization is sufficient to guarantee completion of
C>       all 4 matrices.
C> 
C> Calculate the exponent of a skew-symmetric matrix by the usual
C> way of calculating matrix functions. I.e. diagonalize the matrix,
C> apply the function to the eigenvalues, and transform back. In this
C> case that is a sensible approach because the eigenvalues of a
C> skew-symmetric matrix are purely imaginary. As a result the exponents
C> of these eigenvalues are also imaginary numbers with a norm of 1. 
C> Therefore the reconstruction of the final unitary matrix can be
C> completed without loss of accuracy. Hence this approach is not 
C> bitten by the issues that affect exponentiation of most matrices [1].
C>
C> The main downside of this approach is that it relies on matrix
C> diagonalization which is not the most scalable method in any case.
C> However, in this situation it is even worse in that the parallel
C> equivalent of the LAPACK ZGEEV routine is not available in 
C> ScaLAPACK. Hence we are forced for now to do the diagonalization
C> in serial. A parallel diagonalization routine for skew-Hermitian
C> matrices is available in the Elemental library which is written
C> in C++. So we cannot use that here but we should be able to use
C> that in NWChemEx. 
C>
C> Another issue with ZGEEV is that in some cases (e.g. Hartree-Fock)
C> the virtual-virtual block of the input matrix is 0. This seems to
C> throw ZGEEV and it produces garbage for some of the eigenvectors.
C> So in reality g_exp is not strictly exp(g_x) only part of it is.
C> This means that the result has to be used with some care, in
C> particular if it is used to transform the orbitals a subsequent
C> orthonormalization is recommended.
C>
C> Furthermore, the matrix exponentiation is part of the linesearch
C> algorithm. Hence we could do the diagonalization only once, keep
C> the result and evaluate exp(alpha*X) followed by a reconstruction
C> for all the points on the linesearch. This should reduce the 
C> dependence of the overall performance on the performance of the 
C> diagonalization routine by at least a factor 2 (every linesearch
C> takes at least 2 steps). 
C>
C> References
C>
C> [1] C. Moler, C. van Loan, "Nineteen dubious ways to compute the
C>     exponential of a matrix, twenty-five years later." 
C>     SIAM reviews (2003) 45, pp. 3-49. Doi: 10.1137/S00361445024180
C>
#include "errquit.fh"
#include "global.fh"
c
      double precision, intent(in) :: alpha
c
      integer, intent(in)    :: g_x
      integer, intent(inout) :: g_t
      integer, intent(inout) :: g_exp
c
      double precision :: dnormX  ! The norm of alpha*X
      double precision :: dnormXn ! The norm of the n-th term in the
                                  ! Taylor series of exp(alpha*X)
c
      double complex, allocatable :: skew(:,:) ! skew symmetric mat
c     double complex, allocatable :: skvl(:,:) ! left eigenvectors
      double complex              :: skvl      ! left eigenvectors
      double complex, allocatable :: skvr(:,:) ! right eigenvectors
      double complex, allocatable :: skev(:)   ! eigenvalues
      double complex, allocatable :: work(:)   ! scratch space
c
      double precision, allocatable :: rwork(:)   ! scratch space
c
      double precision, allocatable, target :: buf(:)
      double precision, pointer :: block(:,:)
c
      double precision, allocatable :: u(:,:) ! the unitary matrix
      double precision              :: dnrm
      double precision, external    :: ddot
c
      integer :: itypx, nrowx, ncolx
      integer :: info ! status flag returned by LAPACK
c
      integer :: iproc, nproc
      integer :: ip, ii, jj, kk ! counters
      integer :: ilo, ihi, jlo, jhi ! local block limits
      integer :: idim, jdim
c
      logical :: o_repeat
c
      character*9 pname
      parameter(pname="wfn1_exp:")
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
c
        call ga_inquire(g_x,itypx,nrowx,ncolx)
        if (nrowx.ne.ncolx) then
          call errquit(pname//" matrix g_x not square",0,UERR)
        endif
        call ga_distribution(g_exp,0,ilo,ihi,jlo,jhi)
        idim = ihi-ilo+1
        jdim = jhi-jlo+1
        call ga_distribution(g_x,0,ilo,ihi,jlo,jhi)
        idim = max(idim,ihi-ilo+1)
        jdim = max(idim,jhi-jlo+1)
c
        allocate(skew(nrowx,nrowx))
c       allocate(skvl(nrowx,nrowx))
        allocate(skvr(nrowx,nrowx))
        allocate(skev(nrowx)      )
        allocate(work(10*nrowx)   )
        allocate(rwork(10*nrowx)  )
        allocate(buf(idim*jdim))
        allocate(u(nrowx,nrowx))
c
        do ip = 0, nproc-1
          call ga_distribution(g_x,ip,ilo,ihi,jlo,jhi)
          block(ilo:ihi,jlo:jhi)=>buf
          call ga_get(g_x,ilo,ihi,jlo,jhi,block,ihi-ilo+1)
          do jj = jlo, jhi
            do ii = ilo, ihi
              skew(ii,jj) = dcmplx(block(ii,jj),0.0d0)
            enddo
          enddo
        enddo
c
        call zgeev('N','V',nrowx,skew,nrowx,skev,skvl,nrowx,
     &             skvr,nrowx,work,10*nrowx,rwork,info)
        if (info.ne.0) then
          write(*,*)'wfn1_exp: Error in ZGEEV: ',info
          call errquit("wfn1_exp: ZGEEV failed",info,UERR)
        endif
c
        do ii = 1, nrowx
          skev(ii) = exp(alpha*skev(ii))
        enddo
c
c       Construct initial guess of U
c
        u = 0.0d0
        do jj = 1, nrowx
          do ii = 1, nrowx
            do kk = 1, nrowx
c
c             Note that the phases of the left and right eigenvectors
c             are independent of eachother. Therefore using both 
c             left and right eigenvectors may not result in a 
c             unitary matrix. Hence only the right eigenvectors are
c             used here.
c
c             Note that the eigenvalues of skew-symmetric matrices
c             come in conjugate pairs. Therefore, when summing over
c             all eigenvectors, the imaginary parts cancel, and we
c             may just use the real part only.
c
              u(ii,jj) = u(ii,jj)
     &                 + real(skvr(ii,kk)*skev(kk)*conjg(skvr(jj,kk)))
            enddo
          enddo
        enddo
c
        do ip = 0, nproc-1
          call ga_distribution(g_exp,ip,ilo,ihi,jlo,jhi)
          block(ilo:ihi,jlo:jhi)=>buf
          do jj = jlo, jhi
            do ii = ilo, ihi
              block(ii,jj) = u(ii,jj) 
            enddo
          enddo
          call ga_put(g_exp,ilo,ihi,jlo,jhi,block,ihi-ilo+1)
        enddo
c
        deallocate(u)
        deallocate(rwork)
        deallocate(work)
        deallocate(skew)
c       deallocate(skvl)
        deallocate(skvr)
        deallocate(skev)
        deallocate(buf)
c
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_exp_taylor(alpha,g_x,g_t,g_exp)
      implicit none
C>
C> \brief Calculate exp(alpha*X) where X is a matrix
C>
C> NOTE: This implementation does not work well when the norm of
C>       the matrix is large. In that case some of the intermediate
C>       terms get exceedingly big causing numerical problems.
C>
C> Calculate the exponent of a matrix by expanding the exponent as a
C> Taylor series. The algoritm uses only matrix-matrix multiplications.
C>
C> The result is returned in the global array g_exp, and g_t is used
C> as a scratch array. However, to avoid copying data the handles in
C> g_t and g_exp will be repeatedly swapped around. I.e. on return
C> g_t and g_exp may refer to different arrays than on input.
C>
C> Note also that the series is evaluated from the last term towards
C> the first term. This is done for reasons of accuracy as the last
C> term is the smallest one, and numerically the most accurate way
C> to sum a set of numbers in finite precision is to sum them in the
C> order of increasing value. The reason is that this preserves the
C> least significant digits as much as possible.
C>
C> Note that despite the above considerations calculating the
C> exponent of a matrix is notoriously difficult [1]. Hence this 
C> method is effective only under special conditions, in particular
C> alpha*||g_x|| should be less than some bound to ensure that the
C> Taylor series has no "hump". Even then the result cannot be expected
C> to be accurate to machine precision, the result may only be expected
C> to be an non-defective matrix.
C>
C> References
C>
C> [1] C. Moler, C. van Loan, "Nineteen dubious ways to compute the
C>     exponential of a matrix, twenty-five years later." 
C>     SIAM reviews (2003) 45, pp. 3-49. Doi: 10.1137/S00361445024180
C>
#include "errquit.fh"
#include "global.fh"
c
      double precision, intent(in) :: alpha
c
      integer, intent(in)    :: g_x
      integer, intent(inout) :: g_t
      integer, intent(inout) :: g_exp
c
      double precision :: dnormX  ! The norm of alpha*X
      double precision :: dnormXn ! The norm of the n-th term in the
                                  ! Taylor series of exp(alpha*X)
c
      integer :: itypx, nrowx, ncolx
      integer :: mxpow ! the maximum power
      integer :: pow   ! the current power
c
      integer :: iproc, nproc
      integer :: ii ! counter
      integer :: ilo, ihi, jlo, jhi ! local block limits
c
      character*9 pname
      parameter(pname="wfn1_exp:")
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call ga_inquire(g_x,itypx,nrowx,ncolx)
      if (nrowx.ne.ncolx) then
        call errquit(pname//" matrix g_x not square",0,UERR)
      endif
      call ga_normf(g_x,dnormX)
      dnormX = alpha*dnormX
      dnormXn = sqrt(1.0d0*nrowx)
      mxpow = 0
      do while (abs(dnormXn).gt.1.0d-20)
        mxpow = mxpow + 1
        dnormXn = dnormXn*dnormX/dble(mxpow)
      enddo
c
      call ga_copy(g_x,g_exp)
      call ga_scale(g_exp,alpha/dble(mxpow))
      call ga_distribution(g_exp,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_exp,ii,ii,ii,ii,1.0d0,1,1.0d0)
      enddo
      do pow = mxpow-1,1,-1
c
c       Inter change the global array handles
c
        ii    = g_exp
        g_exp = g_t
        g_t   = ii
c
c       Generate the summed series from pow to mxpow (apart from a
c       factor)
c
        call ga_dgemm('n','n',nrowx,nrowx,nrowx,alpha/dble(pow),g_x,g_t,
     &                0.0d0,g_exp)
        call ga_distribution(g_exp,iproc,ilo,ihi,jlo,jhi)
        do ii = max(ilo,jlo),min(ihi,jhi)
          call ga_acc(g_exp,ii,ii,ii,ii,1.0d0,1,1.0d0)
        enddo
      enddo
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_tst_exp(alpha,g_x,g_t,g_exp)
      implicit none
C>
C> \brief A test driver for wfn1_exp
C>
C> This module uses wfn1_exp to generate a unitary matrix by taken the
C> exponent of a skew-symmetric matrix. This routine calls through to
C> wfn1_exp and checks whether the resulting matrix is actually unitary.
C>
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
c
      double precision, intent(in) :: alpha
c
      integer, intent(in)    :: g_x
      integer, intent(inout) :: g_t
      integer, intent(inout) :: g_exp
c
      integer :: itypt, nrowt, ncolt
      integer :: iproc, nproc
      integer :: ii ! counter
      integer :: ilo, ihi, jlo, jhi ! local block limits
c
      double precision :: dnormT
c
      character*13 pname
      parameter(pname="wfn1_tst_exp:")
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
c     Check whether X + X^T = 0
c
      call ga_transpose(g_x,g_t)
      call ga_add(1.0d0,g_x,1.0d0,g_t,g_t)
      call ga_norm_infinity(g_t,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_exp: X+X^T: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(g_x)
c       call ga_print(g_exp)
      endif
      call ga_sync
c
      call wfn1_exp(alpha,g_x,g_t,g_exp)
c
c     Check whether U^T U = I
c
      call ga_sync
      call ga_inquire(g_exp,itypt,nrowt,ncolt)
      call ga_dgemm('t','n',nrowt,nrowt,nrowt,1.0d0,g_exp,g_exp,
     &              0.0d0,g_t)
      call ga_distribution(g_t,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_t,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_t,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_exp: TN: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(g_t)
c       call ga_print(g_x)
c       call ga_print(g_exp)
      endif
c
c     Check whether U U^T = I
c
      call ga_dgemm('n','t',nrowt,nrowt,nrowt,1.0d0,g_exp,g_exp,
     &              0.0d0,g_t)
      call ga_distribution(g_t,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_t,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_t,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_exp: NT: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(g_t)
c       call ga_print(g_x)
c       call ga_print(g_exp)
      endif
      end
c
c-----------------------------------------------------------------------
