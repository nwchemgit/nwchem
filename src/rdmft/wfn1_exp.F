c-----------------------------------------------------------------------
c
      subroutine wfn1_exp(alpha,g_x,g_t,g_exp)
      implicit none
C>
C> \brief Calculate exp(alpha*X) where X is a matrix
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
      do while (dnormXn.gt.1.0d-20)
        mxpow = mxpow + 1
        dnormXn = dnormXn*alpha/dble(mxpow)
      enddo
c
      call ga_copy(g_x,g_exp)
      call ga_scale(g_exp,alpha/dble(mxpow))
      call ga_distribution(g_exp,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_exp,ii,ii,ii,ii,1.0d0,0,1.0d0)
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
          call ga_acc(g_exp,ii,ii,ii,ii,1.0d0,0,1.0d0)
        enddo
      enddo
      end
c
c-----------------------------------------------------------------------
