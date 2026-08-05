C>
C> \brief On the fly matrix-vector multiplication
C>
C> This routine multiplies the COSMO matrix with the current guess
C> for the charge vector using a dot-product based algorithm. The matrix
C> is generated on the fly. For performance reasons the routine is
C> replicated data parallel.
C>
      subroutine hnd_cosaxd(nat,x,ax,nefc,efcc,efcs,efczz,efciat,ratm)
      implicit none
#include "global.fh"
#include "msgids.fh"
#include "cosmoP.fh"
#include "cosmo_params.fh"
c
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      double precision efczz(nefc)
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision x(nefc)      !< [Input] Vector `x`
      double precision ax(nefc)     !< [Output] Matrix-vector product
                                    !< `Ax`
c
      double precision zetai, zetaj, zetaij, pi, aij, dij, dum, xj
      double precision zero, one
      parameter (zero=0.0d+00, one=1.0d+00)
c
c
c
c     Introduced a trivial replicated data parallelization of this
c     matrix-vector multiplication
c
      double precision r, d, factor, factor2
      integer i, j
      r (i,j)=sqrt((efcc(1,i)-efcc(1,j))**2+
     1             (efcc(2,i)-efcc(2,j))**2+
     2             (efcc(3,i)-efcc(3,j))**2)
      d (i  )=adiag/sqrt(efcs(i))
c
      pi=acos(-1.0d0)
      call dfill(nefc,0.0d0,ax,1)
c
      if (do_cosmo_model.eq.DO_COSMO_KS) then
        do i=ga_nodeid()+1,nefc,ga_nnodes()
          ax(i) = ax(i) + d(i)*x(i)
          do j=i+1,nefc
            aij = one/r(i,j)
            ax(i) = ax(i) + x(j)*aij
            ax(j) = ax(j) + x(i)*aij
          enddo
        enddo
      else if (do_cosmo_model.eq.DO_COSMO_YK) then
        factor = zeta*dsqrt(ptspatm/(4d0*pi))
        factor2 = 2d0/sqrt(pi)
        do i=ga_nodeid()+1,nefc,ga_nnodes()
          zetai = efczz(i)
          ax(i) = ax(i) + efcs(i)*x(i)
          do j=i+1,nefc
            zetaj = efczz(j)
            zetaij= zeta*zetaj/sqrt(zetai**2+zetaj**2)
            dij=r(i,j)
            if (dij.lt.1.0d-8) then
              aij=factor2*zetaij*(1d0 - (zetaij*dij)**2/3d0)
            else
              aij=erf(zetaij*dij)/dij
            endif
            ax(i)=ax(i)+aij*x(j)
            ax(j)=ax(j)+aij*x(i)
          enddo
        enddo
      endif
c
      call ga_mask_sync(.true.,.false.)
      call ga_dgop(msg_cosaxd,ax(1),nefc,'+')
c
      return
      end
