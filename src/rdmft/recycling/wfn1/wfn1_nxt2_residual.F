C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Construct the residual matrix given a Fock matrix
C>
C> Assume we are given a Fock matrix in a particular basis. Then for
C> every basis function we have that the residual is
C> \f{eqnarray*}{
C>    r_i &=& F\phi_i - \langle E\rangle \phi_i \\\\
C>    \langle E\rangle &=& F_{ii}
C> \f}
C> Based on this we can construct a residual matrix
C> \f{eqnarray*}{
C>   R_{ij} &=& (r_i,r_j) \\\\
C>   R_{ij} &=& \langle\phi_i|FF|\phi_j\rangle
C>   - (\langle E_i\rangle + \langle E_j\rangle)
C>     \langle\phi_i|F|\phi_j\rangle
C>   + \langle E_i\rangle\langle E_j\rangle\langle\phi_i|\phi_j\rangle
C> \f}
C> In this routine we assume that we are working in an orthonormal
C> basis.
C>
C> The hope is that we can find a nearby stationary point by
C> diagonalizing this matrix.
C>
      subroutine wfn1_nxt2_residual(nbf,fock,res)
      implicit none
      integer nbf !< [Input] The number of basis functions
c
      double precision fock(nbf,nbf) !< [Input] The Fock matrix
      double precision res(nbf,nbf)  !< [Output] The residual matrix
c
      integer i,j !< Counters
c
      double precision ei, ej, ddot
c
      do j = 1, nbf
        ej = fock(j,j)
        do i = 1, nbf
          res(i,j) = ddot(nbf,fock(1,i),1,fock(1,j),1)
          ei = fock(i,i)
          res(i,j) = res(i,j) - (ei+ej)*fock(i,j)
        enddo
        res(j,j) = res(j,j) + ej*ej
      enddo
c
      end
C>
C> \brief Construct the full residual vector
C>
C> Construct the residual matrix for every Fock matrix in the Fock
C> vector.
C>
      subroutine wfn1_nxt2_fr(nea,neb,nbf,fock,res)
      implicit none
c
#include "wfn1_nxt.fh"
c
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
      integer nbf !< [Input] The number of basis functions
c
      double precision fock(8*nbf*nbf) !< [Input] The Fock vector
      double precision res(8*nbf*nbf)  !< [Output] The residual vector
c
      integer i, j, ioff
c
      call wfn1_nxt2_residual(nbf,fock(wfn1_nxt_aocc()),
     +                            res(wfn1_nxt_aocc()))
c     ioff = wfn1_nxt_aocc()
c     do j = nea+1, nbf
c       do i = nea+1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = 0.0d0
c         endif
c       enddo
c     enddo
c     do j = 1, nbf
c       do i = 1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = res(ioff+(i-1)+(j-1)*nbf)/1.0d2
c         endif
c       enddo
c     enddo
c
      call wfn1_nxt2_residual(nbf,fock(wfn1_nxt_bocc()),
     +                            res(wfn1_nxt_bocc()))
c     ioff = wfn1_nxt_bocc()
c     do j = neb+1, nbf
c       do i = neb+1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = 0.0d0
c         endif
c       enddo
c     enddo
c     do j = 1, nbf
c       do i = 1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = res(ioff+(i-1)+(j-1)*nbf)/1.0d2
c         endif
c       enddo
c     enddo
c
      call wfn1_nxt2_residual(nbf,fock(wfn1_nxt_anat()),
     +                            res(wfn1_nxt_anat()))
c     ioff = wfn1_nxt_anat()
c     do j = 1, nbf
c       do i = 1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = res(ioff+(i-1)+(j-1)*nbf)/1.0d2
c         endif
c       enddo
c     enddo
c
      call wfn1_nxt2_residual(nbf,fock(wfn1_nxt_bnat()),
     +                            res(wfn1_nxt_bnat()))
c     ioff = wfn1_nxt_bnat()
c     do j = 1, nbf
c       do i = 1, nbf
c         if (i.ne.j) then
c           res(ioff+(i-1)+(j-1)*nbf) = res(ioff+(i-1)+(j-1)*nbf)/1.0d2
c         endif
c       enddo
c     enddo
c
      end
C>
C> @}

