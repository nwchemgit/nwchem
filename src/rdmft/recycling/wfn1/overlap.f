C>
C> \brief Calculate the LdR decomposition of the overlap matrix
C>
C> This subroutine is based on the function `wfn1_overlap_bo`. 
C> Whereas `wfn1_overlap_bo` only calculates the overlap between
C> two determinants, this routine also calculates the transformation
C> matrices of the orbitals as well as the diagonal matrix. These
C> additional matrices are needed to calculate the transition 
C> density matrices between the two determinants involved. 
C>
C> In short, given a non-symmetric overlap matrix \f$S\f$ this
C> subroutine calculates the left and right transformation matrices
C> \f$L\f$ and \f$R\f$ and the diagonal \f$d\f$ such that
C> \f{eqnarray*}{
C>   d &=& L\cdot S \cdot R
C> \f}
C> In addition the sign of the overlap is also provided. See [1] for
C> as discussion of this approach.
C> 
C> This routine is essentially a copy of the `subroutine determ` in
C> GAMESS-UK, written by J. Verbeek and J.H. van Lenthe.
C>
C> See also [2].
C>
C> ### References ###
C>
C> [1] Jacob Verbeek, Joop H. van Lenthe,
C>     "On the evaluation of non-orthogonal matrix elements",
C>     J. Mol. Struct. (Theochem) <b>229</b>, 115-137 (1991), DOI:
C>     <a href="https://doi.org/10.1016/0166-1280(91)90141-6">
C>     10.1016/0166-1280(91)90141-6</a>.
C>
C> [2] Joop H. van Lenthe, Gabriel G. Balint-Kurti,
C>     "The valence-bond self-consistent field method (VB-SCF):
C>     Theory and test calculations", J. Chem. Phys. <b>78</b>,
C>     5699-5713 (1983), DOI:
C>     <a href="https://doi.org/10.1063/1.445451">
C>     10.1063/1.445451</a>.
C>
      subroutine wfn1_overlap_ldr2(ne,orbov,lmat,rmat,dvec,ipar,tmp,m2)
      implicit none
c
#include "errquit.fh"
c
      integer ne           !< [Input] The number of particles
      integer ipar         !< [Output] The parity
c
      double precision orbov(ne,ne) !< [Input] The overlaps of extended
                                    !< orbitals
      double precision lmat(ne,ne)  !< [Output] The left transformation
                                    !< matrix
      double precision rmat(ne,ne)  !< [Output] The right transformation
                                    !< matrix
      double precision dvec(ne)     !< [Output] The eigenvalues
c
      double precision tmp(ne,ne) !< [Scratch]
      double precision m2(ne,*)  !< [Scratch]
c
      double precision det    !< The value of the matrix determinant
      double precision pivo   !< The value of the pivot
      double precision cridep !< Criterion for singularity
      parameter (cridep = 1.0d-13)
c
      integer i, j, k, l   !< Counters
      integer irank        !< The rank of the matrix done so far
c
      call dcopy(ne*ne,orbov,1,tmp,1)
      irank = 0
      ipar  = 1
      det   = 1.0d0
      if (ne.eq.0) return
      call wfn1_piv(tmp,ne,ne,1,pivo,ipar)
      if (dabs(pivo).lt.cridep) then
         irank = 0
         det   = 0.0d0
         return
      end if
      irank = 1
      do 40 i = 1,ne-1
         det = det * tmp(i,i)
         dvec(i) = tmp(i,i)
         do 30 j = i+1,ne
c.....
c.....form i-th row of r
c.....
            tmp(i,j) = -tmp(i,j)/tmp(i,i)
c.....
c.....adapt s-matrix to r
c.....
            do 10 k = i+1,ne
               tmp(k,j) = tmp(k,j) + tmp(k,i)*tmp(i,j)
10          continue
c.....
c.....adapt r-matrix to itself
c.....
            do 20 l = 1,i-1
               tmp(l,j) = tmp(l,j) + tmp(l,i)*tmp(i,j)
20          continue
30       continue
         call wfn1_piv(tmp,ne,ne,i+1,pivo,ipar)
         if (dabs(pivo).lt.cridep) then
            det = 0.0d0
            return
         end if
         irank = irank + 1
40    continue
      det = det * tmp(ne,ne) * ipar
      dvec(ne) = tmp(ne,ne)
c
      call dfill(ne*ne,0.0d0,lmat,1)
      call dfill(ne*ne,0.0d0,rmat,1)
      do i = 1, ne
        lmat(i,i) = 1.0d0
        rmat(i,i) = 1.0d0
        do j = i+1, ne
          rmat(i,j) = tmp(i,j)
        enddo
        do j = 1, i-1
          lmat(i,j) = tmp(i,j)
        enddo
      enddo
c
      return
      end
C>
C> \brief Swap rows and columns to get maximum values on the diagonal
C> of the matrix
C>
C> This routine swap rows and columns of a matrix to place the 
C> maximum element of the remaining matrix block on the diagonal.
C>
C> This subroutine was obtained from GAMESS-UK, with permission by
C> Joop H. van Lenthe (Jan, 2014).
C>
      subroutine wfn1_piv(a,nr,nc,n,pivo,ipar)
      implicit none
c
      integer nr !< [Input] The number of rows
      integer nc !< [Input] The number of columns
      integer n  !< [Input] The rank of the first element of the 
                 !< remaining block
      double precision a(nr,nc) !< [In/Output] The matrix to operate on
      double precision pivo     !< [Output] The pivot found
      integer ipar              !< [Output] The parity found
c
c     Local
c
      integer irt !< The current best row number
      integer ict !< The current best column number
      integer i, j, k, l !< Counters
c
      double precision aa !< Temporary variable
c
      pivo = 0.0d0
      irt = n
      ict = n
c
c...  search pivot
c
      do j = n,nc
         do i = n,nr
            if (dabs(pivo).lt.dabs(a(i,j))) then
               pivo = a(i,j)
               irt = i
               ict = j
            end if
         enddo ! i
      enddo ! j
c
      if (irt.ne.n) then
c
c...     permute rows
c
         do k = 1,nc
            aa       = a(irt,k)
            a(irt,k) = a(n,k)
            a(n,k)   = aa
         enddo
         ipar    =-ipar
c
      end if
c
      if (ict.ne.n) then
c
c...     permute columns
c
         do l = 1,nr
            aa       = a(l,ict)
            a(l,ict) = a(l,n)
            a(l,n)   = aa
         enddo
         ipar    =-ipar
c
      end if
      return
      end
