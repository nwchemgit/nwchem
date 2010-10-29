      subroutine sym_get_cart_op(geom, opnum, matrix)
C$Id$
      implicit none
c
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom, opnum
      double precision matrix(3,4)
      double precision tmp(3,4)
      integer i, j, k
      integer opoff
c
c     Return the matrix representation of the operator that
c     acts on CARTESIAN coordinates.
c
c     First 3 columns correspond to the point group operator
c     and the final column is the translation
c
c     OP * r = r'(*) = matrix(*,1:3)*r(1:3) + matrix(*,4)
c
c     To account for solid state fractional coords
c
c     1) For non-translational bit
c
c     R(cartesian) = A * R(fractional) * Ainv
c
c     2) For translation bit
c
c     t(cartesian) = A * t(fractional)
c
c
**      This gets the operator for fractional coords
**      do j = 1, 4
**         do i = 1, 3
**           matrix(i,j) = sym_ops((opnum-1)*3+i,j,geom)
**        enddo
**      enddo
c
c     Molecules are stored in cartesians
c
      opoff = (opnum-1)*3
      if (isystype(geom) .eq. 0) then
        do j = 1, 4
           do i = 1, 3
              matrix(i,j) = sym_ops(opoff+i,j,geom)
           enddo
        enddo
        return
      endif
c
c     Get the matrix from the array premultiplying by A
c
      do j = 1, 4
         do i = 1, 3
            tmp(i,j) = 0.0d0
            do k = 1, 3
               tmp(i,j) = tmp(i,j) +
     $              amatrix(i,k,geom)*
     $              sym_ops(opoff+k,j,geom)
            enddo
         enddo
      enddo
c
c     Post multiply the non-translation bit by Ainv
c
      do j = 1, 4
         do i = 1, 3
            matrix(i,j) = 0.0d0
            do k = 1, 3
               matrix(i,j) = matrix(i,j) +
     $              tmp(i,k) *
     $              amatrix_inv(k,j,geom)
            enddo
         enddo
      enddo
c
c     Copy over the translation piece
c
      do i = 1, 3
         matrix(i,4) = tmp(i,4)
      enddo
c
      end
