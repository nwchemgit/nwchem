      subroutine starcd(c,ss,nbci,ncshi,noshi)
C$Id$
      implicit none
#include "errquit.fh"
      integer nbci, ncshi, noshi
      real *8 ss(*), c(nbci,nbci)
      real *8  a(100)
c.......................................................................
c
c     this routine provides a set of schmidt
c     orthogonalized start vectors.
c
c.......................................................................
      integer i, j, k, l, m
      integer nshti, ixx, klo, khi, indxlm, jm1
      real *8 sum
c
c     inline indexing function for triangular matrices.
c
      integer index
      index(i,j) = max0(i,j)*(max0(i,j)-1)/2 + min0(i,j)
c
c.......................................................................
c     uneducated guess of trial vectors.
c.......................................................................
c
      nshti = ncshi + noshi
      if (nshti.eq.0) return
      if (nshti.gt.nbci) then
        write(6,*)' more atomic shells than basis ',
     &      'functions in atomic startup'
        call errquit('atomscf:startcd: fatal error',0,
     &       UNKNOWN_ERR)
      endif
      ixx = nbci/nshti
      khi = 0
      do j = 1 , nshti
         klo = khi + 1
         khi = khi + ixx
         do k = klo , khi
            c(k,j) = 1.0d0
         enddo
      enddo
c.......................................................................
c
c     orthogonalize.
c
c.......................................................................
      do j = 1 , nshti
c.......................................................................
c
c     take scalar products with preceding vectors of same symm.
c     first do a(*) = ss(*)*c(*,j).
c
c.......................................................................
         do 50 l = 1 , nbci
            sum = 0.0d0
            do 40 m = 1 , nbci
               indxlm = index(l,m)
               sum = sum + ss(indxlm)*c(m,j)
 40         continue
            a(l) = sum
 50      continue
c.......................................................................
c
c     if first vector of symmetry, normalize directly.
c
c.......................................................................
         if (j.ne.1) then
c.......................................................................
c
c     c(*,k)*a(*)
c
c.......................................................................
            jm1 = j - 1
            do k = 1 , jm1
               sum = 0.0d0
               do l = 1 , nbci
                  sum = sum + c(l,k)*a(l)
               enddo
c.......................................................................
c
c     multiply c(*,k) by the scalar product and subtract from c(*,j)
c
c.......................................................................
               do l = 1 , nbci
                  c(l,j) = c(l,j) - c(l,k)*sum
               enddo
            enddo
c.......................................................................
c
c     c(*,j) is orthogonal to all preceding vectors. normalize.
c
c.......................................................................
            do l = 1 , nbci
               sum = 0.0d0
               do m = 1 , nbci
                  indxlm = index(l,m)
                  sum = sum + ss(indxlm)*c(m,j)
               enddo
               a(l) = sum
            enddo
         end if
         sum = 0.0d0
         do l = 1 , nbci
            sum = sum + c(l,j)*a(l)
         enddo
         sum = sqrt(sum)
         do l = 1 , nbci
            c(l,j) = c(l,j)/sum
         enddo
      enddo
      return
      end
