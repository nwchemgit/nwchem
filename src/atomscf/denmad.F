      subroutine denmad(d,c,ns,nb,occ,nrow)
C$Id$
      implicit none
c.......................................................................
c     make actual density matrices.
c.......................................................................
      integer ns, nb, nrow
      real *8  d(*), c(nb,nb), occ
c
      integer i, j, k
      integer icount, klast
      real *8 sum
c
      icount = 0
      do i = 1 , nb
         do j = 1 , i
            icount = icount + 1
            klast = nrow + ns - 1
            sum = 0.0d0
            do k = nrow , klast
               sum = sum + c(i,k)*c(j,k)
            enddo
            sum = sum*occ
            if (i.ne.j) sum = sum + sum
            d(icount) = d(icount) + sum
         enddo
      enddo
      return
      end
