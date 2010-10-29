      subroutine orderd(amat,nrow,ncol,ind,vec)
C$Id$
      implicit none
c.......................................................................
c     this routine sorts a set of column vectors in amat(*,*) according
c     to increasing values in vec(*). care is taken that the output orde
c     of the vectors is as close to the input order as possible.
c.......................................................................
      integer nrow, ncol
      integer ind(*)
      real *8  amat(nrow,ncol), vec(*)
c.......................................................................
c     set tolerance for test, then determine for each element how
c     many smaller elements vec(*) contains.
c.......................................................................
      real *8 tol, test, scra
      integer i, j
      integer indi, itest, im1, ilow, iindx
c
      tol = 1.0d-10
      do i = 1 , ncol
         test = vec(i) - tol
         indi = 1
         do j = 1 , ncol
            if (vec(j).lt.test) indi = indi + 1
         enddo
         ind(i) = indi
      enddo
c.......................................................................
c
c     establish order in degeneracies of the ordering vector.
c.......................................................................
      do i = 2 , ncol
         itest = ind(i)
         im1 = i - 1
         do j = 1 , im1
            if (ind(j).eq.itest) itest = itest + 1
         enddo
         ind(i) = itest
      enddo
c.......................................................................
c
c     ind(*) contains ordering indices for amat(*,*). sort following
c     input order as far as this is correct.
c.......................................................................
      ilow = 1
 90   do i = 1 , ncol
         if (i.ne.ind(i)) go to 110
      enddo
      go to 130
c.......................................................................
c
c     input order wrong. swap present vector into correct col.
c.......................................................................
 110  iindx = ind(i)
      do j = 1 , nrow
         scra = amat(j,i)
         amat(j,i) = amat(j,iindx)
         amat(j,iindx) = scra
      enddo
      scra = vec(i)
      vec(i) = vec(iindx)
      vec(iindx) = scra
      ind(i) = ind(iindx)
      ind(iindx) = iindx
      if (ind(i).ne.i) go to 110
c.......................................................................
c
c     amat(*,*) is ordered through i. go back to checking order.
c.......................................................................
      ilow = i + 1
      if (ilow.lt.ncol) go to 90
 130  return
      end
