      subroutine orderd(amat,nrow,ncol,imemb,nblock,ind,vec,iflag)
C$Id: orderd.F,v 1.2 1995-02-02 23:09:04 d3g681 Exp $
      implicit none
c.......................................................................
c     this routine sorts a set of column vectors in amat(*,*) according
c     to increasing values in vec(*). care is taken that the output orde
c     of the vectors is as close to the input order as possible.
c.......................................................................
      integer nrow, ncol, nblock, iflag
      integer ind(*), imemb(*)
      real *8  amat(nrow,ncol), vec(*)
c.......................................................................
c     set tolerance for test, then determine for each element how
c     many smaller elements vec(*) contains.
c.......................................................................
      real *8 tol, test, scra
      integer i, j
      integer icount, indi, itest, im1, ilow, index
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
c     if desired,scan ind(*) to determine the number of different values
c     in vec(*) and the number of elements of each value.
c.......................................................................
      if (iflag.le.99) then
         do i = 1 , ncol
            imemb(i) = 0
         enddo
         do i = 1 , ncol
            indi = ind(i)
            imemb(indi) = imemb(indi) + 1
         enddo
         icount = 0
         do i = 1 , ncol
            if (imemb(i).ne.0) then
               icount = icount + 1
               imemb(icount) = imemb(i)
            end if
         enddo
         nblock = icount
      end if
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
 110  index = ind(i)
      do j = 1 , nrow
         scra = amat(j,i)
         amat(j,i) = amat(j,index)
         amat(j,index) = scra
      enddo
      scra = vec(i)
      vec(i) = vec(index)
      vec(index) = scra
      ind(i) = ind(index)
      ind(index) = index
      if (ind(i).ne.i) go to 110
c.......................................................................
c
c     amat(*,*) is ordered through i. go back to checking order.
c.......................................................................
      ilow = i + 1
      if (ilow.lt.ncol) go to 90
 130  return
      end
