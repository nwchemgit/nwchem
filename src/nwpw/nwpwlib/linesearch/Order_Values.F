c     $Id$                       
*     *****************************
*     *                           *
*     *       Order_Values        *
*     *                           *
*     *****************************
*
*   this subroutine makes f(indx(1)) < f(indx(2)) < f(indx(3)) < ....
*   via a bubble sort
*   
*   Entry - n,f
*   Exit - returns indx
*

      subroutine Order_Values(n,f,indx)
      implicit none
      integer n
      real*8 f(*)
      integer indx(*)

*     ****** local variables *****
      integer i,j,idum

      do i=1,n
         indx(i) = i
      end do
      do i=1,n-1
        do j=i+1,n
           if (f(indx(j)).lt.f(indx(i))) then
              idum    = indx(i)
              indx(i) = indx(j)
              indx(j) = idum
           end if
         end do
       end do
            
      return 
      end

