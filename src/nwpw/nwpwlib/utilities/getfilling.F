      subroutine getfilling(gammapoint,n,filling)
      implicit none
      logical gammapoint
      integer n
      integer filling(4,*)

*     ***** local variables ****
      integer h,i,j,k
      integer f,f1,f2,f3,f4

      f=1
      filling(1,f) = 0
      filling(2,f) = 0
      filling(3,f) = 0
      filling(4,f) = 0

      if (gammapoint) then
         do k=1,5000
         do j=0,k
         do i=0,j
         do h=1,2
            f = f+1
            if (f.gt.n) go to 100

            f1 =  i
            if (i.eq.0) then
              f2 = j
              f3 = k
            else
               if (mod(j,2).eq.0) then
                  f2 = j/2
               else
                  f2 = -j/2
               end if
               if (mod(k,2).eq.0) then
                  f3 = k/2
               else
                  f3 = -k/2
               end if
            end if

            if (h.eq.1) f4 =  2
            if (h.eq.2) f4 = -2

            filling(1,f) = f1
            filling(2,f) = f2
            filling(3,f) = f3
            filling(4,f) = f4
         end do
         end do
         end do
         end do
  100    continue

      else
         do k=2,5000
         do j=1,k
         do i=1,j
         do h=1,2
            f = f+1
            if (f.gt.n) go to 200

            if (mod(i,2).eq.0) then
               f1 =  i/2
            else
               f1 = -i/2
            end if
            if (mod(j,2).eq.0) then
               f2 = j/2
            else
               f2 = -j/2
            end if
            if (mod(k,2).eq.0) then
               f3 = k/2
            else
               f3 = -k/2
            end if
            if (h.eq.1) f4 =  2
            if (h.eq.2) f4 = -2

            filling(1,f) = f1
            filling(2,f) = f2
            filling(3,f) = f3
            filling(4,f) = f4
         end do
         end do
         end do
         end do
  200    continue
      end if

      return
      end
c $Id$
