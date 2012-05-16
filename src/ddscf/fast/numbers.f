*
* $Id$
*
      program numbers
      implicit integer (a-z)
c
      maxdiff = 0
      do i = 7,2048
         n = i
         call fastj_optimize_n(n,'upperbound')
         if (abs(n-i) .gt. maxdiff) then
            maxdiff = abs(n-i)
            write(6,*) i, n, maxdiff
         endif
      enddo
c
      write(6,*) ' upperbound ', maxdiff
c
      maxdiff = 0
      do i = 7,2048
         n = i
         call fastj_optimize_n(n,'nearby')
         if (abs(n-i) .gt. maxdiff) then
            maxdiff = abs(n-i)
            write(6,*) i, n, maxdiff
         endif
      enddo
      write(6,*) ' nearby ', maxdiff
c
      end
      subroutine fastj_optimize_n(n, mode)
      implicit none
      integer n
      character*(*) mode
c
      integer n2, n3, n5, ii2, ii3, ii5, i2, i3, i5
      double precision test, diff
c
      n = n + 1                 ! Temporarily to force to power of radix
c
      n2 = log(dble(n))/log(2.0d0) + 1
      n3 = log(dble(n))/log(3.0d0) + 1
      n5 = log(dble(n))/log(5.0d0) + 1
      ii2 = 0
      ii3 = 0
      ii5 = 0
      diff = n
      do i2 = 0, n2
         do i3 = 0, n3
            do i5 = 0, n5
               test = (2.0d0**i2)*(3.0d0**i3)*(5.0d0**i5) + 1d-6
               if (test.ge.n .or. mode.ne.'upperbound') then
                  if (abs(test-n) .lt. diff) then
                     diff = abs(test - n)
                     ii2 = i2
                     ii3 = i3
                     ii5 = i5
                  endif
               endif
            enddo
         enddo
      enddo
c
      n = (2.0d0**ii2)*(3.0d0**ii3)*(5.0d0**ii5)  - 1
c
      end
      
