      subroutine tce_sortacc_4(unsorted,sorted,a,b,c,d,
     1                         i,j,k,l,factor)
c
c $Id: tce_sortacc.F 26500 2014-12-10 05:05:57Z jhammond $
c
      implicit none
#include "util.fh"
#include "tce_main.fh"
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision sorted(a*b*c*d)
      double precision unsorted(a*b*c*d)
      double precision factor
c
c     Plain original algorithm
c
      cpusecs(7) = cpusecs(7) - util_cpusec()
      cpusecs(57) = cpusecs(57) - util_wallsec()
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      if ((l.eq.4).or.(l.eq.3)) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j3 = 1,c
         id(3) = j3
         do j4 = 1,d
          id(4) = j4
          ia = id(4)+jd(4)*(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))
          ib = id(l)+jd(l)*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          sorted(ib) = sorted(ib) + unsorted(ia) * factor
         enddo
        enddo
       enddo
      enddo
      else if (l.eq.2) then
      do j1 = 1,a
       id(1) = j1
       do j3 = 1,c
        id(3) = j3
        do j2 = 1,b
         id(2) = j2
         do j4 = 1,d
          id(4) = j4
          ia = id(4)+jd(4)*(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))
          ib = id(l)+jd(l)*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          sorted(ib) = sorted(ib) + unsorted(ia) * factor
         enddo
        enddo
       enddo
      enddo
      else
      do j2 = 1,b
       id(2) = j2
       do j3 = 1,c
        id(3) = j3
        do j1 = 1,a
         id(1) = j1
         do j4 = 1,d
          id(4) = j4
          ia = id(4)+jd(4)*(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))
          ib = id(l)+jd(l)*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          sorted(ib) = sorted(ib) + unsorted(ia) * factor
         enddo
        enddo
       enddo
      enddo
      endif
      cpusecs(7) = cpusecs(7) + util_cpusec()
      cpusecs(57) = cpusecs(57) + util_wallsec()
      return
      end
