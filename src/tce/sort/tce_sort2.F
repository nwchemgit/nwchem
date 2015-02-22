      subroutine tce_sort_2(unsorted,sorted,a,b,i,j,factor)
c
c $Id$
c
      implicit none
c#include "util.fh"
c#include "tce_main.fh"
      integer a,b
      integer i,j
      integer id(2),jd(2),ia,ib,j1,j2
      double precision sorted(a*b)
      double precision unsorted(a*b)
      double precision factor
c      cpusecs(7) = cpusecs(7) - util_cpusec()
c      cpusecs(57) = cpusecs(57) - util_wallsec()
      ia = 0
      jd(1) = a
      jd(2) = b
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        ia = ia + 1
        ib = id(j)+jd(j)*(id(i)-1)
        sorted(ib) = unsorted(ia) * factor
       enddo
      enddo
c      cpusecs(7) = cpusecs(7) + util_cpusec()
c      cpusecs(57) = cpusecs(57) + util_wallsec()
      return
      end
