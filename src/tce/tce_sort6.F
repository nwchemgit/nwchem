      subroutine tce_sort_6(unsorted,sorted,a,b,c,d,e,f,
     1                      i,j,k,l,m,n,factor)
c
c $Id$
c
      implicit none
c#include "global.fh"
c#include "mafdecls.fh"
c#include "sym.fh"
#include "util.fh"
c#include "stdio.fh"
c#include "errquit.fh"
#include "tce_main.fh"
      integer a,b,c,d,e,f
      integer i,j,k,l,m,n
      integer id(6),jd(6),ia,ib,j1,j2,j3,j4,j5,j6
      double precision sorted(a*b*c*d*e*f)
      double precision unsorted(a*b*c*d*e*f)
      double precision factor
      double precision cpu
c
c     Plain original algorithm
c
      cpusecs(7) = cpusecs(7) - util_cpusec()
      cpusecs(57) = cpusecs(57) - util_wallsec()
      cpu = -util_cpusec()
      ia = 0
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      if ((n.eq.6).or.(n.eq.5)) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j3 = 1,c
         id(3) = j3
         do j4 = 1,d
          id(4) = j4
          do j5 = 1,e
           id(5) = j5
           do j6 = 1,f
            id(6) = j6
            ia = id(6)+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     1         *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))
            ib = id(n)+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     1         *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
            sorted(ib) = unsorted(ia) * factor
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.4) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j3 = 1,c
         id(3) = j3
         do j5 = 1,e
          id(5) = j5
          do j4 = 1,d
           id(4) = j4
           do j6 = 1,f
            id(6) = j6
            ia = id(6)+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     1         *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))
            ib = id(n)+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     1         *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
            sorted(ib) = unsorted(ia) * factor
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.3) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j4 = 1,d
         id(4) = j4
         do j5 = 1,e
          id(5) = j5
          do j3 = 1,c
           id(3) = j3
           do j6 = 1,f
            id(6) = j6
            ia = id(6)+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     1         *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))
            ib = id(n)+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     1         *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
            sorted(ib) = unsorted(ia) * factor
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.2) then
      do j1 = 1,a
       id(1) = j1
       do j3 = 1,c
        id(3) = j3
        do j4 = 1,d
         id(4) = j4
         do j5 = 1,e
          id(5) = j5
          do j2 = 1,b
           id(2) = j2
           do j6 = 1,f
            id(6) = j6
            ia = id(6)+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     1         *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))
            ib = id(n)+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     1         *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
            sorted(ib) = unsorted(ia) * factor
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else
      do j2 = 1,b
       id(2) = j2
       do j3 = 1,c
        id(3) = j3
        do j4 = 1,d
         id(4) = j4
         do j5 = 1,e
          id(5) = j5
          do j1 = 1,a
           id(1) = j1
           do j6 = 1,f
            id(6) = j6
            ia = id(6)+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     1         *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))
            ib = id(n)+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     1         *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
            sorted(ib) = unsorted(ia) * factor
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      cpusecs(7) = cpusecs(7) + util_cpusec()
      cpusecs(57) = cpusecs(57) + util_wallsec()
      cpu = cpu + util_cpusec()
c     write(*,*) cpu
      return
      end
