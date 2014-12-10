      subroutine tce_sortacc_0(unsorted,sorted,factor)
c
c $Id$
c
      implicit none
      double precision sorted
      double precision unsorted
      double precision factor
c
c     Plain original algorithm
c
      sorted = sorted + unsorted * factor
      return
      end
      subroutine tce_sortacc_2(unsorted,sorted,a,b,
     1                         i,j,factor)
c
c $Id$
c
      implicit none
#include "util.fh"
#include "tce_main.fh"
      integer a,b
      integer i,j
      integer id(2),jd(2),ia,ib,j1,j2
      double precision sorted(a*b)
      double precision unsorted(a*b)
      double precision factor
c
c     Plain original algorithm
c
      cpusecs(7) = cpusecs(7) - util_cpusec()
      cpusecs(57) = cpusecs(57) - util_wallsec()
      ia = 0
      jd(1) = a
      jd(2) = b
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        ia = ia + 1
        ib = id(j)+jd(j)*(id(i)-1)
        sorted(ib) = sorted(ib) + unsorted(ia) * factor
       enddo
      enddo
      cpusecs(7) = cpusecs(7) + util_cpusec()
      cpusecs(57) = cpusecs(57) + util_wallsec()
      return
      end
      subroutine tce_sortacc_4(unsorted,sorted,a,b,c,d,
     1                         i,j,k,l,factor)
c
c $Id$
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
      subroutine tce_sortacc_6(unsorted,sorted,a,b,c,d,e,f,
     1                         i,j,k,l,m,n,factor)
c
c $Id$
c
      implicit none
#include "util.fh"
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
            sorted(ib) = sorted(ib) + unsorted(ia) * factor
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
            sorted(ib) = sorted(ib) + unsorted(ia) * factor
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
            sorted(ib) = sorted(ib) + unsorted(ia) * factor
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
            sorted(ib) = sorted(ib) + unsorted(ia) * factor
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
            sorted(ib) = sorted(ib) + unsorted(ia) * factor
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
      subroutine tce_sortacc_8(unsorted,sorted,a,b,c,d,e,f,g,h,
     1                         i,j,k,l,m,n,o,p,factor)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce_main.fh"
      integer a,b,c,d,e,f,g,h
      integer i,j,k,l,m,n,o,p
      integer id(8),jd(8),ia,ib,j1,j2,j3,j4,j5,j6,j7,j8
      double precision sorted(a*b*c*d*e*f*g*h)
      double precision unsorted(a*b*c*d*e*f*g*h)
      double precision factor
c
c     Plain original algorithm
c
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      jd(7) = g
      jd(8) = h
      if ((p.eq.8).or.(p.eq.7)) then
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
            do j7 = 1,g
             id(7) = j7
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (p.eq.6) then
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
           do j7 = 1,g
            id(7) = j7
            do j6 = 1,f
             id(6) = j6
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (p.eq.5) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j3 = 1,c
         id(3) = j3
         do j4 = 1,d
          id(4) = j4
          do j6 = 1,f
           id(6) = j6
           do j7 = 1,g
            id(7) = j7
            do j5 = 1,e
             id(5) = j5
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (p.eq.4) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j3 = 1,c
         id(3) = j3
         do j5 = 1,e
          id(5) = j5
          do j6 = 1,f
           id(6) = j6
           do j7 = 1,g
            id(7) = j7
            do j4 = 1,d
             id(4) = j4
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (p.eq.3) then
      do j1 = 1,a
       id(1) = j1
       do j2 = 1,b
        id(2) = j2
        do j4 = 1,d
         id(4) = j4
         do j5 = 1,e
          id(5) = j5
          do j6 = 1,f
           id(6) = j6
           do j7 = 1,g
            id(7) = j7
            do j3 = 1,c
             id(3) = j3
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (p.eq.2) then
      do j1 = 1,a
       id(1) = j1
       do j3 = 1,c
        id(3) = j3
        do j4 = 1,d
         id(4) = j4
         do j5 = 1,e
          id(5) = j5
          do j6 = 1,f
           id(6) = j6
           do j7 = 1,g
            id(7) = j7
            do j2 = 1,b
             id(2) = j2
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
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
          do j6 = 1,f
           id(6) = j6
           do j7 = 1,g
            id(7) = j7
            do j1 = 1,a
             id(1) = j1
             do j8 = 1,h
              id(8) = j8
              ia = id(8)+jd(8)*(id(7)-1+jd(7)
     1           *(id(6)-1+jd(6)*(id(5)-1+jd(5)*(id(4)-1+jd(4)
     2           *(id(3)-1+jd(3)*(id(2)-1+jd(2)*(id(1)-1)))))))
              ib = id(p)+jd(p)*(id(o)-1+jd(o)
     1           *(id(n)-1+jd(n)*(id(m)-1+jd(m)*(id(l)-1+jd(l)
     2           *(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))))
              sorted(ib) = sorted(ib) + unsorted(ia) * factor
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      return
      end
