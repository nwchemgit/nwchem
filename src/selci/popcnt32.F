      integer function selci_popcnt32(i)
*
* $Id: popcnt32.F,v 1.12 2003-10-17 22:58:17 carlfahl Exp $
*
#if defined(LINUX) || defined(LINUX64)
#define INTEGER_1 integer*1
#else
#define INTEGER_1 byte
#endif
#if defined(CRAY_T3D) || defined(CRAY_T3E) 
      integer i
      selci_popcnt32 = popcnt(i)
#elif defined(EXT_INT)
#include "cbits.fh"
      INTEGER_1 i(8)
c
c     pop count of 64 bit word for sun ... the routine inibit must
c     be called once before popcnt64 or leadz64 is used ... this
c     computes ipopcb and leadzb
c
      ind1 = i(1)
      ind2 = i(2)
      ind3 = i(3)
      ind4 = i(4)
      ind5 = i(5)
      ind6 = i(6)
      ind7 = i(7)
      ind8 = i(8)
      selci_popcnt32 = ipopcb(ind1) + ipopcb(ind2) +
     $     ipopcb(ind3) + ipopcb(ind4) + ipopcb(ind5) +
     $     ipopcb(ind6) + ipopcb(ind7) + ipopcb(ind8) 
c
#else
#include "cbits.fh"
      INTEGER_1 i(4)
c
c     pop count of 32 bit word for sun ... the routine inibit must
c     be called once before popcnt32 or leadz32 is used ... this
c     computes ipopcb and leadzb
c
      ind1 = i(1)
      ind2 = i(2)
      ind3 = i(3)
      ind4 = i(4)
      selci_popcnt32 = ipopcb(ind1) + ipopcb(ind2) +
     $     ipopcb(ind3) + ipopcb(ind4)
c
#endif
      end
      subroutine selci_vpop32(n, input, output)
#if defined(CRAY_T3D) || defined(CRAY_T3E) 
      integer n
      integer input(*)
      integer output(*)
      integer i
c
      do i = 1,n
        output(i) = popcnt(input(i))
      enddo
#elif defined (EXT_INT)
      implicit integer(a-z)
#include "cbits.fh"
c
      INTEGER_1 input(*)
      dimension output(*)
c
      i4 = 1
      do 10 i = 1,n
         output(i) = ipopcb(input(  i4)) + ipopcb(input(1+i4)) +
     $               ipopcb(input(2+i4)) + ipopcb(input(3+i4)) +
     $               ipopcb(input(4+i4)) + ipopcb(input(5+i4)) +
     $               ipopcb(input(6+i4)) + ipopcb(input(7+i4)) 
         i4 = i4 + 8
 10   continue
#else
      implicit integer(a-z)
#include "errquit.fh"
#include "cbits.fh"
c
      INTEGER_1 input(*)
      dimension output(*)
c
      i4 = 1
      do 10 i = 1,n
         output(i) = ipopcb(input(  i4)) + ipopcb(input(1+i4)) +
     $               ipopcb(input(2+i4)) + ipopcb(input(3+i4))
         i4 = i4 + 4
 10   continue
c
#endif
      end
      integer function selci_leadz32(i)
#if defined(CRAY_T3D) || defined(CRAY_T3E) 
      integer i
      selci_leadz32 = leadz(i)
#elif defined (EXT_INT)
c     
#include "cbits.fh"
c     
      INTEGER_1 i(8)
c     
      integer longw
      integer*2 shortw(4)
      equivalence (longw,shortw)
      logical olittleendian
c     
      shortw(1) = 0
      longw = 1
      olittleendian = (shortw(1) .eq. 1)
c     
c     leading zero bit count for 64 bit word on sun ... must
c     call inibit once before using
c     
      if (olittleendian) then
         ind = i(8)
         m = leadzb(ind)
         if (m.eq.8) then
            ind = i(7)
            m = m + leadzb(ind)
            if (m.eq.16) then
               ind = i(6)
               m = m + leadzb(ind)
               if (m.eq.24) then
                  ind = i(5)
                  m = m + leadzb(ind)
                  if (m.eq.32) then
                     ind = i(4)
                     m = m + leadzb(ind)
                     if (m.eq.40) then
                        ind = i(3)
                        m = m + leadzb(ind)
                        if (m.eq.48) then
                           ind = i(2)
                           m = m + leadzb(ind)
                           if (m.eq.56) then
                               ind = i(1)
                               m = m + leadzb(ind)
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
         selci_leadz32 = m
      else
         ind = i(1)
         m = leadzb(ind)
         if (m.eq.8) then
            ind = i(2)
            m = m + leadzb(ind)
            if (m.eq.16) then
               ind = i(3)
               m = m + leadzb(ind)
               if (m.eq.24) then
                  ind = i(4)
                  m = m + leadzb(ind)
                  if (m.eq.32) then
                     ind = i(5)
                     m = m + leadzb(ind)
                     if (m.eq.40) then
                        ind = i(6)
                        m = m + leadzb(ind)
                        if (m.eq.48) then
                           ind = i(7)
                           m = m + leadzb(ind)
                           if (m.eq.56) then
                               ind = i(8)
                               m = m + leadzb(ind)
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
         selci_leadz32 = m
      endif
c     
#else 
c     
#include "cbits.fh"
c     
      INTEGER_1 i(4)
c     
      integer longw
      integer*2 shortw(4)
      equivalence (longw,shortw)
      logical olittleendian
c     
      shortw(1) = 0
      longw = 1
      olittleendian = (shortw(1) .eq. 1)
c     
c     leading zero bit count for 32 bit word on sun ... must
c     call inibit once before using
c     
      if (olittleendian) then
         ind = i(4)
         m = leadzb(ind)
         if (m.eq.8) then
            ind = i(3)
            m = m + leadzb(ind)
            if (m.eq.16) then
               ind = i(2)
               m = m + leadzb(ind)
               if (m.eq.24) then
                  ind = i(1)
                  m = m + leadzb(ind)
               endif
            endif
         endif
         selci_leadz32 = m
      else
         ind = i(1)
         m = leadzb(ind)
         if (m.eq.8) then
            ind = i(2)
            m = m + leadzb(ind)
            if (m.eq.16) then
               ind = i(3)
               m = m + leadzb(ind)
               if (m.eq.24) then
                  ind = i(4)
                  m = m + leadzb(ind)
               endif
            endif
         endif
         selci_leadz32 = m
      endif
c     
#endif
c     
      end
#if defined(CRAY_T3D) || defined(CRAY_T3E) 
      subroutine inibit
      end
#elif defined (EXT_INT)
      subroutine inibit
#include "cbits.fh"
      INTEGER_1 ib(8)
      integer i
      equivalence (ib, i)
      integer bytes_per_int
      logical obigendian
c
#include "bitops.fh"
c
c     initialise lookup arrays for faster pop-count and leading zero co
c
      call util_pack_info(bytes_per_int, obigendian)
      if (bytes_per_int .ne. 8) call errquit
     $     ('sorry ... this selected CI piece is ported to 64 bit',0,
     &       CAPMIS_ERR)
c
      i = 0
      do 10 mm = -128,127,1
         if (obigendian) then
            ib(1) = mm
         else
            ib(8) = mm
         endif
c
c     do pop count using generic 64 bit code
c
         mask = 1
         iset = 0
         do 20 m = 0,63
            if (iand(ishft(i,-m), mask) .ne. 0) iset = iset + 1
 20      continue
         ipopcb(mm) = iset
c
c     do leadz using generic 64 bit code .... note 0 should give 8 not
c
         mask = 1
         do 30 m = 63,0,-1
            if (iand(ishft(i,-m),mask).ne.0) goto 40
 30      continue
 40      leadzb(mm)= 63-m
c
 10   continue
      leadzb(0) = 8
c
c      write(6,*) ' ipopcb '
c      write(6,1) (i,ipopcb(i),i=-128,127)
c 1    format(1x,2i8)
c      write(6,*) ' leadzb '
c      write(6,1) (i,leadzb(i),i=-128,127)
c
      end
#else
      subroutine inibit
#include "cbits.fh"
      INTEGER_1 ib(4)
      integer i
      equivalence (ib, i)
      integer bytes_per_int
      logical obigendian
#ifdef SOLARIS
c     needed with -stackvar option
      save ib, i
#endif
c
#include "bitops.fh"
c
c     initialise lookup arrays for faster pop-count and leading zero co
c
      call util_pack_info(bytes_per_int, obigendian)
      if (bytes_per_int .ne. 4) call errquit
     $     ('sorry ... this selected CI has been ported to 32 bit',0,
     &       CAPMIS_ERR)
c
      i = 0
      do 10 mm = -128,127,1
         if (obigendian) then
            ib(1) = mm
         else
            ib(4) = mm
         endif
c
c     do pop count using generic 32 bit code
c
         mask = 1
         iset = 0
         do 20 m = 0,31
            if (iand(ishft(i,-m), mask) .ne. 0) iset = iset + 1
 20      continue
         ipopcb(mm) = iset
c
c     do leadz using generic 32 bit code .... note 0 should give 8 not
c
         mask = 1
         do 30 m = 31,0,-1
            if (iand(ishft(i,-m),mask).ne.0) goto 40
 30      continue
 40      leadzb(mm)= 31-m
c
 10   continue
      leadzb(0) = 8
c
c      write(6,*) ' ipopcb '
c      write(6,1) (i,ipopcb(i),i=-128,127)
c 1    format(1x,2i8)
c      write(6,*) ' leadzb '
c      write(6,1) (i,leadzb(i),i=-128,127)
c
      end
#endif
*mdc*else
*      subroutine inibit
*      end
*      integer function selci_popcnt32(i)
*      include 'bitops.h'
*c
*c     return bit population count of a 32 bit word
*c
*c     this should work on most 32 bit machines ... but slowly
*c
*      mask = 1
*      iset = 0
*      do 10 m = 0,31
*         if (iand(ishft(i,-m), mask) .ne. 0) iset = iset + 1
* 10   continue
*c
*      selci_popcnt32 = iset
*c
*      end
*      subroutine selci_vpop32(n, input, output)
*      implicit integer(a-z)
*c
*      dimension input(*), output(*)
*c
*      do 10 i = 1,n
*         output(i) = selci_popcnt32(input(i))
* 10   continue
*c
*      end
*      integer function selci_leadz32(i)
*      include 'bitops.h'
*c
*c     return leading zero bit count of 32 bit argument
*c
*c     this should work on most 32 bit machines ... but slowly
*c
*      mask = 1
*      do 10 m = 31,0,-1
*         if (iand(ishft(i,-m),mask).ne.0) goto 20
* 10   continue
*c
* 20   selci_leadz32 = 31-m
*c
*      end
*mdc*endif
