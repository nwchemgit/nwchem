c***********************************************************************
c     rgener
c***********************************************************************
      subroutine drdy_rgener(string,iend,istrt)
c
c     called by:
c           read5
c     calls:
c           drdy_rline
c
      implicit none ! double precision (a-h,o-z)
#include "drdyP.fh"
      character * 80 string
*::local:
      integer istrt
      integer isect
      integer iend
c
c     read in first keyword in this section
c
      call drdy_rline(string,istrt,isect,iend)
c
      do while (isect.eq.0.and.iend.eq.0)
c TITLE
         if (string(istrt:istrt+4).eq.'title') then
*            READ(5,500) ctitle
            READ(lu_drdy_rline,500) ctitle
c RESTART
* should never appear in nwchem input deck
         else if (string(istrt:istrt+5).eq.'restar') then
              lgs(8) = 1
c RE-START
         else if (string(istrt:istrt+6).eq.'re-star') then
              lgs(8) = 1
c ATOMS
         else if (string(istrt:istrt+4).eq.'atoms') then
            call drdy_ratoms(string,istrt)
         else
            write(6,1200) string(istrt:80)
            call errquit('drdy_rgener: fatal error ',911)
         end if
         call drdy_rline(string,istrt,isect,iend)
      end do
c
500   format(80a1)
1200  format(3x,'Error:  the following is not a valid keyword in the',
     *       ' general section',/A80)
c
      return
      end
