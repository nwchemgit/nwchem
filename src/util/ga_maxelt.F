      subroutine ga_maxelt(g_a, value)
C     $Id: ga_maxelt.F,v 1.9 2000-11-03 00:43:27 edo Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
      integer g_a
      double precision value
c     
c     Determine the absolute value of the element with largest 
c     absolute magnitude
c     
      integer ilo, ihi, jlo, jhi
      integer ga_maxelt_doit
      external ga_maxelt_doit
c     
c     Each process first determines max of data local to it
c     
      call ga_sync()
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      value = 0.0d0    ! for nodes with no data
      if (ihi.gt.0 .and. ihi.ge.ilo) then
         call ga_access_callback_release(
     $        g_a, ilo, ihi, jlo, jhi,ga_maxelt_doit, value, 
     $        0, 0, 0, 0, 0,0)
      endif
c     
c     Finally a global max via message passing
c     
      call ga_sync()
      call ga_dgop(Msg_MxElt, value, 1, 'absmax')
      call ga_sync()
c     
      end
      integer function ga_maxelt_doit(
     $     g_a, ilo, ihi, jlo, jhi, a, ld, value, 
     $     idum2, idum3, idum4, idum5, idum6,idum7)
      implicit none
      integer g_a, ilo, ihi, jlo, jhi, ld, 
     $     idum2, idum3, idum4, idum5, idum6,idum7
      double precision a(ilo:ilo+ld-1,jlo:jhi), value
      integer i, j
c
      value = 0.0d0
      do j = jlo, jhi
         do i = ilo, ihi
            value = max(abs(a(i,j)),value)
         end do
      end do
c
      ga_maxelt_doit = 0    ! no update necessary
c
      end
