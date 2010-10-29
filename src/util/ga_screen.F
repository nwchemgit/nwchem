       subroutine ga_screen(g_a, value)
C$Id$
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer g_a
       double precision value
c
c    
c     set all elements whose absolute value is less than
c     value to a hard zero
c
       integer ilo, ihi, jlo, jhi
       integer ga_screenb
       external ga_screenb
c
       call ga_sync()
       call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
       if (ihi.gt.0 .and. jhi.gt.0) then
         call ga_access_callback_release(
     $        g_a, ilo, ihi, jlo, jhi,ga_screenb, value, 
     $        0, 0, 0, 0, 0,0)
       endif
c
       call ga_sync()
c
       end
      integer function ga_screenb(
     $     g_a, ilo, ihi, jlo, jhi, a, ld, value, 
     $     idum2, idum3, idum4, idum5, idum6,idum7)
      implicit none
      integer g_a, ilo, ihi, jlo, jhi, ld, 
     $     idum2, idum3, idum4, idum5, idum6,idum7
      double precision a(ilo:ilo+ld-1,jlo:jhi), value
      integer i, j
      double precision zero
      parameter (zero = 0.0d0)
c
      do j = jlo, jhi
         do i = ilo, ihi
            if (abs(a(i,j)) .lt. value) a(i,j) = zero
         end do
      end do
c
      ga_screenb = 1    ! no update necessary
c
      end
