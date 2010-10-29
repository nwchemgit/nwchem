      subroutine stpr_datestp(iunit,msg)
* $Id$
c
c prints date stamp and message "msg" to unit "iunit"
c
      implicit none
c::passed
      integer iunit
      character*(*) msg
c::local
      character*26 datedate
c
      call util_date(datedate)
      write(iunit,*)datedate,' ',msg
      call util_flush(iunit)
      end
