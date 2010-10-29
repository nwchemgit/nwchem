c
c $Id$
c

      SUBROUTINE read_lin(input,ifield,istart,inum_char,irecord,lend)

      implicit none

! This subroutine reads in an input line

      integer i,ifield,iswitch,input,irec_len,istart,inum_char

      character*1 irecord

      logical lcomm,lend

      parameter(irec_len=100)

      dimension istart(irec_len),inum_char(irec_len)
      dimension irecord(irec_len)

      ifield  = 0
      iswitch = 0
      lcomm   = .false.
      lend    = .false.

      read(input,20,end=100)(irecord(i),i=1,irec_len)
20    format(100A1)

      do i=1,irec_len
       if(irecord(i).eq.'#'.and.ifield.eq.0)lcomm = .true.
       if(irecord(i).eq.' ')then
        iswitch = 0
       elseif(iswitch.eq.0)then
        call chg_case(irecord(i))
        ifield = ifield + 1
        istart(ifield) = i
        inum_char(ifield) = 1
        iswitch = 1
       else
        call chg_case(irecord(i))
        inum_char(ifield) = inum_char(ifield) + 1
       endif
      enddo

      goto 200

100   continue

      lend = .true.

200   continue

      return

      END
