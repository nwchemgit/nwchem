      subroutine input_title(rtdb)
C$Id: input_title.F,v 1.3 1995-02-02 23:15:19 d3g681 Exp $
      implicit none
#include "rtdb.fh"
#include "inp.fh"
c
      integer rtdb
      character*255 title
      character*80 blank, underline
      integer lent, i
      integer center
c
      if (.not. inp_read()) goto 1000
      if (.not. inp_line(title)) goto 1000
c
c     Shuffle title down to remove leading blanks
c
      lent = inp_strlen(title)
      do center = 1, lent
         if (title(center:center) .ne. ' ') goto 10
      enddo
 10   do i = center, lent
         title(i-center+1:i-center+1) = title(i:i)
      enddo
      title(lent-center+2:) = ' '
c     
c     Write out centered and underlined title
c
      lent = inp_strlen(title)
      center = (80 - lent) / 2
      blank = ' '
      do i = 1, lent
         underline(i:i) = '-'
      enddo
      write(6,1) blank(1:center),title(1:lent)
      write(6,1) blank(1:center),underline(1:lent)
      write(6,*)
 1    format(1x,a,a)
c
c     Store title in the database
c
      if (.not. rtdb_cput(rtdb, 'title', 1, title))
     $     call errquit('input_title: put of title failed', lent)
c
      return
c
 1000 call errquit('input_title: title ; <title>', 0)
      end
