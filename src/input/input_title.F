      subroutine input_title(rtdb)
C$Id: input_title.F,v 1.7 2001-11-28 18:40:41 bert Exp $
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "stdio.fh"
c
      integer rtdb
      character*255 title
      character*255 blank, underline
      integer lent, i
      integer center
c
c     Detect extraneous input after title directive ... probable error
c
      call inp_set_field(0)
      if (inp_n_field() .gt. 2) goto 1000
      if (.not. inp_a(title)) call errquit
     *      ('input_title: title is not set',555)
      if (.not. inp_a(title)) goto 1000
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
      center = (255 - lent) / 2
      blank = ' '
      do i = 1, lent
         underline(i:i) = '-'
      enddo
      write(LuOut,1) blank(1:center),title(1:lent)
      write(LuOut,1) blank(1:center),underline(1:lent)
      write(LuOut,*)
 1    format(1x,a,a)
c
c     Store title in the database
c
      if (.not. rtdb_cput(rtdb, 'title', 1, title))
     $     call errquit('input_title: put of title failed', lent)
c
      return
c
 1000 call errquit('input_title: title <string title>', 0)
      end
