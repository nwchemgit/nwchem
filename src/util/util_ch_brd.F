      subroutine util_char_ga_brdcst(type, string)
      implicit none
#include "inp.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer type
      character*(*) string
c
      integer maxs, ls, i
      parameter (maxs = 256)
      integer c(maxs)
c
      if (ga_nodeid() .eq. 0) then
         ls = inp_strlen(string)
         if (ls .gt. maxs) call errquit
     $        ('util_char_ga_brdcst: lazy robert: string too large', ls)
c     
         do i = 1,ls
            c(i) = ichar(string(i:i))
         enddo
      endif
c     
      call ga_brdcst(type, ls, mitob(1), 0)
      if (ls .gt.len(string)) call errquit
     $     ('util_char_ga_brdcst: receiving string too small',ls)
      call ga_brdcst(type, c, mitob(ls), 0)
c
      string = ' ' 
      do i = 1,ls
         string(i:i) = char(c(i))
      enddo
c
      end

      
