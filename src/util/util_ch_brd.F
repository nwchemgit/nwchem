      subroutine util_char_ga_brdcst(type, string, origin)
*
* $Id: util_ch_brd.F,v 1.4 2003-10-17 22:58:35 carlfahl Exp $
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer type
      character*(*) string
      integer origin
c
      integer maxs, ls, i
      parameter (maxs = 256)
      integer c(maxs)
c
      if (ga_nodeid() .eq. origin) then
         ls = inp_strlen(string)
         if (ls .gt. maxs) call errquit
     $        ('util_char_ga_brdcst: lazy robert: string too large', ls,
     &       GA_ERR)
c     
         do i = 1,ls
            c(i) = ichar(string(i:i))
         enddo
      endif
c     
      call ga_brdcst(type, ls, mitob(1), origin)
      if (ls .gt.len(string)) call errquit
     $     ('util_char_ga_brdcst: receiving string too small',ls,
     &       GA_ERR)
      if (ls .gt. 0) call ga_brdcst(type, c, mitob(ls), origin)
c
      string = ' ' 
      do i = 1,ls
         string(i:i) = char(c(i))
      enddo
c
      end

      
