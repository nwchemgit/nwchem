      subroutine util_char_ga_brdcst(type, string, origin)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer type
      character*(*) string
      integer origin
c
      integer mitob1
      integer maxs, ls, i
      parameter (maxs = 256)
      integer c(maxs)
c
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
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
      call ga_brdcst(type, ls, mitob1, origin)
      if (ls .gt.len(string)) call errquit
     $     ('util_char_ga_brdcst: receiving string too small',ls,
     &       GA_ERR)
      if (ls .gt. 0) call ga_brdcst(type, c, 
     A     MA_sizeof(MT_INT,ls,MT_BYTE),origin)
c
      string = ' ' 
      do i = 1,ls
         string(i:i) = char(c(i))
      enddo
c
      end

      
