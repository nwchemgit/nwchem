      subroutine gacrea(n,m,name,ihand)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
c
#if defined(DFT_GASTAT)
#include "gastat.fh"
#endif
c
      integer n,m,me
      logical ok
      character*(*) name
      integer ihand,  i
      ok=ga_create(MT_DBL, n, m, name, n, 1, ihand)
      if(.not.ok)  call ga_error(' ga_create failed ',2)
c
#if defined(DFT_GASTAT)
      me = ga_nodeid()
      if(me.eq.0)then
c
c       Tabulate ga_create info.
c
        itab = itab + 1
        table(itab,1) = ihand      
        table(itab,2) = n
        table(itab,3) = m
        table(itab,4) = 1
        names(itab) = name
c
c       Calculate high water mark.
c        
        high = 0
        do i = 1, itab
          if(table(i,4).eq.1)high = high + table(i,2)*table(i,3) 
        enddo
        if(high.gt.highest)then
          highest = high
          do i = 1, itab
            if(table(i,4).eq.1)then
              table(i,5) = 2
            else
              table(i,5) = 0
            endif
          enddo
        endif
      endif
#endif
c
      return
      end
      subroutine gadest(ihand)
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
c
#if defined(DFT_GASTAT)
#include "gastat.fh"
#endif
c
      logical ok
      integer ihand, i, me
      ok=ga_destroy(ihand)
      if(.not.ok)  call ga_error(' ga_destroy failed ',2)
c
#if defined(DFT_GASTAT)
      me = ga_nodeid()
      if(me.eq.0)then
c
c       Tabulate ga_destroy info.
c
        do i = itab, 1, -1
          if(table(i,1).eq.ihand)then
            table(i,4) = 0
            goto 10
          endif
        enddo
   10   continue
      endif
#endif
c
      return
      end
      subroutine gadump
      implicit none
c
#if defined(DFT_GASTAT)
#include "gastat.fh"
c
      integer i, j, ic
      write(6,*)'                GA_STATS                    '
      write(6,*)' At the highwater mark of: ',highest,' MT_DBLs, '
      write(6,*)' the following global arrays were allocated: '
      ic = 0
      do i = 1, itab
        if(table(i,5).eq.2)then
          ic = ic + 1
          write(6,*)ic,names(i),table(i,2),table(i,3)
        endif
      enddo
#endif
c
      return
      end
