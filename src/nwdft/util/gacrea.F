      subroutine gacrea(geom,basis,n,m,name,ihand,how)
C$Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer n,m,me
      integer geom,basis
      logical ok
      character*(*) name,how
      integer ihand,  i
      ok=.false.
      if(how.eq.'row')then
*ga:1:0
        ok=ga_create(MT_DBL, n, m, name, n, 0, ihand)
      elseif(how.eq.'column')then
        ok=ga_create(MT_DBL, n, m, name, 0, m, ihand)
      elseif(how.eq.'atom')then
        ihand=ga_create_atom_blocked(geom, basis, name)
        ok=.true.
      endif

      if(.not.ok)  call ga_error('gacrea: ga_create failed ',2)
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
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
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
*#include "gastat.fh"
c
      integer i, j, ic
      write(LuOut,*)'                GA_STATS                    '
      write(LuOut,*)' At the highwater mark of: ',highest,' MT_DBLs,'
      write(LuOut,*)' the following global arrays were allocated: '
      ic = 0
      do i = 1, itab
        if(table(i,5).eq.2)then
          ic = ic + 1
          write(LuOut,*)ic,names(i),table(i,2),table(i,3)
        endif
      enddo
#endif
c
      return
      end
