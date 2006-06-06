      subroutine tce_hash(hash,key,offset)
c
c $Id: tce_hash.F,v 1.2 2006-06-06 17:10:59 kowalski Exp $
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer hash(*)
      integer key
      integer offset
      integer length
      integer i
      integer less,middle,more
c
      length = hash(1)
      less = 1
      more = length
   10 continue
        if (more-less.le.4) then
          middle = -1
          do i = less, more
            if (hash(i+1).eq.key) middle = i
          enddo
          if (middle.eq.-1) then
            write(*,*) 'key = ',key
            do i = 1,hash(1)
              write(*,*) i,hash(i+1)
            enddo
            call errquit('tce_hash: key not found',0,UNKNOWN_ERR)
          endif
          goto 20
        endif
        middle = (less + more) / 2
        if (hash(middle+1).eq.key) then
          goto 20
        else if (hash(middle+1).gt.key) then
          more = middle
        else
          less = middle
        endif
        goto 10
   20 continue
      offset = hash(length+middle+1)
      return
      end
c
c
c
c
c
c
c
      subroutine tce_hash_n(hash,key,offset)
c
c $Id: tce_hash.F,v 1.2 2006-06-06 17:10:59 kowalski Exp $
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer hash(*)
      integer key
      integer offset
      integer length
      integer i
c
      length = hash(1)
       offset = -1
       do i = 1,length
        if(key.eq.hash(i+1)) then
         offset=hash(length+i+1)
         goto 10
        end if
       end do
 10    continue
       if(offset.eq.-1) then
         call errquit('tce_hash_n: key not found',0,UNKNOWN_ERR)
       end if
c
      return
      end
