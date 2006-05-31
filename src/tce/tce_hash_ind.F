      subroutine tce_hash_ind(k_offset_a,key,indexc)
cccx      subroutine tce_hash_ind(hash,key,indexc)
cccx      subroutine tce_hash_ind(hash,key,offset,indexc)
c
c $Id: tce_hash_ind.F,v 1.1 2006-05-31 23:16:00 kowalski Exp $
c
      implicit none
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "tce.fh"
      integer key
      integer length
      integer i
      integer less,middle,more
      integer indexc
      integer k_offset_a
c
      length=int_mb(k_offset_a)
      less = 1
      more = length
   10 continue
        if (more-less.le.4) then
          middle = -1
          do i = less, more
            if (int_mb(k_offset_a+i).eq.key) middle = i
          enddo
          if (middle.eq.-1) then
            write(*,*) 'key = ',key
            do i=1,int_mb(k_offset_a) 
             write(*,*) i,int_mb(k_offset_a+i)
            enddo
            call errquit('tce_hash: key not found',0,UNKNOWN_ERR)
          endif
          goto 20
        endif
        middle = (less + more) / 2
        if (int_mb(k_offset_a+middle).eq.key) then 
          goto 20
        else if (int_mb(k_offset_a+middle).gt.key) then
          more = middle
        else
          less = middle
        endif
        goto 10
   20 continue
      indexc = middle
      return
      end
