c
c $Id$
c

      block data initial_smd_message_data
      implicit none
#include "smd_message_data.fh"
c     
      data smd_message_n /0/
      data smd_message_string /smd_message_nmax*" "/

      end

      subroutine smd_message_flush(un)

      implicit none
#include "smd_message_data.fh"
#include "stdio.fh"
#include "global.fh"
      integer un
c
      integer i
      do i=1,smd_message_n
         if(ga_nodeid().eq.0)
     >   write(un,*) smd_message_string(i)
         smd_message_string(i) = " "
      end do
      if(smd_message_n.eq.smd_message_nmax) then
        if(ga_nodeid().eq.0) then
        write(luout,*) "end of buffer reached"
        write(luout,*) "messages could be lost"
        end if
      end if
      smd_message_n = 0
      return

      end

      subroutine smd_message_post(buffer)

      implicit none
#include "smd_message_data.fh"
#include "global.fh"
      character*(*) buffer
c
      integer i
      if(smd_message_n.eq.smd_message_nmax) return
      smd_message_n = smd_message_n + 1
      i = smd_message_n
      smd_message_string(i) = buffer 
      return
      end
