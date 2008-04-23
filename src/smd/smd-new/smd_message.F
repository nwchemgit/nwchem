c
c $Id: smd_message.F,v 1.1 2008-04-23 16:48:29 marat Exp $
c

      block data initial_smd_message_data
      implicit none
#include "smd_message_data.fh"
c     
      data smd_message_n /0/
      data smd_message_string /(smd_message_nmax-1)*" ",
     >                         "message buffer overfull"/


      end

      subroutine smd_message_flush(un)

      implicit none
#include "smd_message_data.fh"
#include "global.fh"
      integer un
c
      integer i
      do i=1,smd_message_n
         if(ga_nodeid().eq.0)
     >   write(un,*) smd_message_string(i)
         smd_message_string(i) = " "
      end do
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
      smd_message_n = smd_message_n + 1
      if(smd_message_n.ge.smd_message_nmax) then
        smd_message_n = smd_message_nmax
        return
      end if
      i = smd_message_n
      smd_message_string(i) = buffer 
      return
      end
