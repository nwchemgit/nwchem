      subroutine speech_input(rtdb)
*
* $Id: speech_input.F,v 1.2 1997-11-04 10:01:45 d3e129 Exp $
*
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
      integer rtdb
c
c     speech <integer hostname> <integer portnumber>
c
c     Identify the server for spoken output.  If the hostname is "off"
c     then disable speech.
c
      character*255 hostname
      integer port
      logical status
c
      if (.not. inp_a(hostname)) call errquit
     $     ('speech_input: failed reading the hostname', 0)
c
      if (inp_compare(.false.,hostname,'off')) then
         status = rtdb_delete(rtdb,'speech:hostname')
         status = rtdb_delete(rtdb,'speech:port')
         return
      endif
c         
      if (.not. inp_i(port)) call errquit
     $     ('speech_input: failed reading the port',0)
c
c     Initialize speech here so that input modules can speak
c
      call util_speak_init(hostname, port)
c
      if (.not. rtdb_cput(rtdb,'speech:hostname',1,hostname))
     $     call errquit('speech_input: failed rtdb_cput hostname',0)
      if (.not. rtdb_put(rtdb,'speech:port',mt_int,1,port))
     $     call errquit('speech_input: failed rtdb_put port',0)
c
      end
