      subroutine ana_input(irtdb)
c
c $Id: ana_input.F,v 1.1 1998-12-29 18:01:03 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item,rmsfil
c
      integer irmsd
c
      irmsd=0
c
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,'Analysis Input Module',40,.true.)
      endif
c
    1 continue
c
      if(.not.inp_read()) call errquit('ana_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
c     solute root mean square devation
c
      if(inp_compare(.false.,'rmsdev',item)) then
      irmsd=1
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'initial',item)) then
      irmsd=2
      elseif(inp_compare(.false.,'average',item)) then
      irmsd=3
      else
      rmsfil=item
      endif
      goto 2
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,9999) item
 9999 format('ana_input: unknown field ',a)
      call errquit('ana_input: error reading input file',0)
      endif
c
      if(.not.rtdb_put(irtdb,'ana:irmsd',mt_int,1,irmsd))
     + call errquit('ana_input: rtdb_put irmsdfailed',0)
      if(irmsd.eq.1) then
      if(.not.rtdb_cput(irtdb,'ana:rmsfil',1,rmsfil))
     + call errquit('ana_input: rtdb_put failed',0)
      endif
c
      return
      end
