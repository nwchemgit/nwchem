      subroutine ana_input(irtdb)
c
c $Id: ana_input.F,v 1.2 1999-01-05 21:46:09 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item,rmsfil,select,fpref
c
      integer irmsd
      integer isel,jsel,numsel
      character*255 atom
c
      irmsd=0
c
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,'Analysis Input Module',40,.true.)
      endif
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,fpref))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      select=fpref(1:index(fpref,' ')-1)//'.sel '
      open(unit=69,file=select(1:index(select,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=69)
      numsel=0
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
      if(inp_compare(.false.,'select',item)) then
      isel=0
      jsel=0
      atom=' '
      if(inp_i(isel)) then
      if(.not.inp_i(jsel)) jsel=isel
      endif
      if(.not.inp_a(atom)) goto 4
    4 continue
      numsel=numsel+1
      write(69,6900) isel,jsel,atom(1:index(atom,' ')-1)
 6900 format(2i7,a)
      if(inp_a(atom)) goto 4
      goto 2
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,1999) item
 1999 format('ana_input: unknown field ',a)
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
      if(numsel.gt.0) then
      close(unit=69)
      else
      close(unit=69,status='delete')
      endif
c
      return
c
 9999 continue
      call errquit('ana_input: error',0)
c
      return
      end
