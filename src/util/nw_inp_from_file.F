      logical function nw_inp_from_file(rtdb, filename)
*
* $Id: nw_inp_from_file.F,v 1.6 2000-04-10 12:05:34 hvd Exp $
*
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb
      character*(*) filename
c
      logical input_parse, status
      external input_parse
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('nw_inp_from_file',print_low)) then
            write(LuOut,*) ' Reading NWChem input from ', 
     $           filename(1:inp_strlen(filename))
            call util_flush(LuOut)
         endif
         open(4,file=filename,form='formatted',status='old',err=100)
         call inp_save_state()  ! save state of any inp unit
         call inp_init(4,6)
      endif
c
      status = input_parse(rtdb)
c
      if (ga_nodeid() .eq. 0) then
         call inp_restore_state() ! restore state of any inp unit
         close(4,status='keep')
      endif
c
      nw_inp_from_file = .true.
      return
c
 100  call errquit(filename,0)
c
      end
