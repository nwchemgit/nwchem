      logical function nw_inp_from_file(rtdb, filename)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
c
      integer rtdb
      character*(*) filename
c
      logical input_parse, status
      external input_parse
c
      if (util_print(print_low,'nw_inp_from_file')) then
         write(6,*) ' Reading NWChem input from ', 
     $        filename(1:inp_strlen(filename))
         call util_flush(6)
      endif
      open(4,file=filename,form='formatted',status='old',err=100)
      call inp_save_state() ! save state of any inp unit
      call inp_init(4,6)
      status = input_parse(rtdb)
      call inp_restore_state() ! restore state of any inp unit
      close(4,status='keep')
      nw_inp_from_file = .true.
      return
c
 100  call errquit(filename,0)
c
      end
