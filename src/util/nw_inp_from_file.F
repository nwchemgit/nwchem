      logical function nw_inp_from_character(rtdb, string)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"

      integer rtdb
      character*(*) string
c
      logical nw_inp_from_file
      external nw_inp_from_file
      character*(255) filename

      call util_file_name('temp.nw', .true., .false., filename)

      write (6,*) 'filename=',filename
c
      if (ga_nodeid() .eq. 0) then
         open(1,file=filename,form='formatted',status='unknown',err=100)
         write(1,*,err=101) string
         close(1,err=102)
      endif
      nw_inp_from_character = nw_inp_from_file(rtdb,filename)
      return
c
 100  call errquit('nw_inp_from_char: error opening file',0, INPUT_ERR)
 101  call errquit('nw_inp_from_char: error writing file',0, INPUT_ERR)
 102  call errquit('nw_inp_from_char: error closing file',0, INPUT_ERR)
c
      end
      logical function nw_inp_from_file(rtdb, filename)
*
* $Id: nw_inp_from_file.F,v 1.9 2003-10-17 22:58:35 carlfahl Exp $
*
      implicit none
#include "errquit.fh"
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
 100  call errquit(filename,0, INPUT_ERR)
c
      end
