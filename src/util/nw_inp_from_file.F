      logical function nw_inp_from_character(rtdb, string)
      implicit none
#include "rtdb.fh"
#include "global.fh"
      integer rtdb
      character*(*) string
c
      logical nw_inp_from_file
      external nw_inp_from_file
      character*8 filename
      data filename/'temp.nw'/
      
c
      if (ga_nodeid() .eq. 0) then
         open(1,file=filename,form='formatted',status='unknown',err=100)
         write(1,*,err=101) string
         close(1,err=102)
      endif
      nw_inp_from_character = nw_inp_from_file(rtdb,filename)
      return
c
 100  call errquit('nw_inp_from_char: error opening file',0)
 101  call errquit('nw_inp_from_char: error writing file',0)
 102  call errquit('nw_inp_from_char: error closing file',0)
c
      end
      logical function nw_inp_from_file(rtdb, filename)
*
* $Id: nw_inp_from_file.F,v 1.7 2000-11-08 22:07:09 d3g681 Exp $
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
