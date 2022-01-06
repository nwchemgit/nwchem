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
      integer util_sgroup_mygroup
      external util_sgroup_mygroup
      character*(255) filename
      character*(255) unique_name
      integer number

C     This is bad, not 100% sure to be unique, since could be subgroup
      if (ga_pgroup_get_world() .ne. ga_pgroup_get_default()) then
        number = util_sgroup_mygroup()
      else
        number = 0
      end if

      write(unique_name,737) "temp",number,".nw"
 737       format(A4,I4.4,A3)
      call util_file_name(unique_name(1:inp_strlen(unique_name)),
     &            .true., .false., filename)

c      write (6,*) 'filename=',filename
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
* $Id$
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
c         if (util_print('nw_inp_from_file',print_low)) then
c            write(LuOut,*) ' Reading NWChem input from ', 
c     $           filename(1:inp_strlen(filename))
c            call util_flush(LuOut)
c         endif
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
      subroutine pop_inp_string()
*
* $Id: pop_inp_string.F 19707 2010-10-29 17:59:36Z d3y133 $
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      logical lopen
      if (ga_nodeid() .eq. 0) then
c        call inp_restore_state() ! restore state of any inp unit
         inquire(unit=4,opened=lopen)
         if (lopen) then
           close(4,err=102)
         endif
      endif
      return
c
 102  call errquit('pop_inp_string: error closing file',0, INPUT_ERR)
      end
      logical function push_inp_string(string)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"

      character*(*) string
c
      call pop_inp_string()
      if (ga_nodeid() .eq. 0) then
         open(4,form='formatted',status='scratch',err=100)
         write(4,'(a)',err=101) string
         rewind(4,err=101)
c        call inp_save_state()  ! save state of any inp unit
         call inp_init(4,LuOut)
      endif
c
      push_inp_string = inp_read()
      return
c
 100  call errquit('push_inp_string: error opening file',0, INPUT_ERR)
 101  call errquit('push_inp_string: error writing file',0, INPUT_ERR)
c
      end
