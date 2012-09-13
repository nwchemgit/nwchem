C> \ingroup task
C> @{
C>
C> \brief Read the input a the shell task directive
C>
C> NWChem has a facility to invoke any system command in a shell. The 
C> corresponding input is
C> ~~~~
C> task shell [(<integer-range processor = 0>|all)] <string command>
C> ~~~~
C> This routine reads that input to get the code ready to invoke the shell 
C> command.
C>
      subroutine task_shell_input(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer rtdb !< [Input] The RTDB handle
c
c     Read input of the shell task directive
c
c     task shell [(<integer-range processor = 0>|all)] <string command>
c
      character*1023 command
      character*8 test
      integer range(3)
c
      if (inp_n_field().lt.3 .or. inp_n_field().gt.4)
     $     call errquit('task shell: invalid syntax',0, INPUT_ERR)
c
      if (inp_n_field() .eq. 4) then
         if (.not. inp_irange(range(1),range(2),range(3))) then
            if (inp_a(test)) then
               if (inp_compare(.false.,test,'all')) then
                  range(1) = 0
                  range(2) = ga_nnodes()-1
                  range(3) = 1
                  goto 10
               endif
            endif
         else
            goto 10
         endif
         call errquit('task shell: bad format for processor range',0,
     &       INPUT_ERR)
c
 10      call inp_clear_err()
         if (range(1).lt.0 .or. range(1).ge.ga_nnodes() .or.
     $        range(2).lt.0 .or. range(2).ge.ga_nnodes() .or.
     $        range(3).lt.1)
     $        call errquit('task shell: invalid process range',0,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb, 'task:shell range', mt_int, 3, range))
     $        call errquit('task shell: rtdb corrupt',0, RTDB_ERR)
      endif
c
      if (.not. inp_a(command))
     $     call errquit('task shell: failed reading command',0,
     &       INPUT_ERR)
c     
      if (.not. rtdb_cput(rtdb,'task:shell command',1,command))
     $     call errquit('task shell: rtdb corrupt',0, RTDB_ERR)
c
      end
C> @}
c $Id$
