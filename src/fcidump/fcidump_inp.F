      subroutine fcidump_input(rtdb)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer rtdb
c
c
      integer num_dirs, ind, mlen, iline, n
      parameter (num_dirs=3)
c
      logical molecorb
c      
      character*15 dirs(num_dirs)
      character*255 test
      data dirs /'orbitals','freeze',
     &     'end'
     &     /
c
c     Grab the current input line and dump the functional specification
c     onto the RTDB. This way external programs can pick the actual
c     functional definition up as it was entered.
c
 10   if (.not. inp_read()) call errquit
     $     ('fcidump_inp: inp_read failed', 1, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('dft_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('dft_input: unknown directive', 0,
     &     INPUT_ERR)
c     
      goto (100,  200,
     $     9999) ind

      call errquit('fcidump_int: unimplemented directive', ind,
     &       INPUT_ERR)
c
c     orbitals atomic/molecular atomic is the default
c             
 100  if (inp_a(test)) then
         if (inp_compare(.false., test, 'molecular')) then
            molecorb=.true.
         else if (inp_compare(.false., test, 'atomic')) then
            molecorb=.false.
         else
            write(LuOut,*)' fcidump_inp:molecular value not recognized'
            call errquit('missing keyword for molecular ',0,INPUT_ERR)
         endif
         if (.not. rtdb_put(rtdb, 'fcidump:mos', mt_log, 1, molecorb))
     &        call errquit('fcidump_inp: rtdb_put failed',
     &        1200    , RTDB_ERR)
      else
         write(LuOut,*)' fcidump_inp: molecular value not found',
     T        test,'@'
         call errquit('missing keyword for molecular ',0,INPUT_ERR)
      endif
      goto 10
 200  continue
      call freeze_input(rtdb,'fcidump')
      goto 10
 9999 continue

      return
      end
