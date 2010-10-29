      subroutine cckohn_input(rtdb)
c
c     $Id$
c
c     The main purpose of this subroutine is to parse the input
c     and place any "found" parameters into the rtdb.
c
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c#include "geom.fh"
#include "mafdecls.fh"
c#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
c     
c     read input for CCOHN calculation.  All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
c
      logical status
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 24)
      character*12 dirs(num_dirs)
      character*255 test
      data dirs /'cckohn', 'grid', 'noprint', 'print', 'end',
     *     19*''/
c
c     Any CCKOHN input causes calculation to be marked as not converged
c
      status = rtdb_delete(rtdb,'cckohn:converged')
c
c     Lebedev on by default
c
 10   if (.not. inp_read()) call errquit
     $     ('cckohn_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('cckohn_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('cckohn_input: unknown directive', 0,
     &       INPUT_ERR)
c     
      goto (100,  200,  300, 400, 9999) ind
      call errquit('dft_input: unimplemented directive', ind,
     &       INPUT_ERR)
c
c     cckohn;
c             
 100  continue
c     
      goto 10
c
c     grid;
c             
 200  call grid_input(rtdb, 'cckohn')
c
      goto 10
c     
c     noprint;
c     
 300  call util_print_input(rtdb, 'cckohn')
c
      goto 10
c
c     print;
c
 400  call util_print_input(rtdb, 'cckohn')
c
      goto 10
c
c     sym [off|on] (defaults to on)
c
c 1900 if (.not. inp_a(test)) test = 'on'
c      status = inp_compare(.false.,test,'on')
c      if (.not. rtdb_put(rtdb,'dft:skeleton', MT_LOG, 1, status))
c     &   call errquit('dft_input: rtdb_put failed', 1900)
c      goto 10
c
c     adapt [off|on] (defaults to on)
c
c 2000 if (.not. inp_a(test)) test = 'on'
c      status = inp_compare(.false.,test,'on')
c      if (.not. rtdb_put(rtdb,'dft:adapt', MT_LOG, 1, status))
c     &   call errquit('dft_input: rtdb_put failed', 2000)
c      goto 10
c
c     end;
c             
 9999 continue
      return 
c       
      end
