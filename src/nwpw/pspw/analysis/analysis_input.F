      subroutine analysis_input(rtdb)
c
c $Id: analysis_input.F,v 1.2 2001-11-11 01:29:07 bylaska Exp $
c
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 3)
      integer lorb
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'psp_filename:',
     >            'max_orbital:',
     >            'end'/

      lorb = -1

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'analysis_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >      'analysis_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'analysis_input: unknown directive', 0)


      goto ( 100, 200,
     >      9999) ind
      call errquit(
     >   'analysis_input: unimplemented directive', ind)


c
c  psp_filename
c
 100  if (.not. inp_a(test))
     >     call errquit(
     >          'analysis_input: failed to read psp_filename', 0)

      call analysis_expansion_coef(test,lorb,rtdb)
      lorb = -1
      goto 10
 
c
c  lorb
c
 200  if (.not. inp_i(lorb))
     >        call errquit(
     >        'analysis_input: lorb read failed', 0)
      goto 10

 9999 continue

      return
      end
