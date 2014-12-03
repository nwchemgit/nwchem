c
c $Id$
c
      subroutine analysis_input(rtdb)
      implicit none
      integer rtdb

c
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 5)
      integer lorb
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'psp_filename:',
     >            'psp_filename',
     >            'max_orbital:', 
     >            'max_orbital',
     >            'end'/

      lorb = -1

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'analysis_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >      'analysis_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'analysis_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,100, 200,200,
     >      9999) ind
      call errquit(
     >   'analysis_input: unimplemented directive', ind, INPUT_ERR)


c
c  psp_filename
c
 100  if (.not. inp_a(test))
     >     call errquit(
     >          'analysis_input: failed to read psp_filename', 0,
     &       INPUT_ERR)

      call analysis_expansion_coef(test,lorb,rtdb)
      lorb = -1
      goto 10
 
c
c  lorb
c
 200  if (.not. inp_i(lorb))
     >        call errquit(
     >        'analysis_input: lorb read failed', 0, INPUT_ERR)
      goto 10

 9999 continue

      return
      end
