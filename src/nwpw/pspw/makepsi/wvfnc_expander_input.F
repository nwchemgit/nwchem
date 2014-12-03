*
* $Id$
*

      subroutine wvfnc_expander_input(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 7)
      integer NMAX              ! maximum wavfunctions
      parameter (NMAX = 100)
 
      character*30 dirs(num_dirs)
      character*255 test

      character*50 new_wavefunction_filename
      character*50 old_wavefunction_filename
      integer ngrid(3)
      data dirs / 'old_wavefunction_filename:',
     >            'new_wavefunction_filename:',
     >            'new_ngrid:',
     >            'old_wavefunction_filename',
     >            'new_wavefunction_filename',
     >            'new_ngrid',
     >            'end'/



*     ***** initializations ****
      old_wavefunction_filename = ' '
      new_wavefunction_filename = ' '
      call util_file_prefix('movecs',old_wavefunction_filename)
      call util_file_prefix('movecs',new_wavefunction_filename)
      ngrid(1) = 0
      ngrid(2) = 0
      ngrid(3) = 0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'wvfnc_expander_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'wvfnc_expander_input: unknown directive', 0,
     &       INPUT_ERR)


      goto ( 100,  200, 300, 100, 200, 300,
     >      9999) ind
      call errquit(
     >   'wvfnc_expander_input: unimplemented directive', ind,
     &       INPUT_ERR)



c
c old_wavefunction_filename
c
 100  if (.not. inp_a(old_wavefunction_filename))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c new_wavefunction_filename
c
 200  if (.not. inp_a(new_wavefunction_filename))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0, INPUT_ERR)
      goto 10

 300  if (.not. inp_i(ngrid(1)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0, INPUT_ERR)
      if (.not. inp_i(ngrid(2)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0, INPUT_ERR)
      if (.not. inp_i(ngrid(3)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0, INPUT_ERR)
      goto 10
 

*     ***** add wavefunction to rtdb ****
 9999 continue
      if (old_wavefunction_filename.eq.' ') 
     >  call errquit(
     >    'wvfnc_expander_input: nead a old_wavefunction_filename', 0,
     &       INPUT_ERR)
      if (new_wavefunction_filename.eq.' ') 
     >  call errquit(
     >  'wvfnc_expander_input: nead a new_wavefunction_filename', 0,
     &       INPUT_ERR)



      ind = index(old_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,old_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)


      ind = index(new_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,new_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'xpndr:ngrid',mt_int,3,ngrid))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_put failed', 0, RTDB_ERR)
     

      return
      end
