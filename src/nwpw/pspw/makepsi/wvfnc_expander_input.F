      subroutine wvfnc_expander_input(rtdb)
c
c $Id: wvfnc_expander_input.F,v 1.1 2001-08-30 00:39:46 edo Exp $
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
      parameter (num_dirs = 4)
      integer NMAX              ! maximum wavfunctions
      parameter (NMAX = 100)
 
      character*30 dirs(num_dirs)
      character*255 test

      character*25 new_wavefunction_filename
      character*25 old_wavefunction_filename
      integer ngrid(3)
      data dirs / 'old_wavefunction_filename:',
     >            'new_wavefunction_filename:',
     >            'new_ngrid:',
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
     >           'wvfnc_expander_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'wvfnc_expander_input: unknown directive', 0)


      goto ( 100,  200, 300, 
     >      9999) ind
      call errquit(
     >   'wvfnc_expander_input: unimplemented directive', ind)



c
c old_wavefunction_filename
c
 100  if (.not. inp_a(old_wavefunction_filename))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read keyword', 0)
      goto 10


c
c new_wavefunction_filename
c
 200  if (.not. inp_a(new_wavefunction_filename))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0)
      goto 10

 300  if (.not. inp_i(ngrid(1)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0)
      if (.not. inp_i(ngrid(2)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0)
      if (.not. inp_i(ngrid(3)))
     >     call errquit(
     >      'wvfnc_expander_input: failed to read ', 0)
      goto 10
 

*     ***** add wavefunction to rtdb ****
 9999 continue
      if (old_wavefunction_filename.eq.' ') 
     >  call errquit(
     >    'wvfnc_expander_input: nead a old_wavefunction_filename', 0)
      if (new_wavefunction_filename.eq.' ') 
     >  call errquit(
     >  'wvfnc_expander_input: nead a new_wavefunction_filename', 0)



      ind = index(old_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'xpndr:old_wavefunction_filename',
     >                    1,old_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0)


      ind = index(new_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'xpndr:new_wavefunction_filename',
     >                    1,new_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_cput failed', 0)

      if (.not. rtdb_put(rtdb,'xpndr:ngrid',mt_int,3,ngrid))
     >     call errquit(
     >     'wvfnc_expander_input: rtdb_put failed', 0)
     

      return
      end
