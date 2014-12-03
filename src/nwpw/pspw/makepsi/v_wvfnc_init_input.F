*
* $Id$
*

      subroutine v_wvfnc_init_input(rtdb)
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

      character*25 cell_name
      character*50 v_wavefunction_filename
      integer nup
      integer ndown
      logical restricted
      data dirs / 'cell_name:',
     >            'v_wavefunction_filename:',
     >            'restricted',
     >            'unrestricted',
     >            'up_filling:',
     >            'down_filling:',
     >            'end'/



*     ***** initializations ****
      cell_name             = 'cell_default'
      call util_file_prefix('vmovecs',v_wavefunction_filename)
c     v_wavefunction_filename = ' '

      restricted            = .true.
      nup                   = 0
      ndown                 = 0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'v_wvfnc_init_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'wvfnc_init_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'v_wvfnc_init_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 600,
     >      9999) ind
      call errquit(
     >      'v_wvfnc_init_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'v_wvfnc_init_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c wavefunction_filename
c
 200  if (.not. inp_a(v_wavefunction_filename))
     >     call errquit(
     >          'v_wvfnc_init_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c  restricted
c 
 300  restricted = .true.
      goto 10

c
c  unrestricted
c 
 400  restricted = .false.
      goto 10

c
c  up_filling
c 
 500  if (.not. inp_i(nup))
     >     call errquit(
     >          'v_wvfnc_init_input: failed to read nup', 0, INPUT_ERR)

      goto 10
 
c
c  down_filling
c 
 600  if (.not. inp_i(ndown))
     >     call errquit(
     >          'v_wvfnc_init_input: failed to read ndown', 0,
     &       INPUT_ERR)

      goto 10
 
 

*     ***** add wavefunction to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('v_wvfnc_init_input: nead a cell_name', 0,
     &       INPUT_ERR)
      if (v_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'v_wvfnc_init_input: nead a wavefunction_filename', 0,
     &       INPUT_ERR)


      ind = index(cell_name,' ') - 1
      if (.not. rtdb_cput(rtdb,'v_wvfnc:cell_name',1,cell_name(1:ind)))
     >     call errquit(
     >     'v_wvfnc_init_input: rtdb_cput failed', 0, RTDB_ERR)

      ind = index(v_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'v_wvfnc:v_wavefunction_filename',
     >                    1,v_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'v_wvfnc_init_input: rtdb_cput failed', 0, RTDB_ERR)


      if (.not. rtdb_put(rtdb,'v_wvfnc:nup',
     >                    mt_int,1,nup))
     >     call errquit(
     >     'v_wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'v_wvfnc:ndown',
     >                    mt_int,1,ndown))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'v_wvfnc:restricted',
     >                    mt_log,1,restricted))
     >     call errquit(
     >     'v_wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      return
      end
