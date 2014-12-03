*
* $Id$
*


      subroutine wvfnc_init_input(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 17)
      integer NMAX              ! maximum wavfunctions
      parameter (NMAX = 501)
 
      character*22 dirs(num_dirs)
      character*255 test

      character*50 cell_name,wavefunction_filename
      integer   up_filling(4,NMAX),nup
      integer down_filling(4,NMAX),ndown
      logical restricted
      integer i
      data dirs / 'cell_name:',
     >            'wavefunction_filename:',
     >            'restricted',
     >            'unrestricted',
     >            'up_filling:',
     >            'down_filling:',
     >            'up_electrons:',
     >            'down_electrons:',
     >            'restricted_electrons:',
     >            'cell_name',
     >            'wavefunction_filename',
     >            'up_filling',
     >            'down_filling',
     >            'up_electrons',
     >            'down_electrons',
     >            'restricted_electrons',
     >            'end'/



*     ***** initializations ****
      cell_name             = 'cell_default'

      call util_file_prefix('movecs',wavefunction_filename)
c     wavefunction_filename = ' '

      restricted            = .true.
      nup                   = 0
      ndown                 = 0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'wvfnc_init_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'wvfnc_init_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'wvfnc_init_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500,600,700,800,900,
     >       100,  200,           500,600,700,800,900,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'wvfnc_init_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c wavefunction_filename
c
 200  if (.not. inp_a(wavefunction_filename))
     >     call errquit(
     >          'wvfnc_init_input: failed to read keyword', 0,
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
     >          'wvfnc_init_input: failed to read nup', 0, INPUT_ERR)

      do i=1,nup
         if (.not. inp_read()) 
     >        call errquit(
     >           'wvfnc_init_input: inp_read failed', 0, INPUT_ERR)

           if (.not. inp_i(up_filling(1,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(up_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(up_filling(2,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(up_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(up_filling(3,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(up_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(up_filling(4,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(up_filling) failed', 0,
     &       INPUT_ERR)
      end do
      goto 10
 
c
c  down_filling
c 
 600  if (.not. inp_i(ndown))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0, INPUT_ERR)

      do i=1,ndown
         if (.not. inp_read()) 
     >        call errquit(
     >           'wvfnc_init_input: inp_read failed', 0, INPUT_ERR)

           if (.not. inp_i(down_filling(1,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(down_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(down_filling(2,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(down_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(down_filling(3,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(down_filling) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(down_filling(4,i)))
     >     call errquit(
     >          'wvfnc_init_input: inp_i(down_filling) failed', 0,
     &       INPUT_ERR)
      end do
      goto 10
 
c
c  up_electrons
c 
 700  if (.not. inp_i(nup))
     >     call errquit(
     >          'wvfnc_init_input: failed to read nup', 0, INPUT_ERR)
      call getfilling(.true.,nup,up_filling)
      goto 10
 
 
c
c  down_electrons
c 
 800  if (.not. inp_i(ndown))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0, INPUT_ERR)
      call getfilling(.true.,ndown,down_filling)
      goto 10
 
c
c  restricted_electrons
c 
 900  if (.not. inp_i(nup))
     >     call errquit(
     >          'wvfnc_init_input: failed to read nup', 0, INPUT_ERR)
      nup=nup/2
      call getfilling(.true.,nup,up_filling)
      goto 10


*     ***** add wavefunction to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('wvfnc_init_input: nead a cell_name', 0, INPUT_ERR)
      if (wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'wvfnc_init_input: nead a wavefunction_filename', 0,
     &       INPUT_ERR)


      ind = index(cell_name,' ') - 1
      if (.not. rtdb_cput(rtdb,'wvfnc:cell_name',1,cell_name(1:ind)))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_cput failed', 0, RTDB_ERR)

      ind = index(wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'wvfnc:wavefunction_filename',
     >                    1,wavefunction_filename(1:ind)))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_cput failed', 0, RTDB_ERR)


      if (.not. rtdb_put(rtdb,'wvfnc:nup',
     >                    mt_int,1,nup))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'wvfnc:ndown',
     >                    mt_int,1,ndown))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'wvfnc:up_filling',
     >                    mt_int,(4*nup+1),up_filling))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'wvfnc:down_filling',
     >                    mt_int,(4*ndown+1),down_filling))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      if (.not. rtdb_put(rtdb,'wvfnc:restricted',
     >                    mt_log,1,restricted))
     >     call errquit(
     >     'wvfnc_init_input: rtdb_put failed', 0, RTDB_ERR)

      return
      end
