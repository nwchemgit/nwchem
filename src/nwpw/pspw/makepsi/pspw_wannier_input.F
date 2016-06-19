*
* $Id$
*

      subroutine pspw_wannier_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 7)
      character*30 dirs(num_dirs)
      data dirs / 'old_wavefunction_filename',
     >            'new_wavefunction_filename',
     >            'vectors',
     >            'maxiter',
     >            'time_step',
     >            'algorithm',
     >            'end'/

      character*255 test
      character*50  new_wavefunction_filename
      character*50  old_wavefunction_filename
      integer it
      double precision tstep

*     ***** initializations ****

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'pspw_wannier_input: inp_read failed', 0,0)
      if (.not. inp_a(test))
     >     call errquit(
     >      'pspw_wannier_input: failed to read keyword', 1,0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'pspw_wannier_input: unknown directive', 2,0)


      goto ( 100,  200, 300, 400, 500, 600,
     >      9999) ind
      call errquit(
     >   'pspw_wannier_input: unimplemented directive', ind,0)



c
c old_wavefunction_filename
c
 100  if (.not. inp_a(old_wavefunction_filename))
     >     call errquit(
     >      'pspw_wannier_input: failed to read keyword', 3,0)

      ind = index(old_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'wannier:input_wavefunction_filename',
     >                    1,old_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'pspw_wannier_input: rtdb_cput failed', 7,0)
      goto 10


c
c new_wavefunction_filename
c
 200  if (.not. inp_a(new_wavefunction_filename))
     >     call errquit(
     >      'pspw_wannier_input: failed to read ', 4,0)

      ind = index(new_wavefunction_filename,' ') - 1
      if (.not. rtdb_cput(rtdb,'wannier:output_wavefunction_filename',
     >                    1,new_wavefunction_filename(1:ind)))
     >     call errquit(
     >     'pspw_wannier_input: rtdb_cput failed', 8,0)

      goto 10


c
c     vectors; get input vectors, output vectors, initial guess, etc.
c
 300  call nwpw_vectors_input(rtdb,'wannier')

      goto 10


c
c maxiter
c
 400  if (.not.inp_i(it)) it = 500
      if (.not.rtdb_put(rtdb,'wannier:maxiter',mt_int,1,it))
     >  call errquit('nwpw_input: error writing to rtdb',9,0)

      goto 10

c
c time_step
c
 500  if (inp_f(tstep)) then
        if (.not.rtdb_put(rtdb,'wannier:time_step',mt_dbl,1,tstep))
     >  call errquit('nwpw_input: error writing to rtdb',10,0)
      end if

      goto 10
c
c algorithm
c
 600  if (inp_a(new_wavefunction_filename)) then
         if (inp_compare(.false.,new_wavefunction_filename,
     >                  'Damle-Lin-Ying')) then
            if (.not.rtdb_put(rtdb,'nwpw:wannier2',mt_log,1,.true.))
     >        call errquit('nwpw_input: error writing to rtdb',601,0)
         else
            if (.not.rtdb_put(rtdb,'nwpw:wannier2',mt_log,1,.false.))
     >        call errquit('nwpw_input: error writing to rtdb',602,0)
         end if
      end if

      goto 10


*     ***** add wavefunction to rtdb ****
 9999 continue

      return
      end
