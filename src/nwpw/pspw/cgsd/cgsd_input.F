      subroutine cgsd_input(rtdb)
c
c $Id: cgsd_input.F,v 1.1 2001-08-30 00:25:11 edo Exp $
c
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer rtdb
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 17)
 
      character*30 dirs(num_dirs)
      character*255 test, id

      data dirs / 'cell_name:',
     >            'formatted_filename:',
     >            'input_wavefunction_filename:',
     >            'output_wavefunction_filename:',
     >            'fake_mass:',
     >            'time_step:',
     >            'loop:',
     >            'tolerances:',
     >            'energy_cutoff:',
     >            'wavefunction_cutoff:',
     >            'ewald_rcut:',
     >            'ewald_ncut:',
     >            'exchange_correlation:',
     >            'fractional_coordinates',
     >            'mulliken',
     >            'allow_translation',
     >            'end'/

       character*30 cell_name
       character*30 input_wavefunction_filename
       character*30 output_wavefunction_filename
       character*30 exchange_correlation
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut
       double precision tolerances(3),ecut,wcut,fe
       logical frac_coord,mulliken,allow_translation
       


*     ***** initializations ****
      cell_name             = ' '

      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)
c     input_wavefunction_filename  = ' '
c     output_wavefunction_filename = ' '

      exchange_correlation         = 'vosko'
      frac_coord                   = .false.
      mulliken                     = .false.
      allow_translation            = .false.
      fake_mass = 400000.0d0
      time_step = 5.8d0
      loop(1) = 10
      loop(2) = 1
      tolerances(1) = 1.0d-9
      tolerances(2) = 1.0d-9
      tolerances(3) = 1.0d-4
      ecut=9000.0d0
      wcut=9000.0d0
      rcut = 0.0d0
      ncut = 0
      npsp = 0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cgsd_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cgsd_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cgsd_input: unknown directive', 0)


      goto ( 100, 300, 400, 500, 600, 700, 800, 900,
     >       1000,1100,1200,1300,1400,1500,1600,1700,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cgsd_input: failed to read cell_name', 0)
      goto 10


c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cgsd_input: failed to read psp_filename', 0)
      npsp = npsp + 1
      id = 'cgsd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cgsd_input: rtdb_cput failed', 0)
      goto 10


c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cgsd_input: failed to read input_wavefunction', 0)
      goto 10

c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cgsd_input: failed to read output_wavefunction', 0)
      goto 10


c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'cgsd_input: failed to fake mass', 0)
      goto 10
 
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'cgsd_input: failed to time step', 0)
      goto 10
 
 
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'cgsd_input: failed to read loop(1)', 0)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'cgsd_input: failed to read loop(2)', 0)
      goto 10
 
 
c
c  tolerances
c 
c
c  tolerances
c
 900  if (.not. inp_f(fe)) goto 901
      tolerances(1) = fe
      if (.not. inp_f(fe)) goto 901
      tolerances(2) = fe
      if (.not. inp_f(fe)) goto 901
      tolerances(3) = fe
 901  goto 10


c
c  energy_cutoff
c 
 1000 if (.not. inp_f(ecut))
     >     call errquit(
     >          'cgsd_input: failed to read ndown', 0)
      goto 10

c
c  energy_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cgsd_input: failed to read ndown', 0)
      goto 10

c
c  ewald_rcut
c 
 1200 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cgsd_input: failed to read rcut', 0)
      goto 10

c
c  ewald_rcut
c 
 1300 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cgsd_input: failed to read ncut', 0)
      goto 10
c
c  exchange_correlation
c 
 1400  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cgsd_input: failed to read exchange_correlation', 0)
      goto 10
c
c  fractional_coordinates
c 
 1500 frac_coord = .true.
      goto 10

c
c  Mulliken
c 
 1600 mulliken = .true.
      goto 10
c
c  allow_translation
c 
 1700 allow_translation = .true.
      goto 10


*     ***** add wavefunction to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('cgsd_input: nead a cell_name', 0)
      if (input_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cgsd_input: nead an input wavefunction_filename', 0)
      if (output_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cgsd_input: nead an output wavefunction_filename', 0)
      if (npsp.eq.0)
     >  call errquit(
     >   'cgsd_input: nead to enter formatted pseudopotentials', 0)


      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))

      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

*     *********************************
*     **** used by task_save_state ****
*     *********************************
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cg_pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))


      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))

      value = rtdb_put(rtdb,'cgsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      value = rtdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
      value = rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,
     >                       allow_translation)
      value = rtdb_put(rtdb,'cgsd:npsp',     mt_int,1,npsp)
      value = rtdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      value = rtdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      value = rtdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
      value = rtdb_put(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      value = rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
      value = rtdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
      value = rtdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      return
      end
