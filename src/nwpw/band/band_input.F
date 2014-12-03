      subroutine band_input(rtdb)
c
c $Id$
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer rtdb
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 32)
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'cell_name:',
     >            'cell_name',
     >            'vectors',
     >            'fake_mass:',
     >            'fake_mass',
     >            'time_step:',
     >            'time_step',
     >            'loop:',
     >            'loop',
     >            'tolerances:',
     >            'tolerances',
     >            'energy_cutoff:',
     >            'energy_cutoff',
     >            'wavefunction_cutoff:',
     >            'wavefunction_cutoff',
     >            'ewald_rcut:',
     >            'ewald_rcut',
     >            'ewald_ncut:',
     >            'ewald_ncut',
     >            'exchange_correlation:',
     >            'exchange_correlation',
     >            'xc',
     >            'mulliken',
     >            'allow_translation',
     >            'zone_name:',
     >            'zone_name',
     >            'mult',
     >            'dft',
     >            'odft',
     >            'restricted',
     >            'unrestricted',
     >            'end'/

       character*50 cell_name,zone_name
       character*50 input_wavefunction_filename
       character*50 output_wavefunction_filename
       character*50 exchange_correlation
       double precision fake_mass,time_step,rcut
       integer loop(2),ncut
       double precision tolerances(3),ecut,wcut,fe
       logical mulliken,allow_translation
       integer ispin,multiplicity
       

*     *****************************************
*     ***** add wavefunction names to rtdb ****
*     *****************************************
      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)

      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'band:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

*     **** used by task_save_state ****
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'band:input vectors',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'band:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))


*     ***** other initializations ****
      cell_name             = 'cell_default'
      zone_name             = 'zone_default'
      exchange_correlation         = 'vosko'
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
      ispin = 1
      multiplicity = 1
c      spin_orbit=.false.
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'band_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'band_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'band_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,100, 200, 300,300, 400,400, 500,500, 600,600, 
     >       700,700, 800,800, 900,900, 1000,1000, 1100,1100,1100,
     >       1200,1300, 1400,1400,
     >       1500,1600,1700,1600,1700,
     >      9999) ind
      call errquit(
     >      'band_input: unimplemented directive', ind, INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'band_input: failed to read cell_name', 0, INPUT_ERR)
      goto 10


c
c  vectors
c 
 200  write(*,*) "vectors_band_input(rtdb) not written yet"
      goto 10

c
c  fake_mass
c 
 300  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'band_input: failed to fake mass', 0, INPUT_ERR)
      goto 10
 
c
c  time_step
c 
 400  if (.not. inp_f(time_step))
     >     call errquit(
     >          'band_input: failed to time step', 0, INPUT_ERR)
      goto 10
 
 
c
c  loop
c 
 500  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'band_input: failed to read loop(1)', 0, INPUT_ERR)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'band_input: failed to read loop(2)', 0, INPUT_ERR)
      goto 10
 
c
c  tolerances
c 
 600  if (.not. inp_f(fe)) goto 601
      tolerances(1) = fe
      if (.not. inp_f(fe)) goto 601
      tolerances(2) = fe
      if (.not. inp_f(fe)) goto 601
      tolerances(3) = fe
 601  goto 10


c
c  energy_cutoff
c 
 700  if (.not. inp_f(ecut))
     >     call errquit(
     >          'band_input: failed to read ndown', 0, INPUT_ERR)
      goto 10

c
c  wavefunction_cutoff
c 
 800  if (.not. inp_f(wcut))
     >     call errquit(
     >          'band_input: failed to read ndown', 0, INPUT_ERR)
      goto 10

c
c  ewald_rcut
c 
 900  if (.not. inp_f(rcut))
     >     call errquit(
     >          'band_input: failed to read rcut', 0, INPUT_ERR)
      goto 10

c
c  ewald_rcut
c 
 1000 if (.not. inp_i(ncut))
     >     call errquit(
     >          'band_input: failed to read ncut', 0, INPUT_ERR)
      goto 10
c
c  exchange_correlation
c 
 1100  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'band_input: failed to read exchange_correlation', 0,
     &       INPUT_ERR)
      goto 10

c
c  Mulliken
c 
 1200 mulliken = .true.
      goto 10
c
c  allow_translation
c 
 1300 allow_translation = .true.
      goto 10

c     
c  zone_name
c     
 1400 if (.not. inp_a(zone_name))
     >     call errquit(
     >          'band_input: failed to read zone_name', 0, INPUT_ERR)
      goto 10
c
c  multiplicity
c
 1500 if (.not. inp_i(multiplicity))
     >     call errquit(
     >          'band_input: failed to read mult', 0, INPUT_ERR)
      if (multiplicity.lt.1) then
         multiplicity = 1
         ispin        = 1
      else if (multiplicity.gt.1) then
         ispin        = 2
      end if
      goto 10
c
c  dft/restricted
c
 1600 ispin = 1
      goto 10
c
c  odft/unrestricted
c
 1700 ispin = 2
      goto 10



*     ***** add wavefunction to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('band_input: nead a cell_name', 0, INPUT_ERR)


      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'band:cell_name',1,cell_name(1:ind))

      ind = index(zone_name,' ') - 1
      value = rtdb_cput(rtdb,'band:zone_name',1,zone_name(1:ind))


      value = rtdb_put(rtdb,'band:mulliken',mt_log,1,
     >                       mulliken)
      value = rtdb_put(rtdb,'band:allow_translation',mt_log,1,
     >                       allow_translation)
      value = rtdb_put(rtdb,'band:fake_mass',mt_dbl,1,fake_mass)
      value = rtdb_put(rtdb,'band:time_step',mt_dbl,1,time_step)
      value = rtdb_put(rtdb,'band:loop',mt_int,2,loop)
      value = rtdb_put(rtdb,'band:tolerances',mt_dbl,3,tolerances)
      value = rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
      value = rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
      value = rtdb_put(rtdb,'band:rcut',mt_dbl,1,rcut)
      value = rtdb_put(rtdb,'band:ncut',mt_int,1,ncut)

      value = rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
      value = rtdb_put(rtdb,'band:mult',mt_int,1,multiplicity)

      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'band:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      return
      end
