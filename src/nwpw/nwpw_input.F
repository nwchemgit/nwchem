c
c $Id: nwpw_input.F,v 1.5 2002-01-04 23:19:59 bylaska Exp $
c
      subroutine nwpw_input(rtdb)
      implicit none

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c 
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 50)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 'simulation_cell',
     >            'psp_generator',
     >            'psp_formatter',
     >            'wavefunction_initializer',
     >            'wavefunction_converter',
     >            'steepest_descent',
     >            'conjugate_gradient',
     >            'Car-Parrinello',
     >            'v_wavefunction_initializer',
     >            'wavefunction_expander',
     >            'fractional_to_cartesian:',
     >            'cartesian_to_fractional:',
     >            'analysis',
     >            'cell_optimize',
     >            'qmmm',
     >            'dplot',
     >            'band',
     >            'brillouin_zone',
     >            'cell_name:',
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
     >            'mulliken',
     >            'allow_translation',
     >            'mult',
     >            'dft',
     >            'odft',
     >            'restricted',
     >            'unrestricted',
     >            'vectors',
     >            'xc',
     >            'cell_name',
     >            'fake_mass',
     >            'time_step',
     >            'loop',
     >            'tolerances',
     >            'energy_cutoff',
     >            'wavefunction_cutoff',
     >            'ewald_rcut',
     >            'ewald_ncut',
     >            'exchange_correlation',
     >            'end'/

       logical      value
       character*30 cell_name
       character*30 input_wavefunction_filename
       character*30 output_wavefunction_filename
       character*30 exchange_correlation
       double precision fake_mass,time_step,rcut
       integer loop(2),ncut
       double precision tolerances(3),ecut,wcut,fe
       logical mulliken,allow_translation
       integer ispin,multiplicity
       

*     ***** initializations ****
      call nwpw_default(rtdb)

c
c     Any input causes convergence flag to be deleted
c
c
c
cedo#ifdef HPUX
cedo      call errquit(' nwpw not running yet on HP ',0)
cedo#endif
 10   continue
      if (.not. inp_read()) 
     >     call errquit('nwpw_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit('nwpw_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('nwpw_input: unknown directive', 0)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800,
     >       900, 1000,1100,1200,1300,1400,1500,1600,
     >       1700,1800,1900,2000,2100,2200,2300,2400,
     >       2500,2600,2700,2800,2900,3000,3100,3200,
     >       3300,3400,3500,3400,3500,3600,3000,
     >       1900,2200,2300,2400,2500,2600,2700,2800,
     >       2900,3000,
     >       9999) ind
      call errquit('nwpw_input: unimplemented directive', ind)


c     
c  simululation_cell 
c     
 100  call simulation_cell_input(rtdb)
      goto 10


c
c psp_generator
c
 200  call psp_generator_input(rtdb) 
      goto 10


c
c  psp_formatter
c 
 300  call psp_formatter_input(rtdb)
      goto 10


c
c  wavefunction_initializer
c
 400  call wvfnc_init_input(rtdb)
      goto 10


c
c  wavefunction_converter
c
 500  write(*,*) "wavefunction converter"
      write(*,*) "not finished"
      goto 10


c
c  steepest_descent
c
 600  call cpsd_input(rtdb)
      goto 10


c
c  conjugate gradient
c
 700  call cgsd_input(rtdb)
      goto 10


c
c  Car-Parrinello
c
 800  call cpmd_input(rtdb)
      goto 10

c
c  wavefunction_initializer
c
 900  call v_wvfnc_init_input(rtdb)
      goto 10


 1000 continue
       call wvfnc_expander_input(rtdb)
      goto 10

c
c fractional_to_cartesian
c
 1100 call fractional_to_cartesian(rtdb) 
      goto 10

c
c cartesian_to_fractional
c
 1200 call cartesian_to_fractional(rtdb) 
      goto 10

c
c analysis
c
 1300 call analysis_input(rtdb) 
      goto 10

c
c cell_optimize
c
 1400 call cell_optimize_input(rtdb) 
      goto 10

c
c qmmm
c
 1500 call pspw_qmmm_input(rtdb) 
      goto 10

c
c dplot
c
 1600 call pspw_dplot_input(rtdb) 
      goto 10

c
c band
c
 1700 call band_input(rtdb) 
      goto 10

c
c brillouin_zone
c
 1800 call brillouin_zone_input(rtdb) 
      goto 10

**************************************************
****************** cgsd/band input ***************
**************************************************
c     
c  cell_name
c     
 1900 if (.not. inp_a(cell_name))
     >     call errquit(
     >          'nwpw_input: failed to read cell_name', 0)
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))
      value = value.and.
     >        rtdb_cput(rtdb,'band:cell_name',1,cell_name(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1900)
      goto 10

c
c  input_wavefunction_filename
c 
 2000 if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'nwpw_input: failed to read input_wavefunction', 0)
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
     >   .and.rtdb_cput(rtdb,'band:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
*     *********************************
*     **** used by task_save_state ****
*     *********************************
      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cg_pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
     >   .and.rtdb_cput(rtdb,'band:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2000)
      goto 10

c
c  output_wavefunction_filename
c 
 2100 if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'nwpw_input: failed to read output_wavefunction', 0)
      ind = index(output_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cgsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
     >   .and.rtdb_cput(rtdb,'band:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2100)
      goto 10
c
c  fake_mass
c 
 2200 if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'nwpw_input: failed to fake mass', 0)
      value = rtdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
     >   .and.rtdb_put(rtdb,'band:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2200)
      goto 10
c
c  time_step
c 
 2300 if (.not. inp_f(time_step))
     >     call errquit(
     >          'nwpw_input: failed to time step', 0)
      value = rtdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
     >   .and.rtdb_put(rtdb,'band:time_step',mt_dbl,1,time_step)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2300)
      goto 10
c
c  loop
c 
 2400 if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'nwpw_input: failed to read loop(1)', 0)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'nwpw_input: failed to read loop(2)', 0)
      value = rtdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
     >   .and.rtdb_put(rtdb,'band:loop',mt_int,2,loop)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2400)
      goto 10
c
c  tolerances
c
 2500 if (.not. inp_f(fe)) fe=1.0d-7
      tolerances(1) = fe
      if (.not. inp_f(fe)) fe=tolerances(1)
      tolerances(2) = fe
      if (.not. inp_f(fe)) fe=0.0d0
      tolerances(3) = fe
 2501 value = rtdb_put(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
     >   .and.rtdb_put(rtdb,'band:tolerances',mt_dbl,3,tolerances)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2500)
      goto 10
c
c  energy_cutoff
c 
 2600 if (.not. inp_f(ecut))
     >     call errquit(
     >          'nwpw_input: failed to read ndown', 0)
      value = rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2600)
      goto 10

c
c  wavefunction_cutoff
c 
 2700 if (.not. inp_f(wcut))
     >     call errquit(
     >          'nwpw_nwpwt: failed to read ndown', 0)
      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2700)
      goto 10
c
c  ewald_rcut
c 
 2800 if (.not. inp_f(rcut))
     >     call errquit(
     >          'nwpw_input: failed to read rcut', 0)
      value = rtdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
     >   .and.rtdb_put(rtdb,'band:rcut',mt_dbl,1,rcut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2800)
      goto 10
c
c  ewald_rcut
c 
 2900 if (.not. inp_i(ncut))
     >     call errquit(
     >          'nwpw_input: failed to read ncut', 0)
      value = rtdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
     >   .and.rtdb_put(rtdb,'band:ncut',mt_int,1,ncut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2900)
      goto 10
c
c  exchange_correlation
c 
 3000  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'nwpw_input: failed to read exchange_correlation', 0)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
     >   .and.rtdb_cput(rtdb,'band:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3000)
      goto 10
c
c  Mulliken
c 
 3100 mulliken = .true.
      value = rtdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'band:mulliken',mt_log,1,
     >                       mulliken)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3100)
      goto 10
c
c  allow_translation
c 
 3200 allow_translation = .true.
      value = rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,
     >                       allow_translation)
     >   .and.rtdb_put(rtdb,'band:allow_translation',mt_log,1,
     >                       allow_translation)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3200)
      goto 10
c
c  multiplicity
c 
 3300 if (.not. inp_i(multiplicity))
     >     call errquit(
     >          'nwpw_input: failed to read mult', 0)
      if (multiplicity.lt.1) then
         multiplicity = 1
         ispin        = 1
      else if (multiplicity.gt.1) then
         ispin        = 2
      end if
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'cgsd:mult',mt_int,1,multiplicity)
     >   .and.rtdb_put(rtdb,'band:mult',mt_int,1,multiplicity)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3300)
      goto 10
c
c  dft/restricted
c 
 3400 ispin = 1
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3400)
      goto 10
c
c  odft/unrestricted
c 
 3500 ispin = 2
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3500)
      goto 10
c
c     vectors; get input vectors, output vectors, initial guess, etc.
c
 3600 call nwpw_vectors_input(rtdb,'cg_pspw')

      goto 10


************************************************************************
c
c END/Finished
c
 9999 continue
      return
      end
