c
c $Id: nwpw_input.F,v 1.54 2009-03-24 18:04:27 bylaska Exp $
c
      subroutine nwpw_input(rtdb)
      implicit none
      integer rtdb

c 
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 93)
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
     >            'dplot','band_dplot',
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
     >            'wannier',
     >            'geometry_optimize',
     >            'pseudopotential_libraries',
     >            'pseudopotentials',
     >            'dos',
     >            'monkhorst-pack',
     >            'zone_name',
     >            'zone_structure_name',
     >            'dos-grid',
     >            'dos_filename',
     >            'efield',
     >            'occupations',
     >            'expand_cell',
     >            'mapping',
     >            'cg',
     >            'lmbfgs',
     >            'scf',
     >            'rotation',
     >            'translation',
     >            'integrate_mult_l',
     >            'balance',
     >            'nobalance',
     >            'bo_steps',
     >            'bo_time_step',
     >            'bo_algorithm',
     >            'bo_fake_mass',
     >            'scaling',
     >            'xyz_filename',
     >            'np_dimensions',
     >            'spin_orbit',
     >            'spin-orbit',
     >            'cutoff',
     >            'virtual',
     >            'fractional_orbitals',
     >            'tetrahedron',
     >            'smear',
     >            'lcao_skip',
     >            'lcao',
     >            'ewald_ngrid',
     >            'auxiliary_potentials',
     >            'Nose-Hoover',
     >            'band_structure',
     >            'end'/

       double precision kb
       parameter (kb=3.16679d-6)
       logical      value,move,mulliken_kawai,nose,nosers
       character*30 cell_name,zone_name
       character*30 input_wavefunction_filename
       character*30 output_wavefunction_filename
       character*30 exchange_correlation
       double precision fake_mass,time_step,rcut
       integer loop(2),ncut,cell_expand(3)
       double precision tolerances(3),ecut,wcut,fe,ff,scaling(2)
       double precision alpha,emin,emax
       double precision Pe,Pr,Te,Tr
       integer npoints
       logical mulliken,allow_translation
       integer ispin,multiplicity
       integer nx,ny,nz,minimizer
       

*     ***** initializations ****
      call nwpw_default(rtdb)
      ispin = 1

c
c     Any input causes convergence flag to be deleted
c
c
c
cedo#ifdef HPUX
cedo      call errquit(' nwpw not running yet on HP ',0, CAPMIS_ERR)
cedo#endif
 10   continue
      if (.not. inp_read()) 
     >     call errquit('nwpw_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit('nwpw_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('nwpw_input: unknown directive', 0,
     &       INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800,
     >       900, 1000,1100,1200,1300,1400,1500,1600,1650,
     >       1700,1800,1900,2000,2100,2200,2300,2400,
     >       2500,2600,2700,2800,2900,3000,3100,3200,
     >       3300,3400,3500,3400,3500,3600,3000,
     >       1900,2200,2300,2400,2500,2600,2700,2800,
     >       2900,3000,3700,3800,3900,3900,4000,1805,1806,1807,1808,
     >       4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,
     >       5200,5300, 5400,5500,5600,5650,5700,5800,5900,6000,6000,
     >       2750,6090,6100,1809,6200,6300,6400,6500,1500,6600,1807,
     >       9999) ind
      call errquit('nwpw_input: unimplemented directive', ind,
     &       INPUT_ERR)


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
c dplot
c
 1650 call band_dplot_input(rtdb) 
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


c
c Monkhorst-Pack
c
 1805 if (.not. inp_i(nx)) nx = 1
      if (.not. inp_i(ny)) ny = 1
      if (.not. inp_i(nz)) nz = 1
      if (.not. inp_a(zone_name)) zone_name    = 'zone_default'
      
      call monkhorst_pack_set(rtdb,zone_name,nx,ny,nz) 
      goto 10


c
c  zone_name
c
 1806 if (.not. inp_a(zone_name))
     >     call errquit(
     >          'nwpw_input: failed to read zone_name', 0, INPUT_ERR)

      ind = index(zone_name,' ') - 1
      value = rtdb_cput(rtdb,'band:zone_name',1,zone_name(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1806, RTDB_ERR)

      goto 10

c
c  zone_structure_name
c
 1807 if (.not. inp_a(zone_name))
     >     call errquit(
     >          'nwpw_input: failed to read zone_name', 0, INPUT_ERR)

      ind = index(zone_name,' ') - 1
      value = rtdb_cput(rtdb,'band_structure:zone_name',1,
     >                  zone_name(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1806, RTDB_ERR)

      goto 10

c
c DOS-Grid
c
 1808 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 2
      if (.not. inp_i(cell_expand(2))) cell_expand(2) = 2
      if (.not. inp_i(cell_expand(3))) cell_expand(3) = 2
      if (.not. inp_a(zone_name)) zone_name    = 'structure_default'

      call dos_states_set(rtdb,zone_name,cell_expand(1), 
     >                                   cell_expand(2),
     >                                   cell_expand(3))
      value = rtdb_put(rtdb,'band:dos-grid',mt_int,3,cell_expand)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1808, RTDB_ERR)
      goto 10

c
c tetrahedron
c
 1809 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 2
      if (.not. inp_i(cell_expand(2))) cell_expand(2) = 2
      if (.not. inp_i(cell_expand(3))) cell_expand(3) = 2
      if (.not. inp_a(zone_name)) zone_name    = 'zone_default'

      call dos_states_set(rtdb,zone_name,cell_expand(1),
     >                                   cell_expand(2),
     >                                   cell_expand(3))
      goto 10


cc
cc DOS-structure
cc
c 1809 value = rtdb_put(rtdb,'band:dos-structure',mt_log,1,.true.)
c      if (.not.value)
c     >  call errquit(
c     >       'nwpw_input: error writing to rtdb', 1809, RTDB_ERR)
c      goto 10
cc
cc band-structure
cc
c 1810 value = rtdb_put(rtdb,'band:dos-structure',mt_log,1,.false.)
c      if (.not.value)
c     >  call errquit(
c     >       'nwpw_input: error writing to rtdb', 1809, RTDB_ERR)
c      goto 10



*******************************************************
****************** cgsd/band/cpsd input ***************
*******************************************************
c     
c  cell_name
c     
 1900 if (.not. inp_a(cell_name))
     >     call errquit(
     >          'nwpw_input: failed to read cell_name', 0, INPUT_ERR)
      ind = index(cell_name,' ') - 1
      value =      rtdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))
     >        .and.rtdb_cput(rtdb,'band:cell_name',1,cell_name(1:ind))
     >        .and.rtdb_cput(rtdb,'cpsd:cell_name',1,cell_name(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1900, RTDB_ERR)
      goto 10

c
c  input_wavefunction_filename
c 
 2000 if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'nwpw_input: failed to read input_wavefunction', 0,
     &       INPUT_ERR)
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
     >        rtdb_cput(rtdb,'pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
     >   .and.rtdb_cput(rtdb,'band:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2000, RTDB_ERR)
      goto 10

c
c  output_wavefunction_filename
c 
 2100 if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'nwpw_input: failed to read output_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(output_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cgsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
     >   .and.rtdb_cput(rtdb,'band:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2100, RTDB_ERR)
      goto 10
c
c  fake_mass
c 
 2200 if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'nwpw_input: failed to fake mass', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
     >   .and.rtdb_put(rtdb,'band:fake_mass',mt_dbl,1,fake_mass)
     >   .and.rtdb_put(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2200, RTDB_ERR)
      goto 10
c
c  time_step
c 
 2300 if (.not. inp_f(time_step))
     >     call errquit(
     >          'nwpw_input: failed to time step', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
     >   .and.rtdb_put(rtdb,'band:time_step',mt_dbl,1,time_step)
     >   .and.rtdb_put(rtdb,'cpsd:time_step',mt_dbl,1,time_step)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2300, RTDB_ERR)
      goto 10
c
c  loop
c 
 2400 if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'nwpw_input: failed to read loop(1)', 0, INPUT_ERR)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'nwpw_input: failed to read loop(2)', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
     >   .and.rtdb_put(rtdb,'band:loop',mt_int,2,loop)
     >   .and.rtdb_put(rtdb,'cpsd:loop',mt_int,2,loop)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2400, RTDB_ERR)
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
     >       'nwpw_input: error writing to rtdb', 2500, RTDB_ERR)
      goto 10
c
c  energy_cutoff
c 
 2600 if (.not. inp_f(ecut))
     >     call errquit(
     >          'nwpw_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2600, RTDB_ERR)
      goto 10

c
c  wavefunction_cutoff
c 
 2700 if (.not. inp_f(wcut))
     >     call errquit(
     >          'nwpw_nwpwt: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2700, RTDB_ERR)
      goto 10

c
c  cutoff
c
 2750 if (.not. inp_f(wcut))
     >     call errquit(
     >          'nwpw_nwpwt: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,2.0d0*wcut)
     >   .and.rtdb_put(rtdb,'band:ecut',mt_dbl,1,2.0d0*wcut)
     >   .and.rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,2.0d0*wcut)
     >   .and.rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,2.0d0*wcut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2750, RTDB_ERR)
      goto 10


c
c  ewald_rcut
c 
 2800 if (.not. inp_f(rcut))
     >     call errquit(
     >          'nwpw_input: failed to read rcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
     >   .and.rtdb_put(rtdb,'band:rcut',mt_dbl,1,rcut)
     >   .and.rtdb_put(rtdb,'cpsd:rcut',mt_dbl,1,rcut)
     >   .and.rtdb_put(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2800, RTDB_ERR)
      goto 10
c
c  ewald_rcut
c 
 2900 if (.not. inp_i(ncut))
     >     call errquit(
     >          'nwpw_input: failed to read ncut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
     >   .and.rtdb_put(rtdb,'band:ncut',mt_int,1,ncut)
     >   .and.rtdb_put(rtdb,'cpsd:ncut',mt_int,1,ncut)
     >   .and.rtdb_put(rtdb,'cpmd:ncut',mt_int,1,ncut)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 2900, RTDB_ERR)
      goto 10
c
c  exchange_correlation
c 
 3000  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'nwpw_input: failed to read exchange_correlation', 0,
     &       INPUT_ERR)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
     >   .and.rtdb_cput(rtdb,'band:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
     >   .and.rtdb_cput(rtdb,'cpsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
     >   .and.rtdb_cput(rtdb,'cpmd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.false.)
      value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.false.)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,1.0d0)

      if ((inp_compare(.false.,exchange_correlation,'lda-sic')).or.
     >    (inp_compare(.false.,exchange_correlation,'pbe96-sic'))) then
       if (.not.inp_f(fe)) fe = 1.0d0
       if (.not.inp_f(ff)) ff = fe
       value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,fe)
       value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,ff)

      end if


      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3000, RTDB_ERR)
      goto 10
c
c  Mulliken
c 
 3100 if (.not.inp_a(zone_name)) zone_name = 'lcao'

      mulliken       = .true.
      mulliken_kawai = .false.
      ind = index(zone_name,' ') - 1
      if (inp_compare(.false.,zone_name,'kawai')) mulliken_kawai=.true.

      if (inp_compare(.false.,zone_name,'off')) then
          mulliken      =.false.
          mulliken_kawai=.false.
      end if

      value = rtdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'band:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'nwpw:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'cpsd:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'cpmd:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'nwpw:mulliken_kawai',mt_log,1,
     >                       mulliken_kawai)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3100, RTDB_ERR)
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
     >       'nwpw_input: error writing to rtdb', 3200, RTDB_ERR)
      goto 10
c
c  multiplicity
c 
 3300 if (.not. inp_i(multiplicity))
     >     call errquit(
     >          'nwpw_input: failed to read mult', 0, INPUT_ERR)
      if (multiplicity.lt.1) then
         multiplicity = 1
         ispin        = 1
      else if (multiplicity.gt.1) then
         ispin        = 2
      end if
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'cgsd:mult',mt_int,1,multiplicity)
     >   .and.rtdb_put(rtdb,'band:mult',mt_int,1,multiplicity)
     >   .and.rtdb_put(rtdb,'cpsd:mult',mt_int,1,multiplicity)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3300, RTDB_ERR)
      goto 10
c
c  dft/restricted
c 
 3400 ispin = 1
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3400, RTDB_ERR)
      goto 10
c
c  odft/unrestricted
c 
 3500 ispin = 2
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'band:ispin',mt_int,1,ispin)
     >   .and.rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3500, RTDB_ERR)
      goto 10
c
c     vectors; get input vectors, output vectors, initial guess, etc.
c
 3600 call nwpw_vectors_input(rtdb,'pspw')

      goto 10

c
c wannier
c
 3700 call pspw_wannier_input(rtdb) 
      goto 10

c
c  geometry_optimize
c 
 3800 move = .true.
      value =     rtdb_put(rtdb,'cgsd:geometry_optimize',mt_log,1,move)
     >       .and.rtdb_put(rtdb,'cpsd:geometry_optimize',mt_log,1,move)
     >       .and.rtdb_put(rtdb,'band:geometry_optimize',mt_log,1,move)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 3800, RTDB_ERR)
      goto 10


c
c  pseudopotential_library
c
 3900 call nwpw_set_library(rtdb)

      goto 10


c
c  DOS - density of states input
c
 4000 if (.not.inp_f(alpha)) alpha = 0.05d0/27.2116d0
      value = rtdb_put(rtdb,'dos:alpha',mt_dbl,1,alpha)

      if (inp_i(npoints))
     >  value = value.and.rtdb_put(rtdb,'dos:npoints',mt_int,1,npoints)
      if (inp_f(emin))
     >  value = value.and.rtdb_put(rtdb,'dos:emin',mt_dbl,1,emin)
      if (inp_f(emax))
     >  value = value.and.rtdb_put(rtdb,'dos:emax',mt_dbl,1,emax)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 4000, RTDB_ERR)

      goto 10




c
c  DOS_FILENAME - density of states input
c
 4100 if (inp_a(zone_name)) then
        ind = index(zone_name,' ') - 1
        if (.not.rtdb_cput(rtdb,'nwpw:dos:filename',1,zone_name(1:ind)))
     >  call errquit('nwpw_input: error writing to rtdb',4100,RTDB_ERR)
      end if

      goto 10


c
c  EFIELD
c
 4200 if (.not.rtdb_put(rtdb,'nwpw:efield',mt_log,1,.true.))
     >  call errquit('nwpw_input: error writing to rtdb',4200,RTDB_ERR)

      goto 10

c
c OCCUPATIONS
c
 4300 call frac_occ_input(rtdb)
      goto 10

c
c  expand_cell
c
 4400 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 1
      if (.not. inp_i(cell_expand(2))) cell_expand(2) = 1
      if (.not. inp_i(cell_expand(3))) cell_expand(3) = 1
      value = rtdb_put(rtdb,'nwpw:cell_expand',mt_int,3,cell_expand)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 4400, RTDB_ERR)
      goto 10

c
c  mapping
c
 4500 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 1
      value = rtdb_put(rtdb,'nwpw:mapping',mt_int,1,cell_expand)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 4500, RTDB_ERR)
      goto 10

c
c  CG minimizer
c
 4600 nx = 1
 4601 if (.not.inp_a(zone_name)) goto 4602
       if (inp_compare(.false.,zone_name,'stiefel'))  nx=4
       if (inp_compare(.false.,zone_name,'grassman')) nx=1
      go to 4601

 4602 if (.not.rtdb_put(rtdb,'nwpw:minimizer',mt_int,1,nx))
     >  call errquit('nwpw_input: error writing to rtdb',4600,RTDB_ERR)

      goto 10

c
c  LMBFGS minimizer
c
 4700 nx=2
 4701 if (.not.inp_a(zone_name)) goto 4702
       if (inp_compare(.false.,zone_name,'stiefel'))  nx=7
       if (inp_compare(.false.,zone_name,'grassman')) nx=2
      go to 4701

 4702 if (.not.rtdb_put(rtdb,'nwpw:minimizer',mt_int,1,nx))
     >  call errquit('nwpw_input: error writing to rtdb',4700,RTDB_ERR)

      goto 10


c
c  SCF minimizer
c
 4800 minimizer = 8 !default is density mixing
 4801 if (.not.inp_a(zone_name)) goto 4802
       if (inp_compare(.false.,zone_name,'potential')) minimizer = 5
       if (inp_compare(.false.,zone_name,'density'))   minimizer = 8
       if (inp_compare(.false.,zone_name,'cg')) then
         if(.not.rtdb_put(rtdb,'nwpw:ks_algorithm',mt_int,1,0))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'rmm-diis')) then
         if(.not.rtdb_put(rtdb,'nwpw:ks_algorithm',mt_int,1,1))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'simple')) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,0))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'broyden')) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,1))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'diis')) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,2))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'anderson')) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,3))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if

       if (inp_compare(.false.,zone_name,'alpha')) then
          if (.not.inp_f(alpha)) 
     >      call errquit('nwpw_input: error reading scf alpha',
     >                   4800,RTDB_ERR)
          if(.not.rtdb_put(rtdb,'nwpw:ks_alpha',mt_dbl,1,alpha))
     >    call errquit('nwpw_input: error writing to rtdb',
     >                 4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'iterations')) then
          if (.not.inp_i(nx)) 
     >      call errquit('nwpw_input: error reading iterations',
     >                   4800,RTDB_ERR)
          if(.not.rtdb_put(rtdb,'nwpw:ks_maxit_orb',mt_int,1,nx))
     >    call errquit('nwpw_input: error writing to rtdb',
     >                 4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'outer_iterations')) then
          if (.not.inp_i(nx)) 
     >      call errquit('nwpw_input: error reading iterations',
     >                   4800,RTDB_ERR)
          if(.not.rtdb_put(rtdb,'nwpw:ks_maxit_orbs',mt_int,1,nx))
     >    call errquit('nwpw_input: error writing to rtdb',
     >                 4800,RTDB_ERR)
       end if

      goto 4801

 4802 if (.not.rtdb_put(rtdb,'nwpw:minimizer',mt_int,1,minimizer))
     >  call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)

      goto 10


c
c  rotation
c
 4900 if (.not.inp_a(zone_name)) goto 4901
       if (inp_compare(.false.,zone_name,'on'))     move = .true.
       if (inp_compare(.false.,zone_name,'yes'))    move = .true.
       if (inp_compare(.false.,zone_name,'.true.')) move = .true.
       if (inp_compare(.false.,zone_name,'T'))      move = .true.
       if (inp_compare(.false.,zone_name,'off'))     move = .false.
       if (inp_compare(.false.,zone_name,'no'))      move = .false.
       if (inp_compare(.false.,zone_name,'.false.')) move = .false.
       if (inp_compare(.false.,zone_name,'F'))       move = .false.
       if (inp_compare(.false.,zone_name,'fix'))     move = .false.
      if (.not.rtdb_put(rtdb,'nwpw:rotation',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',4900,RTDB_ERR)

 4901 continue
      goto 10


c
c  translation
c
 5000 if (.not.inp_a(zone_name)) goto 5001
      if (inp_compare(.false.,zone_name,'on'))     move = .true.
      if (inp_compare(.false.,zone_name,'yes'))    move = .true.
      if (inp_compare(.false.,zone_name,'.true.')) move = .true.
      if (inp_compare(.false.,zone_name,'T'))      move = .true.
      if (inp_compare(.false.,zone_name,'off'))     move = .false.
      if (inp_compare(.false.,zone_name,'no'))      move = .false.
      if (inp_compare(.false.,zone_name,'.false.')) move = .false.
      if (inp_compare(.false.,zone_name,'F'))       move = .false.
      if (inp_compare(.false.,zone_name,'fix'))     move = .false.
      if (.not.rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',5000,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'band:allow_translation',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',5000,RTDB_ERR)

 5001 continue
      goto 10

c
c  integrate_mult_l
c
 5100 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 0
      value = rtdb_put(rtdb,'nwpw:lmax_multipole',mt_int,1,cell_expand)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 5100, RTDB_ERR)
      goto 10

c
c  balance
c
 5200 if (.not.rtdb_put(rtdb,'nwpw:balance',mt_log,1,.true.))
     >  call errquit('nwpw_input: error writing to rtdb',5200,RTDB_ERR)

      goto 10
c
c  nobalance
c
 5300 if (.not.rtdb_put(rtdb,'nwpw:balance',mt_log,1,.false.))
     >  call errquit('nwpw_input: error writing to rtdb',5300,RTDB_ERR)

      goto 10

c
c  bo_steps
c
 5400 if (.not. inp_i(loop(1))) loop(1) = 10
      if (.not. inp_i(loop(2))) loop(2) = 100
      if (.not.rtdb_put(rtdb,'nwpw:bo_steps',mt_int,2,loop))
     >  call errquit('nwpw_input: error writing to rtdb',5400,RTDB_ERR)

      goto 10

c
c  bo_time_step
c
 5500 if (.not. inp_f(wcut)) wcut = 15.0d0
      if (.not.rtdb_put(rtdb,'nwpw:bo_time_step',mt_dbl,1,wcut))
     >  call errquit('nwpw_input: error writing to rtdb',5500,RTDB_ERR)

      goto 10

c
c  bo_algorithm
c
 5600 if (.not.inp_a(zone_name)) nx = 0
      if (inp_compare(.false.,zone_name,'verlet'))          nx = 0
      if (inp_compare(.false.,zone_name,'velocity-verlet')) nx = 1
      if (inp_compare(.false.,zone_name,'leap-frog'))       nx = 2
      if (.not.rtdb_put(rtdb,'nwpw:bo_algorithm',mt_int,1,nx))
     >  call errquit('nwpw_input: error writing to rtdb',5600,RTDB_ERR)

      goto 10

c
c  bo_time_step
c
 5650 if (.not. inp_f(wcut)) wcut = 500.d0
      if (.not.rtdb_put(rtdb,'nwpw:bo_fake_mass',mt_dbl,1,wcut))
     >  call errquit('nwpw_input: error writing to rtdb',5650,RTDB_ERR)

      goto 10

c
c  scaling
c
 5700 if (.not. inp_f(scaling(1))) scaling(1) = 1.0d0
      if (.not. inp_f(scaling(2))) scaling(2) = scaling(1)
      if (.not.(rtdb_put(rtdb,'nwpw:scaling',mt_dbl,2,scaling).and.
     >          rtdb_put(rtdb,'cpmd:scaling',mt_dbl,2,scaling)))
     >  call errquit('nwpw_input: error writing to rtdb',5700,RTDB_ERR)

      goto 10

c
c  xyz_filename
c
 5800 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:xyz_filename',1,cell_name(1:ind))
     >   .and.rtdb_cput(rtdb,'nwpw:xyz_filename',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5800, RTDB_ERR)

      goto 10

c
c  np_dimensions (fft,orbital,k-space) 
c     - only values greater than 0 are meaningfull, less than zero are default
c
 5900 if (.not. inp_i(cell_expand(1))) cell_expand(1) = -1
      if (.not. inp_i(cell_expand(2))) cell_expand(2) = -1
      if (.not. inp_i(cell_expand(3))) cell_expand(3) = -1
      if (.not.rtdb_put(rtdb,'nwpw:np_dimensions',mt_int,3,cell_expand))
     >  call errquit('nwpw_input: error writing to rtdb',5900,RTDB_ERR)
      goto 10

c
c spin orbit option "watch out"
c
 6000 if (.not.inp_a(zone_name)) move = .true.
      if (inp_compare(.false.,zone_name,'on'))     move = .true.
      if (inp_compare(.false.,zone_name,'yes'))    move = .true.
      if (inp_compare(.false.,zone_name,'.true.')) move = .true.
      if (inp_compare(.false.,zone_name,'T'))      move = .true.
      if (inp_compare(.false.,zone_name,'off'))     move = .false.
      if (inp_compare(.false.,zone_name,'no'))      move = .false.
      if (inp_compare(.false.,zone_name,'.false.')) move = .false.
      if (inp_compare(.false.,zone_name,'F'))       move = .false.
      if (inp_compare(.false.,zone_name,'fix'))     move = .false.

      if (.not.rtdb_put(rtdb,'nwpw:spin_orbit',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',6000,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'pspw:spin_orbit',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',6000,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'band:spin_orbit',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',6000,RTDB_ERR)
         
      goto 10

c
c virtual_orbitals 
c
 6090 if (.not. inp_i(nx)) nx = 3
      if (.not. inp_i(ny)) ny = nx
      loop(1) = nx
      loop(2) = ny
      if (.not.rtdb_put(rtdb,'nwpw:excited_ne',mt_int,2,loop))
     >  call errquit('nwpw_input: error writing to rtdb',6090,RTDB_ERR)
      goto 10

c
c fractional_orbitals 
c
 6100 if (.not. inp_i(nx)) nx = 1
      if (.not. inp_i(ny)) ny = nx
      loop(1) = nx
      loop(2) = ny
      if (.not.rtdb_put(rtdb,'nwpw:fractional_orbitals',mt_int,2,loop))
     >  call errquit('nwpw_input: error writing to rtdb',6100,RTDB_ERR)
      goto 10

c
c smear temperature
c
 6200 alpha = 0.001d0/kb !* use same default as gaussian dft
      nz      = 2  !* smeartype  !* use same default as gaussian dft
      loop(1) = 4  !* frac_orb(1)
      loop(2) = 4  !* frac_orb(2)
      if (inp_f(ecut)) then
         alpha = ecut/kb
      else
         call inp_prev_field()
      end if
 6201 if (.not.inp_a(zone_name)) goto 6202
      if (inp_compare(.false.,zone_name,'temperature')) then
         if (.not.inp_f(alpha)) 
     >   call errquit(
     >    'nwpw_input:error reading smear temperature',6200,INPUT_ERR)
      end if
      if (inp_compare(.false.,zone_name,'fermi'))    nz = 1
      if (inp_compare(.false.,zone_name,'gaussian')) nz = 2
     
      if (inp_compare(.false.,zone_name,'orbitals')) then
         if (.not.inp_i(loop(1))) then
            loop(1) = 4
            call inp_prev_field()
         end if
         if (.not.inp_i(loop(2))) then
            loop(2) = loop(1)
            call inp_prev_field()
         end if
      end if
      go to 6201
 6202 if (.not.
     > (rtdb_put(rtdb,'nwpw:fractional_orbitals',mt_int,2,loop)    .and.
     >  rtdb_put(rtdb,'nwpw:fractional_temperature',mt_dbl,1,alpha).and.
     >  rtdb_put(rtdb,'nwpw:fractional_smeartype',mt_int,1,nz)))
     >  call errquit('nwpw_input: error writing to rtdb',6200,RTDB_ERR)

      go to 10

c
c lcao_skip
c
 6300 if (.not.rtdb_put(rtdb,'nwpw:lcao_skip',mt_log,1,.true.))
     >  call errquit('nwpw_input: error writing to rtdb',6300,RTDB_ERR)
      goto 10
c
c lcao
c
 6400 if (.not.rtdb_put(rtdb,'nwpw:lcao_skip',mt_log,1,.false.))
     >  call errquit('nwpw_input: error writing to rtdb',6400,RTDB_ERR)
      goto 10
c
c ewald_ngrid
c
 6500 if (inp_i(cell_expand(1))) then
        if (.not. inp_i(cell_expand(2))) cell_expand(2) = cell_expand(1)
        if (.not. inp_i(cell_expand(3))) cell_expand(3) = cell_expand(2)
        if (.not.rtdb_put(rtdb,'nwpw:ewald_ngrid',mt_int,3,cell_expand))
     >  call errquit('nwpw_input: error writing to rtdb',6500,RTDB_ERR)
      end if
      goto 10

c
c Nose-Hoover
c
 6600 nose   = .true.
      nosers = .true.
      Pe = 100.00d0
      Pr = 100.00d0
      Te = 298.15d0
      Tr = 298.15d0
      if (.not. inp_f(fe)) goto 6601
      Pe = fe
      if (.not. inp_f(fe)) goto 6601
      Te = fe
      if (.not. inp_f(fe)) goto 6601
      Pr = fe
      if (.not. inp_f(fe)) goto 6601
      Tr = fe
      if (.not.inp_a(zone_name)) goto 6601
      if (inp_compare(.false.,zone_name,'restart'))    nosers=.true.
      if (inp_compare(.false.,zone_name,'start'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'clear'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'initialize')) nosers=.false.
 6601 value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
     >   .and.rtdb_put(rtdb,'nwpw:nose',mt_log,1,nose)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:nose_restart',mt_log,1,nosers)
     >   .and.rtdb_put(rtdb,'nwpw:nose_restart',mt_log,1,nosers)
      value = value.and.rtdb_put(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value.and.rtdb_put(rtdb,'nwpw:Pe',mt_dbl,1,Pe)
      value = value.and.rtdb_put(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value.and.rtdb_put(rtdb,'nwpw:Te',mt_dbl,1,Te)
      value = value.and.rtdb_put(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value.and.rtdb_put(rtdb,'nwpw:Pr',mt_dbl,1,Pr)
      value = value.and.rtdb_put(rtdb,'cpmd:Tr',mt_dbl,1,Tr)
      value = value.and.rtdb_put(rtdb,'nwpw:Tr',mt_dbl,1,Tr)
      if (.not.value)
     > call errquit('nwpw_input: writing ', 6600, RTDB_ERR)
      goto 10


************************************************************************
c
c END/Finished
c
 9999 continue
      return
      end
