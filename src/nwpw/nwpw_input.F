c
c $Id$
c
**** checking git hub ***

*     ***********************************
*     *                                 *
*     *           nwpw_input            *
*     *                                 *
*     ***********************************
      subroutine nwpw_input(rtdb)
      implicit none
      integer rtdb

c 
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "errquit.fh"
#include "util.fh"
c#include "nwpwxc.fh"
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 130)
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
     >            'mc_steps',
     >            'bo_time_step',
     >            'bo_algorithm',
     >            'bo_fake_mass',
     >            'scaling',
     >            'xyz_filename',
     >            'ion_motion_filename',
     >            'emotion_filename',
     >            'hmotion_filename',
     >            'omotion_filename',
     >            'eigmotion_filename',
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
     >            'metadynamics',
     >            'pspspin',
     >            'Fei',
     >            'parallel_io',
     >            'efield_grad',
     >            'uterm',
     >            'et',
     >            'cpmd_properties',
     >            'use_grid_cmp',
     >            'lcao_mask',
     >            'tamd',
     >            'fast_erf',
     >            'fmm',
     >            'smooth_cutoff',
     >            'md_xs',
     >            'dipole_motion',
     >            'energy',
     >            'initial_velocities',
     >            'temperature',
     >            'cutoff_boot_wavefunction',
     >            'periodic_dipole',
     >            'translate_vector',
     >            'makehmass2',
     >            'apc',
     >            'socket',
     >            'vfield',
     >            'born',
     >            'symmetry',
     >            'one_electron_guess',
     >            'zone_fft_name',
     >            'dos-fft-grid',
     >            'end'/

      double precision kb
      parameter (kb=3.16679d-6)
      logical      value,move,mulliken_kawai,nose,nosers,not_m
      character*50 cell_name,zone_name,rtdb_name
      character*50 input_wavefunction_filename
      character*50 output_wavefunction_filename
      character*50 exchange_correlation
      character*128 ipname
      character*255 vnames
      double precision fake_mass,time_step,rcut
      integer loop(2),ncut,cell_expand(3)
      double precision tolerances(3),ecut,wcut,fe,ff,scaling(2)
      double precision alpha,emin,emax,agamma(10),Gc
      double precision Pe,Pr,Te,Tr
      integer npoints,mchain,nchain,nga
      logical mulliken,allow_translation
      integer ispin,multiplicity,space_avail
      integer nx,ny,nz,minimizer,h_values,i_values,nuterms
      integer index_start1(2),n1,j,jstart,jlast,jstride

chvd
      double precision a,b,c
      logical hfsr
      integer flag
      logical lstat,out1
      double precision xfac
chvd
       
*     **** external functions ****
      character*7 c_index_name
      external    c_index_name
      logical  nwpw_parse_boolean
      external nwpw_parse_boolean
      logical  nwpwxc_input,nwpwxc_rtdb_store,nwpwxc_has_hfx 
      external nwpwxc_input,nwpwxc_rtdb_store,nwpwxc_has_hfx 
      logical  nwpwxc_has_cam,nwpwxc_rtdb_delete
      external nwpwxc_has_cam,nwpwxc_rtdb_delete
      real*8   nwpwxc_wght_hfx
      external nwpwxc_wght_hfx

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
     >       5200,5300,5400,5400,5500,5600,5650,5700,
     >       5800,5810,5820,5830,5840,5850,
     >       5900,6000,6000,
     >       2750,6090,6100,1809,6200,6300,6400,6500,1510,6600,1807,
     >       6700,6800,6900,7000,4250,7100,7200,7300,7400,7500,7600,
     >       7700,7800,7900,8000,8100,6650,8200,8300,8400,8500,8600,
     >       8700,8800,8900,9000,9010,9020,9030,1810,1811,
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
c auxiliary_potentials
c
 1510 if (.not.rtdb_put(rtdb,'pspw_qmmm_auxon',mt_log,1,.true.))
     >   call errquit(
     >   'auxiliary_potentials:failed writing auxon:',0,
     >    RTDB_ERR)
      call pspw_qmmm_input(rtdb) 
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
      value = rtdb_put(rtdb,'band:tetrahedron-grid',mt_int,3,
     >                 cell_expand)
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1808, RTDB_ERR)
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

c
c  zone_fft_name
c
 1810 if (.not. inp_a(zone_name))
     >     call errquit(
     >          'nwpw_input: failed to read zone_fft_name',0,INPUT_ERR)

      ind = index(zone_name,' ') - 1
      value = rtdb_cput(rtdb,'band_fft:zone_name',1,
     >                  zone_name(1:ind))
      if (.not.value)
     >  call errquit(
     >       'nwpw_input: error writing to rtdb', 1810, RTDB_ERR)

      goto 10
c
c dos-fft-grid
c
 1811 if (.not. inp_i(cell_expand(1))) cell_expand(1) = 2
      if (.not. inp_i(cell_expand(2))) cell_expand(2) = 2
      if (.not. inp_i(cell_expand(3))) cell_expand(3) = 2
      if (.not. inp_a(zone_name)) zone_name = 'zone_fft_default'

      call dos_states_set(rtdb,zone_name,cell_expand(1),
     >                                   cell_expand(2),
     >                                   cell_expand(3))

      goto 10



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
      value = rtdb_cput(rtdb,'cgsd:output_wavefunction_filename',
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
      if (.not. inp_f(fe)) fe=1.0e-4
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
 2600 if (.not.inp_f(ecut))
     > call errquit('nwpw_input: failed to read ecut',0,INPUT_ERR)
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
     > call errquit('nwpw_input: failed to read wcut',0,INPUT_ERR)
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
 2750 if (.not.inp_f(wcut))
     > call errquit('nwpw_input: failed to read wcut',0,INPUT_ERR)
      if (.not.inp_f(ecut)) ecut = 2.0d0*wcut

      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
     >   .and.rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
     >   .and.rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
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
      if (exchange_correlation(1:ind).eq.'new') then
c       call inp_prev_field()
        if (.not.nwpwxc_input()) then
        !call nwpwxc_input(lstat)
        !if(.not.lstat) then
          call errquit("nwpw_input: not a functional",
     +                 0,UERR)
        endif
        if (.not.nwpwxc_rtdb_store(rtdb,"dft")) then
          call errquit("nwpw_input: store functional failed",
     +                 0,UERR)
        endif
        if (nwpwxc_has_hfx()) then
          xfac = nwpwxc_wght_hfx()
          value = rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
          value = value .and.
     +            rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)
          value = value .and.
     +            rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,
     +                     xfac)
          value = value .and.
     +            rtdb_put(rtdb,'band:HFX_parameter',mt_dbl,1,
     +                     xfac)
          if (.not.value)
     +    call errquit("nwpw_input: store HFX A failed",0,UERR)
        else
          value = rtdb_put(rtdb,'pspw:HFX',mt_log,1,.false.)
          value = value .and.
     +            rtdb_put(rtdb,'band:HFX',mt_log,1,.false.)
          if (.not.value)
     +    call errquit("nwpw_input: store HFX B failed",0,UERR)
        endif
        if (nwpwxc_has_cam()) then
          call nwpwxc_get_cam(a,b,c,hfsr)
          if (hfsr) then
            flag = 4
          else
            flag = 2
          endif
          value = rtdb_put(rtdb,'pspw:HFX_screening_type',mt_int,1,flag)
          value = value .and.
     +            rtdb_put(rtdb,'band:HFX_screening_type',mt_int,1,flag)
          value = value .and.
     +            rtdb_put(rtdb,'pspw:HFX_solver_type',mt_int,1,1)
          value = value .and.
     +            rtdb_put(rtdb,'band:HFX_solver_type',mt_int,1,1)
          value = value .and. 
     +            rtdb_put(rtdb,'nwpw:attenuation',mt_dbl,1,c)
          if (.not.value)
     +    call errquit("nwpw_input: store CAMHFX failed",0,UERR)
        endif
        goto 10
      endif

      value = nwpwxc_rtdb_delete(rtdb,"dft") !*** don't need to check return value ****

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
      value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.false.)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'band:HFX_parameter',mt_dbl,1,1.0d0)

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
c  DOS - density of states input, alpha,npoints, emin,emax
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
c  EFIELD_GRAD
c
 4250 if (.not.rtdb_put(rtdb,'nwpw:efield_grad',mt_log,1,.true.))
     >  call errquit('nwpw_input: error writing to rtdb',4250,RTDB_ERR)

      goto 10

c
c occupations
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
       if (inp_compare(.false.,zone_name,'block-cg')) then
         if(.not.rtdb_put(rtdb,'nwpw:ks_algorithm',mt_int,1,-1))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
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
       if (inp_compare(.false.,zone_name,'diis').or.
     >     inp_compare(.false.,zone_name,'pulay').or.
     >     inp_compare(.false.,zone_name,'johnson-pulay')
     >    ) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,2))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'anderson')) then
         if(.not.rtdb_put(rtdb,'nwpw:scf_algorithm',mt_int,1,3))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'noprecondition')) then
         if(.not.rtdb_put(rtdb,'nwpw:precondition',mt_log,1,.false.))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'precondition')) then
         if(.not.rtdb_put(rtdb,'nwpw:precondition',mt_log,1,.true.))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if

       if (inp_compare(.false.,zone_name,'kerker')) then
         if (.not.inp_f(alpha)) 
     >   call errquit('nwpw_input: error reading kerker g0',
     >                4800,RTDB_ERR)
         if(.not.rtdb_put(rtdb,'nwpw:kerker_g0',mt_dbl,1,alpha))
     >   call errquit('nwpw_input: error writing to rtdb',4800,RTDB_ERR)
       end if
       if (inp_compare(.false.,zone_name,'alpha')) then
         if (.not.inp_f(alpha)) 
     >   call errquit('nwpw_input: error reading scf alpha',
     >                4800,RTDB_ERR)
         if(.not.rtdb_put(rtdb,'nwpw:ks_alpha',mt_dbl,1,alpha))
     >   call errquit('nwpw_input: error writing to rtdb',
     >                4800,RTDB_ERR)
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
      move = nwpw_parse_boolean(zone_name,.false.)
c       if (inp_compare(.false.,zone_name,'on'))     move = .true.
c       if (inp_compare(.false.,zone_name,'yes'))    move = .true.
c       if (inp_compare(.false.,zone_name,'.true.')) move = .true.
c       if (inp_compare(.false.,zone_name,'T'))      move = .true.
c       if (inp_compare(.false.,zone_name,'off'))     move = .false.
c       if (inp_compare(.false.,zone_name,'no'))      move = .false.
c       if (inp_compare(.false.,zone_name,'.false.')) move = .false.
c       if (inp_compare(.false.,zone_name,'F'))       move = .false.
c       if (inp_compare(.false.,zone_name,'fix'))     move = .false.
      if (.not.rtdb_put(rtdb,'nwpw:rotation',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',4900,RTDB_ERR)

 4901 continue
      goto 10


c
c  translation
c
 5000 if (.not.inp_a(zone_name)) goto 5001
      move = nwpw_parse_boolean(zone_name,.true.)
c      if (inp_compare(.false.,zone_name,'on'))     move = .true.
c      if (inp_compare(.false.,zone_name,'yes'))    move = .true.
c      if (inp_compare(.false.,zone_name,'.true.')) move = .true.
c      if (inp_compare(.false.,zone_name,'T'))      move = .true.
c      if (inp_compare(.false.,zone_name,'off'))     move = .false.
c      if (inp_compare(.false.,zone_name,'no'))      move = .false.
c      if (inp_compare(.false.,zone_name,'.false.')) move = .false.
c      if (inp_compare(.false.,zone_name,'F'))       move = .false.
c      if (inp_compare(.false.,zone_name,'fix'))     move = .false.
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

      !*** try reading in atom indexes ****
      if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'nwpw_input:failed allocating index_start1',0,MA_ERR)
      n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do

      if (.not.rtdb_put(rtdb,'nwpw:scaling_natms',
     >                  mt_int,1,n1))
     >   call errquit('nwpw_input:nwpw:scaling_natms',0,RTDB_ERR)

      if (n1.gt.0) then
         if (.not.rtdb_put(rtdb,'nwpw:scaling_atoms',
     >                     mt_int,n1,int_mb(index_start1(1))))
     >      call errquit('nwpw_input:nwpw:scaling_atoms',0,RTDB_ERR)
      end if
      if (.not.BA_pop_stack(index_start1(2))) call errquit(
     > 'nwpw_input:failed deallocating index_start1',0,MA_ERR)

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
c  ion_motion_filename
c
 5810 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:ion_motion_filename',1,
     >                  cell_name(1:ind))
     >  .and. rtdb_cput(rtdb,'nwpw:ion_motion_filename',1,
     >                  cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5810, RTDB_ERR)
      goto 10
c
c  emotion_filename
c
 5820 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:emotion_filename',1,cell_name(1:ind))
     >  .and. rtdb_cput(rtdb,'nwpw:emotion_filename',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5820, RTDB_ERR)
      goto 10
c
c  hmotion_filename
c
 5830 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:hmotion_filename',1,cell_name(1:ind))
     >  .and. rtdb_cput(rtdb,'nwpw:hmotion_filename',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5830, RTDB_ERR)
      goto 10
c
c  omotion_filename
c
 5840 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:omotion_filename',1,cell_name(1:ind))
     >  .and. rtdb_cput(rtdb,'nwpw:omotion_filename',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5840, RTDB_ERR)
      goto 10
c
c  eigmotion_filename
c
 5850 if (.not. inp_a(cell_name)) goto 10
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:eigmotion_filename',1,
     >                  cell_name(1:ind))
     >  .and. rtdb_cput(rtdb,'nwpw:eigmotion_filename',1,
     >                  cell_name(1:ind))
      if (.not.value)
     > call errquit('nwpw_input: writing ', 5850, RTDB_ERR)
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
      move = nwpw_parse_boolean(zone_name,move)
c      if (inp_compare(.false.,zone_name,'on'))     move = .true.
c      if (inp_compare(.false.,zone_name,'yes'))    move = .true.
c      if (inp_compare(.false.,zone_name,'.true.')) move = .true.
c      if (inp_compare(.false.,zone_name,'T'))      move = .true.
c      if (inp_compare(.false.,zone_name,'off'))     move = .false.
c      if (inp_compare(.false.,zone_name,'no'))      move = .false.
c      if (inp_compare(.false.,zone_name,'.false.')) move = .false.
c      if (inp_compare(.false.,zone_name,'F'))       move = .false.
c      if (inp_compare(.false.,zone_name,'fix'))     move = .false.

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
 6200 ff = 0.001d0/kb !* use same default as gaussian dft
      alpha = 1.0d0 !* smear mixing parameter
      nz      = 2  !* smeartype  !* use same default as gaussian dft
      if (.not.rtdb_get(rtdb,'nwpw:fractional_orbitals',
     >                  mt_int,2,loop)) then
         loop(1) = 4  !* frac_orb(1)
         loop(2) = 4  !* frac_orb(2)
      end if
      if (.not.rtdb_get(rtdb,'nwpw:fractional_temperature',
     >                  mt_dbl,1,ff)) then
         ff = 0.001d0/kb !* use same default as gaussian dft
      end if
      if (.not.rtdb_get(rtdb,'nwpw:fractional_alpha',
     >                  mt_dbl,1,alpha)) then
         alpha = 1.0d0 !* smear mixing parameter
      end if
      if (.not.rtdb_get(rtdb,'nwpw:fractional_smeartype',
     >                  mt_int,1,nz)) then
         nz      = 2  !* smeartype  !* use same default as gaussian dft
      end if

      if (inp_f(ecut)) then
         ff = ecut/kb
      else
         call inp_prev_field()
      end if
 6201 if (.not.inp_a(zone_name)) goto 6202
      if (inp_compare(.false.,zone_name,'temperature')) then
         if (.not.inp_f(ff)) 
     >   call errquit(
     >    'nwpw_input:error reading smear temperature',6200,INPUT_ERR)
      end if
      if (inp_compare(.false.,zone_name,'alpha')) then
         if (.not.inp_f(alpha)) 
     >   call errquit(
     >    'nwpw_input:error reading smear alpha',6200,INPUT_ERR)
      end if
      if (inp_compare(.false.,zone_name,'fixed'))    nz = -1
      if (inp_compare(.false.,zone_name,'step'))     nz = 0
      if (inp_compare(.false.,zone_name,'fermi'))    nz = 1
      if (inp_compare(.false.,zone_name,'gaussian')) nz = 2
      if (inp_compare(.false.,zone_name,'marzari-vanderbilt')) nz = 4

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
     > (rtdb_put(rtdb,'nwpw:fractional_orbitals',mt_int,2,loop) .and.
     >  rtdb_put(rtdb,'nwpw:fractional_temperature',mt_dbl,1,ff).and.
     >  rtdb_put(rtdb,'nwpw:fractional_alpha',mt_dbl,1,alpha).and.
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
      mchain = 1
      nchain = 1
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
      if (.not.inp_i(mchain)) goto 6601
      if (.not.inp_i(nchain)) goto 6601
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
      value = value.and.rtdb_put(rtdb,'cpmd:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'nwpw:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'cpmd:Nchain',mt_int,1,nchain)
      value = value.and.rtdb_put(rtdb,'nwpw:Nchain',mt_int,1,nchain)
      if (.not.value)
     > call errquit('nwpw_input: writing ', 6600, RTDB_ERR)
      goto 10

c
c energy
c
 6650 nose = .false.
      value =    rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
     >      .and.rtdb_put(rtdb,'nwpw:nose',mt_log,1,nose)
      if (.not.value)
     > call errquit('nwpw_input: writing ', 6650, RTDB_ERR)

      goto 10

c
c  metadynamics
c
 6700 call metadynamics_input(rtdb)
      goto 10


c
c  pspspin up d scale atom_indexes
c
 6800 if (inp_a(zone_name)) then
         value = .true.

         if (inp_compare(.false.,zone_name,'off')) then
            value = rtdb_put(rtdb,'nwpw:pspspin',mt_log,1,.false.)
            value = rtdb_put(rtdb,'nwpw:pspspin_count',mt_int,1,0)
            goto 6810
         end if
         if (inp_compare(.false.,zone_name,'.false.')) then
            value = rtdb_put(rtdb,'nwpw:pspspin',mt_log,1,.false.)
            value = rtdb_put(rtdb,'nwpw:pspspin_count',mt_int,1,0)
            goto 6810
         end if
         if (inp_compare(.false.,zone_name,'F')) then
            value = rtdb_put(rtdb,'nwpw:pspspin',mt_log,1,.false.)
            value = rtdb_put(rtdb,'nwpw:pspspin_count',mt_int,1,0)
            goto 6810
         end if

         if (.not.rtdb_get(rtdb,'nwpw:pspspin_count',mt_int,1,nz)) nz=0
         nz = nz + 1
         

*        **** add m options ****
         not_m = .false.
         if (inp_compare(.false.,zone_name,'not_m')) then
            if (inp_i(ny)) not_m = .true.
            if (.not.inp_a(zone_name)) zone_name = "?????"
         end if

         if (inp_compare(.false.,zone_name,'up') .or.
     >       inp_compare(.false.,zone_name,'alpha')) then
            value = rtdb_put(rtdb,'nwpw:pspspin',mt_log,1,.true.)
            rtdb_name = 'nwpw:pspspin_iamup:'//c_index_name(nz)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_log,1,.true.)
            if (inp_a(cell_name).and.inp_f(emax)) then
               rtdb_name = 'nwpw:pspspin_upscale:'//c_index_name(nz)
               value = value.and.
     >         rtdb_put(rtdb,rtdb_name,mt_dbl,1,emax)
               nx = -1
               if (inp_compare(.false.,cell_name,'s')) nx = 0
               if (inp_compare(.false.,cell_name,'p')) nx = 1
               if (inp_compare(.false.,cell_name,'d')) nx = 2
               if (inp_compare(.false.,cell_name,'f')) nx = 3
               if (inp_compare(.false.,cell_name,'0')) nx = 0
               if (inp_compare(.false.,cell_name,'1')) nx = 1
               if (inp_compare(.false.,cell_name,'2')) nx = 2
               if (inp_compare(.false.,cell_name,'3')) nx = 3
               rtdb_name = 'nwpw:pspspin_upl:'//c_index_name(nz)
               value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nx)
               if (not_m) then
                  rtdb_name = 'nwpw:pspspin_upm:'//c_index_name(nz)
                  value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,ny)
               end if

               space_avail = min(100000,ma_inquire_avail(mt_int))
               if (.not.ma_push_get(mt_int,space_avail,'input scratch',
     >            h_values, i_values) ) call errquit(
     >            'nwpw_input: unable to allocate input scratch space',
     >            space_avail, MA_ERR)
               if (.not.inp_ilist(space_avail,int_mb(i_values),nx))
     >         call errquit('nwpw_input: failed reading integer list',
     >                      0, INPUT_ERR)
               rtdb_name = 'nwpw:pspspin_upions:'//c_index_name(nz)
               value = value.and.
     >         rtdb_put(rtdb,rtdb_name,mt_int,nx,int_mb(i_values))
               if (.not.ma_pop_stack(h_values)) call errquit(
     >         'nwpw_input: unable to deallocate input scratch space',
     >                      0, INPUT_ERR)

            end if
         end if

         if (inp_compare(.false.,zone_name,'down') .or.
     >       inp_compare(.false.,zone_name,'beta')) then
            value = rtdb_put(rtdb,'nwpw:pspspin',mt_log,1,.true.)
            rtdb_name = 'nwpw:pspspin_iamup:'//c_index_name(nz)
            value = value.and.rtdb_put(rtdb,rtdb_name,mt_log,1,.false.)
            if (inp_a(cell_name).and.inp_f(emax)) then
               rtdb_name = 'nwpw:pspspin_downscale:'//c_index_name(nz)
               value = value.and.rtdb_put(rtdb,rtdb_name,mt_dbl,1,emax)
               nx = -1
               if (inp_compare(.false.,cell_name,'s')) nx = 0
               if (inp_compare(.false.,cell_name,'p')) nx = 1
               if (inp_compare(.false.,cell_name,'d')) nx = 2
               if (inp_compare(.false.,cell_name,'f')) nx = 3
               if (inp_compare(.false.,cell_name,'0')) nx = 0
               if (inp_compare(.false.,cell_name,'1')) nx = 1
               if (inp_compare(.false.,cell_name,'2')) nx = 2
               if (inp_compare(.false.,cell_name,'3')) nx = 3
               rtdb_name = 'nwpw:pspspin_downl:'//c_index_name(nz)
               value = value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,nx)
               if (not_m) then
                  rtdb_name = 'nwpw:pspspin_downm:'//c_index_name(nz)
                  value=value.and.rtdb_put(rtdb,rtdb_name,mt_int,1,ny)
               end if

               space_avail = min(100000,ma_inquire_avail(mt_int))
               if (.not.ma_push_get(mt_int,space_avail,'input scratch',
     >            h_values, i_values) ) call errquit(
     >            'nwpw_input: unable to allocate input scratch space',
     >            space_avail, MA_ERR)
               if (.not.inp_ilist(space_avail,int_mb(i_values),nx))
     >         call errquit('nwpw_input: failed reading integer list', 
     >                      0, INPUT_ERR)
               rtdb_name = 'nwpw:pspspin_downions:'//c_index_name(nz)
               value = value.and.
     >         rtdb_put(rtdb,rtdb_name,mt_int,nx,int_mb(i_values))
               if (.not.ma_pop_stack(h_values)) call errquit(
     >         'nwpw_input: unable to deallocate input scratch space',
     >                      0, INPUT_ERR)

            end if
         end if
         value = value.and.
     >           rtdb_put(rtdb,'nwpw:pspspin_count',mt_int,1,nz)

 6810    if (.not.value)
     >    call errquit('nwpw_input: writing ', 6800, RTDB_ERR)

      end if

      goto 10
c
c Fei
c
 6900 if (inp_a(input_wavefunction_filename)) then
         ind = index(input_wavefunction_filename,' ') - 1
         value = rtdb_cput(rtdb,'cpmd:fei_filename',
     >                    1,input_wavefunction_filename(1:ind))
     >      .and.rtdb_cput(rtdb,'nwpw:fei_filename',
     >                    1,input_wavefunction_filename(1:ind))
          if (.not.value)
     >    call errquit('nwpw_input: writing ', 6900, RTDB_ERR)
      end if
      value = .true.
      if (.not.rtdb_put(rtdb,'cpmd:fei',mt_log,1,value))
     > call errquit('nwpw_input: writing fei ', 6900, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'nwpw:fei',mt_log,1,value))
     > call errquit('nwpw_input: writing fei ', 6900, RTDB_ERR)

      goto 10
c
c  parallel_io
c
 7000 if (inp_a(zone_name)) then
         move = .true.
         move = nwpw_parse_boolean(zone_name,move)
c         if (inp_compare(.false.,zone_name,'on'))     move = .true.
c         if (inp_compare(.false.,zone_name,'yes'))    move = .true.
c         if (inp_compare(.false.,zone_name,'.true.')) move = .true.
c         if (inp_compare(.false.,zone_name,'T'))      move = .true.
c         if (inp_compare(.false.,zone_name,'off'))     move = .false.
c         if (inp_compare(.false.,zone_name,'no'))      move = .false.
c         if (inp_compare(.false.,zone_name,'.false.')) move = .false.
c         if (inp_compare(.false.,zone_name,'F'))       move = .false.
c         if (inp_compare(.false.,zone_name,'fix'))     move = .false.
      else
         move = .true.
      end if
      if (.not.rtdb_put(rtdb,'nwpw:parallel_io',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',7000,RTDB_ERR)

      goto 10
c
c  uterm l scale atom_indexes
c  uterm s 1.2 2:3
c  uterm d 1.4 4:20:2
c  uterm off
c
 7100 if (inp_a(zone_name)) then
         value = .true.
         if (.not.rtdb_get(rtdb,'nwpw:nuterms',mt_int,1,nuterms)) then
            nuterms = 0
         end if
         nuterms = nuterms + 1
         
         if ((inp_compare(.false.,zone_name,'off')).or.
     >       (inp_compare(.false.,zone_name,'.false.'))) then
            value = rtdb_put(rtdb,'nwpw:uterm',mt_log,1,.false.)
            value = value.and.
     >              rtdb_put(rtdb,'nwpw:nuterms',mt_int,1,0)
            goto 7110
         end if
         value = rtdb_put(rtdb,'nwpw:uterm',mt_log,1,.true.)
         value = value.and.
     >           rtdb_put(rtdb,'nwpw:nuterms',mt_int,1,nuterms)


         if (inp_f(emax).and.inp_f(emin)) then
               cell_name = 'nwpw:uterm_scale:'//c_index_name(nuterms)
               value = value.and.
     >         rtdb_put(rtdb,cell_name,mt_dbl,1,emax)
               cell_name = 'nwpw:jterm_scale:'//c_index_name(nuterms)
               value = value.and.
     >         rtdb_put(rtdb,cell_name,mt_dbl,1,emin)

               nx = -1
               if (inp_compare(.false.,zone_name,'s')) nx = 0
               if (inp_compare(.false.,zone_name,'p')) nx = 1
               if (inp_compare(.false.,zone_name,'d')) nx = 2
               if (inp_compare(.false.,zone_name,'f')) nx = 3
               if (inp_compare(.false.,zone_name,'0')) nx = 0
               if (inp_compare(.false.,zone_name,'1')) nx = 1
               if (inp_compare(.false.,zone_name,'2')) nx = 2
               if (inp_compare(.false.,zone_name,'3')) nx = 3
               cell_name = 'nwpw:uterm_l:'//c_index_name(nuterms)
               value = value.and.
     >         rtdb_put(rtdb,cell_name,mt_int,1,nx)

               space_avail = min(100000,ma_inquire_avail(mt_int))
               if (.not.ma_push_get(mt_int,space_avail,'input scratch',
     >            h_values, i_values) ) call errquit(
     >            'nwpw_input: unable to allocate input scratch space',
     >            space_avail, MA_ERR)
               if (.not.inp_ilist(space_avail,int_mb(i_values),ny))
     >         call errquit('nwpw_input: failed reading integer list',
     >                      0, INPUT_ERR)
               cell_name = 'nwpw:uterm_ions:'//c_index_name(nuterms)
               value = value.and.
     >         rtdb_put(rtdb,cell_name,mt_int,ny,int_mb(i_values))
               if (.not.ma_pop_stack(h_values)) call errquit(
     >         'nwpw_input: unable to deallocate input scratch space',
     >                      0, INPUT_ERR)
         end if

 7110    if (.not.value)
     >    call errquit('nwpw_input: writing ', 7100, RTDB_ERR)

      end if

      goto 10
c
c  et movecs_a movecs_b ion_a ion_b
c
 7200 if (inp_a(zone_name)) then
         ind = index(zone_name,' ') - 1
         value = rtdb_cput(rtdb,'pspw:et:movecs_a',1,zone_name(1:ind))
         if (.not.value)
     >    call errquit('nwpw_input: writing ',7200,RTDB_ERR)
      end if 
      if (inp_a(zone_name)) then
         ind = index(zone_name,' ') - 1
         value = rtdb_cput(rtdb,'pspw:et:movecs_b',1,zone_name(1:ind))
         if (.not.value)
     >    call errquit('nwpw_input: writing ',7200,RTDB_ERR)
      end if 
      if (inp_a(zone_name)) then
         ind = index(zone_name,' ') - 1
         value = rtdb_cput(rtdb,'pspw:et:ion_a',1,zone_name(1:ind))
         if (.not.value)
     >    call errquit('nwpw_input: writing ',7200,RTDB_ERR)
      end if 
      if (inp_a(zone_name)) then
         ind = index(zone_name,' ') - 1
         value = rtdb_cput(rtdb,'pspw:et:ion_b',1,zone_name(1:ind))
         if (.not.value)
     >    call errquit('nwpw_input: writing ',7200,RTDB_ERR)
      end if 

      goto 10

c
c cpmd_properties
c
 7300 if (.not.inp_a(zone_name)) move = .true.
      move = nwpw_parse_boolean(zone_name,move)
      if (.not.rtdb_put(rtdb,'nwpw:cpmd_properties',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',7300,RTDB_ERR)

      goto 10

c
c use_grid_cmp
c
 7400 if (.not.inp_a(zone_name)) move = .true.
      move = nwpw_parse_boolean(zone_name,move)
      if (.not.rtdb_put(rtdb,'nwpw:use_grid_cmp',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',7400,RTDB_ERR)

      goto 10


c
c  lcao_mask up aorb_indexes
c  lcao_mask down aorb_indexes
c
 7500 if (inp_a(zone_name)) then
         value = .true.

         if (inp_compare(.false.,zone_name,'off')) then
            value = rtdb_put(rtdb,'nwpw:lcao_mask',mt_log,1,.false.)
            goto 7510
         end if
         if (inp_compare(.false.,zone_name,'.false.')) then
            value = rtdb_put(rtdb,'nwpw:lcao_mask',mt_log,1,.false.)
            goto 7510
         end if
         if (inp_compare(.false.,zone_name,'F')) then
            value = rtdb_put(rtdb,'nwpw:lcao_mask',mt_log,1,.false.)
            goto 7510
         end if

         if (inp_compare(.false.,zone_name,'up') .or.
     >       inp_compare(.false.,zone_name,'alpha')) then
            value = rtdb_put(rtdb,'nwpw:lcao_mask',mt_log,1,.true.)

               space_avail = min(100000,ma_inquire_avail(mt_int))
               if (.not.ma_push_get(mt_int,space_avail,'input scratch',
     >            h_values, i_values) ) call errquit(
     >            'nwpw_input: unable to allocate input scratch space',
     >            space_avail, MA_ERR)
               if (.not.inp_ilist(space_avail,int_mb(i_values),nx))
     >         call errquit('nwpw_input: failed reading integer list',
     >                      0, INPUT_ERR)
               value = value.and.
     >         rtdb_put(rtdb,'nwpw:lcao_mask_uporbs',mt_int,
     >                  nx,int_mb(i_values))
               if (.not.ma_pop_stack(h_values)) call errquit(
     >         'nwpw_input: unable to deallocate input scratch space',
     >                      0, INPUT_ERR)

         end if
         if (inp_compare(.false.,zone_name,'down') .or.
     >       inp_compare(.false.,zone_name,'beta')) then
            value = rtdb_put(rtdb,'nwpw:lcao_mask',mt_log,1,.true.)

               space_avail = min(100000,ma_inquire_avail(mt_int))
               if (.not.ma_push_get(mt_int,space_avail,'input scratch',
     >            h_values, i_values) ) call errquit(
     >            'nwpw_input: unable to allocate input scratch space',
     >            space_avail, MA_ERR)
               if (.not.inp_ilist(space_avail,int_mb(i_values),nx))
     >         call errquit('nwpw_input: failed reading integer list',
     >                      0, INPUT_ERR)
               value = value.and.
     >         rtdb_put(rtdb,'nwpw:lcao_mask_downorbs',mt_int,
     >                  nx,int_mb(i_values))
               if (.not.ma_pop_stack(h_values)) call errquit(
     >         'nwpw_input: unable to deallocate input scratch space',
     >                      0, INPUT_ERR)
         end if

 7510    if (.not.value)
     >    call errquit('nwpw_input: writing ', 6800, RTDB_ERR)

      end if

      goto 10


c
c  TAMD
c
 7600 call tamd_input(rtdb)
      goto 10

c
c  fast_erf
c
 7700 if (inp_a(zone_name)) then
         move = .true.
         move = nwpw_parse_boolean(zone_name,move)
      else
         move = .true. 
      end if
      if (.not.rtdb_put(rtdb,'nwpw:fast_erf',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',7700,RTDB_ERR)

      goto 10

c
c  fmm [true/false... [lmax] [lr]]
c
 7800 mchain = 10
      nchain = 1
      if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,move)
         if (.not.inp_i(mchain)) mchain = 10
         if (.not.inp_i(nchain)) nchain = 1
      else
         move = .true.
      end if

 7810 if (.not.rtdb_put(rtdb,'nwpw:fmm',mt_log,1,move))
     >  call errquit('nwpw_input: error writing to rtdb',7800,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'nwpw:fmm_lmax',mt_int,1,mchain))
     >  call errquit('nwpw_input: error writing to rtdb',7800,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'nwpw:fmm_lr',mt_int,1,nchain))
     >  call errquit('nwpw_input: error writing to rtdb',7800,RTDB_ERR)

      goto 10

c
c  smooth_cutoff [true/false ... Afac sigma]
c
 7900 if (.not.inp_a(zone_name)) then
         move = .true.
         scaling(1) = 2.0d0
         scaling(2) = 4.0d0
      else 
         move = nwpw_parse_boolean(zone_name,.true.)
         if (.not.inp_f(scaling(1))) scaling(1) = 2.0d0
         if (.not.inp_f(scaling(2))) scaling(2) = 4.0d0
      end if
      if (move) then
         if (.not.rtdb_put(rtdb,'nwpw:smooth_cutoff',mt_dbl,2,scaling))
     >   call errquit('nwpw_input: error writing to rtdb',7900,RTDB_ERR)
      else
         if (.not.rtdb_delete(rtdb,'nwpw:smooth_cutoff'))
     >   call errquit('nwpw_input: error writing to rtdb',7900,RTDB_ERR)
      end if

      goto 10

c
c  md_xs 
c
 8000 call md_xs_input(rtdb)
      goto 10

c
c dipole_motion
c
 8100 if (inp_a(input_wavefunction_filename)) then
         ind = index(input_wavefunction_filename,' ') - 1
         value = rtdb_cput(rtdb,'nwpw:dipole_motion_filename',
     >                    1,input_wavefunction_filename(1:ind))
          if (.not.value)
     >    call errquit('nwpw_input: writing ', 8100, RTDB_ERR)
      end if
      value = .true.
      if (.not.rtdb_put(rtdb,'nwpw:dipole_motion',mt_log,1,value))
     > call errquit('nwpw_input: writing dipole',8100,RTDB_ERR)

      goto 10
c
c initial_velocities temperature seed
c
 8200 if (inp_f(scaling(1))) then            !*** in units of K ***
         if (.not.inp_i(ispin)) ispin = 494

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities_temperature',
     >                     mt_dbl,1,scaling))
     >   call errquit('cpmd_input: initial temperature',8200,RTDB_ERR)

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities_seed',
     >                     mt_int,1,ispin))
     >   call errquit('cpmd_input: initial temperature',8200,RTDB_ERR)

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities',
     >                  mt_log,1,.true.))
     >   call errquit('cpmd_input: initial velocities',8200,RTDB_ERR)

      else
         value = rtdb_delete(rtdb,'nwpw:init_velocities')
      end if
      goto 10

c
c temperature Tion Pion Telc Pion Pelc start/restart nchain mchain
c - a reordering of the parameters for nose-hoover
c
 8300 nose   = .true.
      nosers = .true.
      Tr = 298.15d0
      Pr = 1200.00d0
      Te = 298.15d0
      Pe = 1200.00d0
      mchain = 1
      nchain = 1
      if (.not. inp_f(fe)) goto 8301
      Tr = fe
      Te = fe
      if (.not. inp_f(fe)) goto 8301
      Pr = fe
      Pe = fe
      if (.not. inp_f(fe)) goto 8301
      Te = fe
      if (.not. inp_f(fe)) goto 8301
      Pe = fe
 8301 if (.not.inp_a(zone_name)) goto 8302
      if (inp_compare(.false.,zone_name,'restart'))    nosers=.true.
      if (inp_compare(.false.,zone_name,'start'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'clear'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'initialize')) nosers=.false.
 8302 if (.not.inp_i(nchain)) goto 8303
      mchain = nchain
      if (.not.inp_i(mchain)) goto 8303

 8303 value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
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
      value = value.and.rtdb_put(rtdb,'cpmd:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'nwpw:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'cpmd:Nchain',mt_int,1,nchain)
      value = value.and.rtdb_put(rtdb,'nwpw:Nchain',mt_int,1,nchain)
      if (.not.value)
     > call errquit('nwpw_input: writing ', 8300, RTDB_ERR)
      goto 10

c
c cutoff_boot_wavefunction .true.||.false
c
 8400 if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,.true.)
      else
         move = .true.
      end if
      if (.not.rtdb_put(rtdb,'nwpw:cutoff_boot_psi',mt_log,1,
     >                  move))
     >  call errquit('nwpw_input: error writing to rtdb',8400,RTDB_ERR)
      goto 10

c
c periodic_dipole  .true.||.false
c
 8500 if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,.true.)
      else
         move = .true.
      end if
      if (.not.rtdb_put(rtdb,'nwpw:periodic_dipole',mt_log,1,
     >                  move))
     >  call errquit('nwpw_input: error writing to rtdb',8500,RTDB_ERR)
      goto 10

c
c translate_vector x y z [geom_name [reorder]]
c
 8600 if (inp_f(tolerances(1))) then
         if (.not.inp_f(tolerances(2))) tolerances(2) = tolerances(1)
         if (.not.inp_f(tolerances(3))) tolerances(3) = tolerances(2)
         if (inp_a(zone_name)) then
            ind = index(zone_name,' ') - 1
            if (.not.rtdb_cput(rtdb,'nwpw:translate_geom_name',
     >                         1,zone_name(1:ind)))
     >      call errquit(
     >       'nwpw_input: error writing to rtdb', 8600, RTDB_ERR)
         end if
         if (inp_a(zone_name)) then
            move = .false.
            if (inp_compare(.false.,zone_name,'reorder'))   move=.true.
            if (inp_compare(.false.,zone_name,'noreorder')) move=.false.
            if (.not.rtdb_put(rtdb,'nwpw:translate_reorder',
     >                        mt_log,1,move))
     >      call errquit('nwpw_input: error writing to rtdb',
     >                    8600,RTDB_ERR)

         end if
         if (.not.rtdb_put(rtdb,'nwpw:translate_vector',
     >                     mt_dbl,3,tolerances))
     >   call errquit('nwpw_input: error writing to rtdb',8600,RTDB_ERR)
      end if
      goto 10


c
c makehmass2  .true.||.false
c
 8700 if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,.true.)
      else
         move = .true.
      end if
      if (.not.rtdb_put(rtdb,'nwpw:makehmass2',mt_log,1,
     >                  move))
     >  call errquit('nwpw_input: error writing to rtdb',8700,RTDB_ERR)
      goto 10

c
c apc Gc gamma1 gamma2 gamma3 ...
c
 8800 if (inp_f(Gc)) then
         nga = 0
         value = .false.
         do while (inp_f(agamma(nga+1)).and.(.not.value))
            nga = nga + 1
            value = (nga.ge.10)
         end do
         if (nga.gt.0) then
         if ((.not.rtdb_put(rtdb,'nwpw_APC:Gc',mt_dbl,1,Gc)) .or.
     >    (.not.rtdb_put(rtdb,'nwpw_APC:nga',mt_int,1,nga))  .or.
     >    (.not.rtdb_put(rtdb,'nwpw_APC:gamma',mt_dbl,nga,agamma))) then
         call errquit('nwpw_input: error writing to rtdb',8800,RTDB_ERR)
         end if
         end if
      end if

      goto 10

c
c socket ipi_client ip:port
c
 8900 if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,.true.)
         if (move) then
            ind = index(zone_name,' ') - 1
            value=rtdb_cput(rtdb,'nwpw:socket_type',1,
     >                      zone_name(1:ind))
            if (inp_a(ipname)) then
               ind = index(ipname,' ') - 1
               value = value.and.
     >                 rtdb_cput(rtdb,'nwpw:socket_ip',1,
     >                           ipname(1:ind))
            end if
         else
            value=rtdb_delete(rtdb,'nwpw:socket_type')
         end if
         if (.not.value)
     >      call errquit(
     >      'nwpw_input: error writing to rtdb', 8900, RTDB_ERR)
      end if
      goto 10

c
c
c  vfield file1.ascii file2.ascii ....
c  vfield off
c
 9000 ncut = 0
      do while (inp_a(zone_name)) 
         move = nwpw_parse_boolean(zone_name,.true.)
         if (move) then
            ind = index(zone_name,' ')
            if (ncut.eq.0) then
               vnames(1:ind) = zone_name(1:ind)
            else 
               vnames = vnames(1:ncut)//zone_name(1:ind)
            end if
            ncut = ncut + ind
            value=rtdb_cput(rtdb,'nwpw:vfield_filenames',
     >                      1,vnames(1:ncut))
         else
            value=rtdb_delete(rtdb,'nwpw:vfield_filenames')
         end if
         if (.not.value)
     >      call errquit(
     >      'nwpw_input: error writing to rtdb',9000,RTDB_ERR)
      end do
      goto 10

c
c  born [dielc] on/off [r1 r2 r3...] 
c
 9010 n1 = 0

      if (inp_f(Gc)) then
         if (.not.rtdb_put(rtdb,'nwpw:born_dielec',mt_dbl,1,Gc))
     >   call errquit('nwpw_input: error writing to rtdb',9010,RTDB_ERR)
      end if

      move = .true.
      if (inp_a(zone_name)) move = nwpw_parse_boolean(zone_name,.true.)

      if (.not.rtdb_put(rtdb,'nwpw:born',mt_log,1,move))
     > call errquit('nwpw_input: error writing to rtdb',9010,RTDB_ERR)

      !**** try reading in atom radii ****
      if (move) then
         if (.not.BA_push_get(mt_dbl,nw_max_atom,'indx_start1',
     >                        index_start1(2),index_start1(1)))
     >   call errquit('nwpw_input:;push index_start1',0,MA_ERR)

         do while (inp_f(Gc)) 
            dbl_mb(index_start1(1)+n1) = Gc/0.529177d0
            n1 = n1 + 1
         end do
         if (n1.gt.0) then
            if (.not.rtdb_put(rtdb,'nwpw:born_vradii',mt_dbl,
     >                        n1,dbl_mb(index_start1(1))))
     >      call errquit('nwpw_input: error writing to rtdb',
     >                   9010,RTDB_ERR)
         end if
         if (.not.BA_pop_stack(index_start1(2))) call errquit(
     >   'nwpw_input:failed deallocating index_start1',0,MA_ERR)
      end if
      goto 10

c
c symmetry  .true.||.false
c
 9020 if (inp_a(zone_name)) then
         move = nwpw_parse_boolean(zone_name,.true.)
      else
         move = .true.
      end if
      if (.not.rtdb_put(rtdb,'nwpw:rho_use_symmetry',mt_log,1,
     >                  move))
     >  call errquit('nwpw_input: error writing to rtdb',9020,RTDB_ERR)
      goto 10

c
c one-electron-guess it_in it_out it_ortho
c
 9030 if (.not.inp_i(nx)) nx = 50
      if (.not.inp_i(ny)) ny = 1
      if (.not.inp_i(nz)) nz = 1
      if (nz.ge.0) then
         if (.not.rtdb_put(rtdb,'nwpw:H1_it_in',mt_int,1,nx))
     >   call errquit('nwpw_input:error writing to rtdb',9030,RTDB_ERR)
         if (.not.rtdb_put(rtdb,'nwpw:H1_it_out',mt_int,1,ny))
     >   call errquit('nwpw_input:error writing to rtdb',9030,RTDB_ERR)
         if (.not.rtdb_put(rtdb,'nwpw:H1_it_ortho',mt_int,1,nz))
     >   call errquit('nwpw_input:error writing to rtdb',9030,RTDB_ERR)
      end if

      goto 10


************************************************************************
c
c END/Finished
c
 9999 continue
      return
      end

