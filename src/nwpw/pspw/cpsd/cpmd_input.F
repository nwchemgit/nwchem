c
c $Id: cpmd_input.F,v 1.4 2001-11-30 17:30:29 bylaska Exp $
c


      subroutine cpmd_input(rtdb)
      implicit none
      integer rtdb
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 27)
 
      character*32 dirs(num_dirs)
      character*255 test, id

      data dirs / 'cell_name:',
     >            'geometry_optimize',
     >            'formatted_filename:',
     >            'input_wavefunction_filename:',
     >            'output_wavefunction_filename:',
     >            'fake_mass:',
     >            'time_step:',
     >            'loop:',
     >            'scaling:',
     >            'energy_cutoff:',
     >            'wavefunction_cutoff:',
     >            'input_v_wavefunction_filename:',
     >            'output_v_wavefunction_filename:',
     >            'ewald_rcut:',
     >            'ewald_ncut:',
     >            'xyz_filename:',
     >            'exchange_correlation:',
     >            'fractional_coordinates',
     >            'Nose-Hoover:',
     >            'emotion_filename:',
     >            'eigmotion_filename:',
     >            'hmotion_filename:',
     >            'omotion_filename:',
     >            'ion_motion_filename:',
     >            'mulliken',
     >            'energy',
     >            'end'/

       character*30 cell_name
       character*30 input_wavefunction_filename
       character*30 output_wavefunction_filename
       character*30 input_v_wavefunction_filename
       character*30 output_v_wavefunction_filename
       character*30 xyz_filename
       character*30 emotion_filename
       character*30 eigmotion_filename
       character*30 hmotion_filename
       character*30 omotion_filename
       character*30 ion_motion_filename
       character*30 exchange_correlation
       logical geometry_optimize,frac_coord,mulliken
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut
       double precision scaling(2),ecut,wcut
       logical          nose
       double precision Pe,Te,fe
       double precision Pr,Tr
        

*     ***** initializations ****
      call cpmd_input_default(rtdb)
      npsp = 0


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cpmd_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cpmd_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cpmd_input: unknown directive', 0)


      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900,
     >       1000,1100,1200,1300,1400,1500,1600,1700,1800,
     >       1900,2000,2100,2200,2300,2400,2500,2600,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind)

c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cpmd_input: failed to read cell_name', 0)
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:cell_name',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing cell_name', 100)
      goto 10
c
c geometry_optimize
c
 200  geometry_optimize = .true.
      value = rtdb_put(rtdb,'cpmd:geometry_optimize',mt_log,1,
     >                       geometry_optimize)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 200)
      goto 10
c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cpmd_input: failed to read psp_filename', 0)
      npsp = npsp + 1
      id = 'cpmd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cpmd_input: rtdb_cput failed', 0)
      value = rtdb_put(rtdb,'cpmd:npsp',     mt_int,1,npsp)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 300)
      goto 10
c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read input_wavefunctions', 0)
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 400)
      goto 10
c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read output_wavefunction', 0)
      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 500)
      goto 10
c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'cpmd_input: failed to read fake_mass', 0)
      value = rtdb_put(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 600)
      goto 10
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'cpmd_input: failed to read time_step', 0)
      value = rtdb_put(rtdb,'cpmd:time_step',mt_dbl,1,time_step)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 700)
      goto 10
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'cpmd_input: failed to read loop', 0)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'cpmd_input: failed to read loop', 0)
      value = rtdb_put(rtdb,'cpmd:loop',mt_int,2,loop)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 800)
      goto 10
c
c  scaling
c 
 900  if (.not. inp_f(scaling(1)))
     >     call errquit(
     >          'cpmd_input: failed to read scaling', 0)
      if (.not. inp_f(scaling(2)))
     >     call errquit(
     >          'cpmd_input: failed to read scaling', 0)
      value = rtdb_put(rtdb,'cpmd:scaling',mt_dbl,2,scaling)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 900)
      goto 10
c
c  energy_cutoff
c 
 1000 if (.not. inp_f(ecut))
     >     call errquit(
     >          'cpmd_input: failed to read ecut', 0)
      value = rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1000)
      goto 10
c
c  wavefunction_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cpmd_input: failed to read wcut', 0)
      value = rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1100)
      goto 10
c
c  input_v_wavefunction_filename
c 
 1200 if (.not. inp_a(input_v_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read input_v_wavefunction', 0)
      ind = index(input_v_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:input_v_wavefunction_filename',
     >                    1,input_v_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1200)

      goto 10
c
c  output_v_wavefunction_filename
c 
 1300 if (.not. inp_a(output_v_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read output_v_wavefunction', 0)
      ind = index(output_v_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:output_v_wavefunction_filename',
     >                    1,output_v_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1300)
      goto 10
c
c  ewald_rcut
c 
 1400 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cpmd_input: failed to read rcut', 0)
      value = rtdb_put(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1400)
      goto 10
c
c  ewald_ncut
c 
 1500 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cpmd_input: failed to read ncut', 0)
      value = rtdb_put(rtdb,'cpmd:ncut',mt_int,1,ncut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1500)
      goto 10
c
c  xyz_filename
c 
 1600 if (.not. inp_a(xyz_filename))
     >     call errquit(
     >          'cpmd_input: failed to read xyz_filename', 0)
      ind = index(xyz_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:xyz_filename',
     >                    1,xyz_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1600)
      goto 10
c
c  exchange_correlation
c 
 1700  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cpmd_input: failed to read exchange_correlation', 0)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1700)
      goto 10
c
c fractional_coordinates
c
 1800 frac_coord = .true.
      value = rtdb_put(rtdb,'cpmd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1800)
      goto 10
c
c Nose-Hoover
c
 1900 nose = .true.
      if (.not. inp_f(fe)) goto 1901   
      Pe = fe
      if (.not. inp_f(fe)) goto 1901   
      Te = fe
      if (.not. inp_f(fe)) goto 1901   
      Pr = fe
      if (.not. inp_f(fe)) goto 1901   
      Tr = fe
 1901 value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value.and.rtdb_put(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value.and.rtdb_put(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value.and.rtdb_put(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value.and.rtdb_put(rtdb,'cpmd:Tr',mt_dbl,1,Tr)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1900)
      goto 10
c
c  emotion_filename
c 
 2000 if (.not. inp_a(emotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read emotion_filename', 0)
      ind = index(emotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:emotion_filename',
     >                    1,emotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2000)
      goto 10
c
c  eigmotion_filename
c 
 2100 if (.not. inp_a(eigmotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read eigmotion_filename', 0)
      ind = index(eigmotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:eigmotion_filename',
     >                    1,eigmotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2100)
      goto 10
c
c  hmotion_filename
c 
 2200 if (.not. inp_a(hmotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read hmotion_filename', 0)
      ind = index(hmotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:hmotion_filename',
     >                    1,hmotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2200)
      goto 10
c
c  omotion_filename
c 
 2300 if (.not. inp_a(omotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read omotion_filename', 0)
      ind = index(omotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:omotion_filename',
     >                    1,omotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2300)
      goto 10
c
c  ion_motion_filename
c 
 2400 if (.not. inp_a(ion_motion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read ion_motion_filename', 0)
      ind = index(ion_motion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:ion_motion_filename',
     >                    1,ion_motion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2400)
      goto 10
c
c Mulliken
c
 2500 mulliken = .true.
      value = rtdb_put(rtdb,'cpmd:mulliken',mt_log,1,
     >                       mulliken)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2500)
      goto 10
c
c energy
c
 2600 nose = .false.
      value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2600)
      goto 10

*     ***** add wavefunction to rtdb ****
 9999 continue
      return
      end

      subroutine cpmd_input_default(rtdb)
      implicit none
      integer rtdb
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

*     ***** local variables *****
      logical value
      integer ind               

      character*30 cell_name
      character*30 input_wavefunction_filename
      character*30 output_wavefunction_filename
      character*30 input_v_wavefunction_filename
      character*30 output_v_wavefunction_filename
      character*30 xyz_filename
      character*30 emotion_filename
      character*30 eigmotion_filename
      character*30 hmotion_filename
      character*30 omotion_filename
      character*30 ion_motion_filename
      character*30 exchange_correlation
      logical geometry_optimize,frac_coord,mulliken
      double precision fake_mass,time_step,rcut
      integer loop(2),npsp,ncut
      double precision scaling(2),ecut,wcut
      logical          nose
      double precision Pe,Te
      double precision Pr,Tr
        

*     **** don't set defaults if they already exist ****
      value = rtdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut)
      if (value) return


*     ***** initializations ****
      cell_name             = 'cell_default'

      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)
c     input_wavefunction_filename  = ' '
c     output_wavefunction_filename = ' '

      call util_file_prefix('vmovecs',input_v_wavefunction_filename)
      call util_file_prefix('vmovecs',output_v_wavefunction_filename)
c     input_v_wavefunction_filename = ' '
c     output_v_wavefunction_filename = ' '

      exchange_correlation = 'vosko'
      geometry_optimize            = .true.
      frac_coord                   = .false.
      mulliken                     = .false.
      fake_mass = 1000.0d0
      time_step = 5.0d0
      loop(1) = 10
      loop(2) = 1
      scaling(1) = 1.0d0
      scaling(2) = 1.0d0
      ecut=9000.0d0
      wcut=9000.0d0
      rcut=0.0d0
      ncut=0
      npsp = 0
      xyz_filename = 'XYZ'
      emotion_filename    = 'EMOTION'
      eigmotion_filename  = 'EIGMOTION'
      hmotion_filename    = 'HMOTION'
      omotion_filename    = 'OMOTION'
      ion_motion_filename = 'MOTION'

      nose = .false.
      Pe = 100.00d0
      Pr = 100.00d0
      Te = 298.15d0
      Tr = 298.15d0

      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:cell_name',1,cell_name(1:ind))

      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(output_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))

      ind = index(input_v_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:input_v_wavefunction_filename',
     >                    1,input_v_wavefunction_filename(1:ind))

      ind = index(output_v_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:output_v_wavefunction_filename',
     >                    1,output_v_wavefunction_filename(1:ind))

      ind = index(xyz_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:xyz_filename',
     >                    1,xyz_filename(1:ind))

      ind = index(exchange_correlation,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      ind = index(emotion_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:emotion_filename',
     >                    1,emotion_filename(1:ind))
      ind = index(eigmotion_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:eigmotion_filename',
     >                    1,eigmotion_filename(1:ind))
      ind = index(hmotion_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:hmotion_filename',
     >                    1,hmotion_filename(1:ind))
      ind = index(omotion_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:omotion_filename',
     >                    1,omotion_filename(1:ind))
      ind = index(ion_motion_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cpmd:ion_motion_filename',
     >                    1,ion_motion_filename(1:ind))

      value = value.and.
     >        rtdb_put(rtdb,'cpmd:geometry_optimize',mt_log,1,
     >                       geometry_optimize)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:mulliken',mt_log,1,
     >                       mulliken)

      value = value.and.
     >        rtdb_put(rtdb,'cpmd:npsp',     mt_int,1,npsp)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:scaling',mt_dbl,2,scaling)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:ncut',mt_int,1,ncut)

      value = value.and.
     >        rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:Tr',mt_dbl,1,Tr)

*     ***** Error writing to RTDB *****
      if (.not.value)
     >  call errquit(
     >       'cpmd_default: error writing to rtdb', 0)

      return
      end

