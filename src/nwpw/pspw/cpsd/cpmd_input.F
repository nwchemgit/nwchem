c
c $Id$
c


      subroutine cpmd_input(rtdb)
      implicit none
      integer rtdb
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "errquit.fh"
c
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 61)
 
      character*32 dirs(num_dirs)
      character*255 test, id

      data dirs / 'cell_name:',
     >            'cell_name',
     >            'geometry_optimize',
     >            'formatted_filename:',
     >            'formatted_filename ',
     >            'input_wavefunction_filename:',
     >            'input_wavefunction_filename ',
     >            'output_wavefunction_filename:',
     >            'output_wavefunction_filename',
     >            'fake_mass:',
     >            'fake_mass',
     >            'time_step:',
     >            'time_step',
     >            'loop:',
     >            'loop',
     >            'scaling:',
     >            'scaling',
     >            'energy_cutoff:',
     >            'energy_cutoff',
     >            'wavefunction_cutoff:',
     >            'wavefunction_cutoff',
     >            'input_v_wavefunction_filename:',
     >            'input_v_wavefunction_filename',
     >            'output_v_wavefunction_filename:',
     >            'output_v_wavefunction_filename',
     >            'ewald_rcut:',
     >            'ewald_rcut',
     >            'ewald_ncut:',
     >            'ewald_ncut',
     >            'xyz_filename:',
     >            'xyz_filename',
     >            'exchange_correlation:',
     >            'exchange_correlation',
     >            'xc',
     >            'fractional_coordinates',
     >            'Nose-Hoover:',
     >            'Nose-Hoover',
     >            'emotion_filename:',
     >            'emotion_filename',
     >            'eigmotion_filename:',
     >            'eigmotion_filename',
     >            'hmotion_filename:',
     >            'hmotion_filename',
     >            'omotion_filename:',
     >            'omotion_filename',
     >            'ion_motion_filename:',
     >            'ion_motion_filename',
     >            'mulliken',
     >            'energy',
     >            'SA_decay',
     >            'Fei',
     >            'Fei_quench',
     >            'pressure',
     >            'cutoff',
     >            'initial_velocities',
     >            'reset_Temperature',
     >            'rotation',
     >            'translation',
     >            'dipole_motion',
     >            'temperature',
     >            'end'/

       character*50 cell_name
       character*50 input_wavefunction_filename
       character*50 output_wavefunction_filename
       character*50 input_v_wavefunction_filename
       character*50 output_v_wavefunction_filename
       character*50 xyz_filename
       character*50 emotion_filename
       character*50 eigmotion_filename
       character*50 hmotion_filename
       character*50 omotion_filename
       character*50 ion_motion_filename
       character*50 exchange_correlation,zone_name
       logical geometry_optimize,frac_coord,mulliken,mulliken_kawai
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut,mchain,nchain
       double precision scaling(2),ecut,wcut
       logical          nose,nosers,move
       double precision Pe,Te,fe
       double precision Pr,Tr
       integer index_start1(2),n1,j,jstart,jlast,jstride

       logical  nwpw_parse_boolean
       external nwpw_parse_boolean
        

*     ***** initializations ****
c      call cpmd_input_default(rtdb)
      npsp = 0


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cpmd_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cpmd_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cpmd_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,100,  200, 300,300, 400,400, 500,500, 600,600, 
     >       700,700,  800,800, 900,900, 1000,1000, 1100,1100,
     >       1200,1200, 1300,1300, 1400,1400, 1500,1500, 1600,1600,
     >       1700,1700,1700, 1800, 1900,1900, 2000,2000, 2100,2100,
     >       2200,2200, 2300,2300, 2400,2400,
     >       2500,2600,2700,2800,2900,3000,1150,3100,3200,
     >       3300,3400,3500,3600,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind,
     &       INPUT_ERR)

c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cpmd_input: failed to read cell_name', 0, INPUT_ERR)
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:cell_name',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing cell_name', 100, RTDB_ERR)
      goto 10
c
c geometry_optimize
c
 200  geometry_optimize = .true.
      value = rtdb_put(rtdb,'cpmd:geometry_optimize',mt_log,1,
     >                       geometry_optimize)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 200, RTDB_ERR)
      goto 10
c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cpmd_input: failed to read psp_filename', 0, INPUT_ERR)
      npsp = npsp + 1
      id = 'cpmd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cpmd_input: rtdb_cput failed', 0, RTDB_ERR)
      value = rtdb_put(rtdb,'cpmd:npsp',     mt_int,1,npsp)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 300, RTDB_ERR)
      goto 10
c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read input_wavefunctions', 0,
     &       INPUT_ERR)
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 400, RTDB_ERR)
      goto 10
c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read output_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 500, RTDB_ERR)
      goto 10
c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'cpmd_input: failed to read fake_mass', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 600, RTDB_ERR)
      goto 10
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'cpmd_input: failed to read time_step', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:time_step',mt_dbl,1,time_step)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 700, RTDB_ERR)
      goto 10
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'cpmd_input: failed to read loop', 0, INPUT_ERR)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'cpmd_input: failed to read loop', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:loop',mt_int,2,loop)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 800, RTDB_ERR)
      goto 10
c
c  scaling
c 
 900  if (.not. inp_f(scaling(1)))
     >     call errquit(
     >          'cpmd_input: failed to read scaling', 0, INPUT_ERR)
      if (.not. inp_f(scaling(2)))
     >     call errquit(
     >          'cpmd_input: failed to read scaling', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:scaling',mt_dbl,2,scaling)

      !*** try reading in atom indexes ****
      if (.not.BA_push_get(mt_int,nw_max_atom,'indx_start1',
     >                    index_start1(2),index_start1(1)))
     >    call errquit(
     > 'cpmd_input:failed allocating index_start1',0,MA_ERR)
      n1 = 0
      do while (inp_irange(jstart,jlast,jstride))
         do j=jstart,jlast,jstride               
          int_mb(index_start1(1)+n1) = j
          n1 = n1+1
         end do
      end do
      value = value.and.rtdb_put(rtdb,'nwpw:scaling_natms',
     >                              mt_int,1,n1)
      if (n1.gt.0) then
         value = value.and.rtdb_put(rtdb,'nwpw:scaling_atoms',
     >                              mt_int,n1,int_mb(index_start1(1)))
      end if

      if (.not.BA_pop_stack(index_start1(2))) call errquit(
     > 'cpmd_input:failed deallocating index_start1',0,MA_ERR)

      if (.not.value)
     > call errquit('cpmd_input: writing ', 900, RTDB_ERR)
      goto 10
c
c  energy_cutoff
c 
 1000 if (.not. inp_f(ecut))
     >     call errquit(
     >          'cpmd_input: failed to read ecut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1000, RTDB_ERR)
      goto 10
c
c  wavefunction_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cpmd_input: failed to read wcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1100, RTDB_ERR)
      goto 10

c
c  cutoff
c
 1150 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cpmd_input: failed to read wcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:wcut',mt_dbl,1,wcut)
      value = rtdb_put(rtdb,'cpmd:ecut',mt_dbl,1,2.0d0*wcut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1150, RTDB_ERR)
      goto 10

c
c  input_v_wavefunction_filename
c 
 1200 if (.not. inp_a(input_v_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read input_v_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(input_v_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:input_v_wavefunction_filename',
     >                    1,input_v_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1200, RTDB_ERR)

      goto 10
c
c  output_v_wavefunction_filename
c 
 1300 if (.not. inp_a(output_v_wavefunction_filename))
     >     call errquit(
     >          'cpmd_input: failed to read output_v_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(output_v_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:output_v_wavefunction_filename',
     >                    1,output_v_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1300, RTDB_ERR)
      goto 10
c
c  ewald_rcut
c 
 1400 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cpmd_input: failed to read rcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1400, RTDB_ERR)
      goto 10
c
c  ewald_ncut
c 
 1500 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cpmd_input: failed to read ncut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpmd:ncut',mt_int,1,ncut)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1500, RTDB_ERR)
      goto 10
c
c  xyz_filename
c 
 1600 if (.not. inp_a(xyz_filename))
     >     call errquit(
     >          'cpmd_input: failed to read xyz_filename', 0, INPUT_ERR)
      ind = index(xyz_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:xyz_filename',
     >                    1,xyz_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1600, RTDB_ERR)
      goto 10
c
c  exchange_correlation
c 
 1700  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cpmd_input: failed to read exchange_correlation', 0,
     &       INPUT_ERR)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1700, RTDB_ERR)
      goto 10
c
c fractional_coordinates
c
 1800 frac_coord = .true.
      value = rtdb_put(rtdb,'cpmd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1800, RTDB_ERR)
      goto 10
c
c Nose-Hoover
c
 1900 nose   = .true.
      nosers = .true.
      Pe = 1200.00d0
      Pr = 1200.00d0
      Te = 298.15d0
      Tr = 298.15d0
      mchain = 1
      nchain = 1
      if (.not. inp_f(fe)) goto 1901   
      Pe = fe
      if (.not. inp_f(fe)) goto 1901   
      Te = fe
      if (.not. inp_f(fe)) goto 1901   
      Pr = fe
      if (.not. inp_f(fe)) goto 1901   
      Tr = fe
      if (.not.inp_a(zone_name)) goto 1901
      if (inp_compare(.false.,zone_name,'restart'))    nosers=.true.
      if (inp_compare(.false.,zone_name,'start'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'clear'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'initialize')) nosers=.false.
      if (.not.inp_i(mchain)) goto 1901
      if (.not.inp_i(nchain)) goto 1901
 1901 value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:nose_restart',mt_log,1,nosers)
      value = value.and.rtdb_put(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value.and.rtdb_put(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value.and.rtdb_put(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value.and.rtdb_put(rtdb,'cpmd:Tr',mt_dbl,1,Tr)
      value = value.and.rtdb_put(rtdb,'cpmd:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'cpmd:Nchain',mt_int,1,nchain)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 1900, RTDB_ERR)
      goto 10
c
c  emotion_filename
c 
 2000 if (.not. inp_a(emotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read emotion_filename', 0,
     &       INPUT_ERR)
      ind = index(emotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:emotion_filename',
     >                    1,emotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2000, RTDB_ERR)
      goto 10
c
c  eigmotion_filename
c 
 2100 if (.not. inp_a(eigmotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read eigmotion_filename', 0,
     &       INPUT_ERR)
      ind = index(eigmotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:eigmotion_filename',
     >                    1,eigmotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2100, RTDB_ERR)
      goto 10
c
c  hmotion_filename
c 
 2200 if (.not. inp_a(hmotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read hmotion_filename', 0,
     &       INPUT_ERR)
      ind = index(hmotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:hmotion_filename',
     >                    1,hmotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2200, RTDB_ERR)
      goto 10
c
c  omotion_filename
c 
 2300 if (.not. inp_a(omotion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read omotion_filename', 0,
     &       INPUT_ERR)
      ind = index(omotion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:omotion_filename',
     >                    1,omotion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2300, RTDB_ERR)
      goto 10
c
c  ion_motion_filename
c 
 2400 if (.not. inp_a(ion_motion_filename))
     >     call errquit(
     >          'cpmd_input: failed to read ion_motion_filename', 0,
     &       INPUT_ERR)
      ind = index(ion_motion_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpmd:ion_motion_filename',
     >                    1,ion_motion_filename(1:ind))
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2400, RTDB_ERR)
      goto 10
c
c Mulliken
c
 2500 if (.not.inp_a(zone_name)) zone_name = 'lcao'
      mulliken       = .true.
      mulliken_kawai = .false.
      ind = index(zone_name,' ') - 1
      if (inp_compare(.false.,zone_name,'kawai')) mulliken_kawai=.true.
      value = rtdb_put(rtdb,'cpmd:mulliken',mt_log,1,
     >                       mulliken)
     >   .and.rtdb_put(rtdb,'nwpw:mulliken_kawai',mt_log,1,
     >                       mulliken_kawai)
      if (.not.value)
     >  call errquit(
     >       'cpmd_input: error writing to rtdb', 2500, RTDB_ERR)
      goto 10

c
c energy
c
 2600 nose = .false.
      value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 2600, RTDB_ERR)
      goto 10
c
c SA_decay
c
 2700 if (inp_f(scaling(1))) then
         scaling(2) = scaling(1)
         if (.not.inp_f(scaling(2))) scaling(2) = scaling(1)
         value = rtdb_put(rtdb,'cpmd:sa_decay',mt_dbl,2,scaling)
      else
         value = rtdb_delete(rtdb,'cpmd:sa_decay')
      end if
      if (.not.value)
     > call errquit('cpmd_input: writing sa_decay', 2700, RTDB_ERR)
      goto 10
c
c Fei
c
 2800 if (inp_a(ion_motion_filename)) then
         ind = index(ion_motion_filename,' ') - 1
         value = rtdb_cput(rtdb,'cpmd:fei_filename',
     >                    1,ion_motion_filename(1:ind))
         value = rtdb_cput(rtdb,'nwpw:fei_filename',
     >                    1,ion_motion_filename(1:ind))
          if (.not.value)
     >    call errquit('cpmd_input: writing ', 2800, RTDB_ERR)
      end if
      value = .true.
      if (.not.rtdb_put(rtdb,'cpmd:fei',mt_log,1,value))
     > call errquit('cpmd_input: writing fei ', 2800, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'nwpw:fei',mt_log,1,value))
     > call errquit('cpmd_input: writing fei ', 2800, RTDB_ERR)
      goto 10
c
c Fei Quench
c
 2900 value = .true.
      if (.not.rtdb_put(rtdb,'cpmd:fei_quench',mt_log,1,value))
     > call errquit('cpmd_input: writing fei quench ', 2900, RTDB_ERR)
      goto 10

c
c pressure
c
 3000 if (.not.inp_a(zone_name)) zone_name = 'on'
      value = .true.
      value = nwpw_parse_boolean(zone_name,value)
c      if (inp_compare(.false.,zone_name,'on'))      value=.true.
c      if (inp_compare(.false.,zone_name,'.true.'))  value=.true.
c      if (inp_compare(.false.,zone_name,'off'))     value=.false.
c      if (inp_compare(.false.,zone_name,'.false.')) value=.false.
      if (.not.rtdb_put(rtdb,'cpmd:pressure',mt_log,1,value))
     > call errquit('cpmd_input: writing pressure ', 3000, RTDB_ERR)
      goto 10

c
c initial_velocities temperature seed
c  rtdb = cpmd:init_velocities_temperature
c
 3100 if (inp_f(scaling(1))) then            !*** in units of K ***
         if (.not.inp_i(n1)) n1 = 494

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities_temperature',
     >                  mt_dbl,1,scaling))
     >   call errquit('cpmd_input: initial temperature',3100,RTDB_ERR)

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities_seed',
     >                  mt_int,1,n1))
     >   call errquit('cpmd_input: initial temperature',3100,RTDB_ERR)

         if (.not.rtdb_put(rtdb,'nwpw:init_velocities',
     >                  mt_log,1,.true.))
     >   call errquit('cpmd_input: initial velocities',3100,RTDB_ERR)

      else
         value = rtdb_delete(rtdb,'nwpw:init_velocities')
      end if
      goto 10

 3200 value = rtdb_delete(rtdb,"nwpw:ke_total")
      value = rtdb_delete(rtdb,"nwpw:kg_total")
      value = rtdb_delete(rtdb,"nwpw:ke_count")
      goto 10


c
c  rotation
c
 3300 if (.not.inp_a(zone_name)) goto 3301
      move = nwpw_parse_boolean(zone_name,.false.)
      if (.not.rtdb_put(rtdb,'nwpw:rotation',mt_log,1,move))
     >  call errquit('cpmd_input: error writing to rtdb',3300,RTDB_ERR)

 3301 continue
      goto 10

c
c  translation
c
 3400 if (.not.inp_a(zone_name)) goto 3401
      move = nwpw_parse_boolean(zone_name,.true.)
      if (.not.rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,move))
     >  call errquit('cpmd_input: error writing to rtdb',3400,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'band:allow_translation',mt_log,1,move))
     >  call errquit('cpmd_input: error writing to rtdb',3400,RTDB_ERR)

 3401 continue
      goto 10

c
c dipole_motion
c
 3500 if (inp_a(ion_motion_filename)) then
         ind = index(ion_motion_filename,' ') - 1
         if (.not.rtdb_cput(rtdb,'nwpw:dipole_motion_filename',
     >                    1,ion_motion_filename(1:ind)))
     >   call errquit('nwpw_input: writing filename',3500,RTDB_ERR)
      end if
      if (.not.rtdb_put(rtdb,'nwpw:dipole_motion',mt_log,1,value))
     > call errquit('cpmd_input: writing dipole_motion',3500,RTDB_ERR)
      goto 10


c
c temperature Tion Pion Telc Pion Pelc start/restart nchain mchain
c - a reordering of the parameters for nose-hoover
c
 3600 nose   = .true.
      nosers = .true.
      Tr =  298.15d0
      Pr = 1200.00d0
      Te =  298.15d0
      Pe = 1200.00d0
      mchain = 1
      nchain = 1
      if (.not. inp_f(fe)) goto 3601
      Tr = fe
      Te = fe
      if (.not. inp_f(fe)) goto 3601
      Pr = fe
      Pe = fe
      if (.not. inp_f(fe)) goto 3601
      Te = fe
      if (.not. inp_f(fe)) goto 3601
      Pe = fe
 3601 if (.not.inp_a(zone_name)) goto 3602
      if (inp_compare(.false.,zone_name,'restart'))    nosers=.true.
      if (inp_compare(.false.,zone_name,'start'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'clear'))      nosers=.false.
      if (inp_compare(.false.,zone_name,'initialize')) nosers=.false.
 3602 if (.not.inp_i(nchain)) goto 3603
      mchain = nchain
      if (.not.inp_i(mchain)) goto 3603

 3603 value = rtdb_put(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value.and.
     >        rtdb_put(rtdb,'cpmd:nose_restart',mt_log,1,nosers)
      value = value.and.rtdb_put(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value.and.rtdb_put(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value.and.rtdb_put(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value.and.rtdb_put(rtdb,'cpmd:Tr',mt_dbl,1,Tr)
      value = value.and.rtdb_put(rtdb,'cpmd:Mchain',mt_int,1,mchain)
      value = value.and.rtdb_put(rtdb,'cpmd:Nchain',mt_int,1,nchain)
      if (.not.value)
     > call errquit('cpmd_input: writing ', 3600, RTDB_ERR)
      goto 10



*     ***** add wavefunction to rtdb ****
 9999 continue
      return
      end

      subroutine cpmd_input_default(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"

*     ***** local variables *****
      logical value
      integer ind               

      character*50 cell_name
      character*50 input_wavefunction_filename
      character*50 output_wavefunction_filename
      character*50 input_v_wavefunction_filename
      character*50 output_v_wavefunction_filename
      character*50 xyz_filename
      character*50 emotion_filename
      character*50 eigmotion_filename
      character*50 hmotion_filename
      character*50 omotion_filename
      character*50 ion_motion_filename
      character*50 exchange_correlation
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

      !exchange_correlation = 'vosko'
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

      call util_file_prefix('xyz',xyz_filename)
      call util_file_prefix('emotion',emotion_filename)
      call util_file_prefix('eigmotion',eigmotion_filename)
      call util_file_prefix('hmotion',hmotion_filename)
      call util_file_prefix('omotion',omotion_filename)
      call util_file_prefix('ion_motion',ion_motion_filename)
c      xyz_filename = 'XYZ'
c      emotion_filename    = 'EMOTION'
c      eigmotion_filename  = 'EIGMOTION'
c      hmotion_filename    = 'HMOTION'
c      omotion_filename    = 'OMOTION'
c      ion_motion_filename = 'MOTION'

      nose = .false.
      Pe = 1200.00d0
      Pr = 1200.00d0
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

c      ind = index(exchange_correlation,' ') - 1
c      value = value.and.
c     >        rtdb_cput(rtdb,'cpmd:exchange_correlation',
c     >                    1,exchange_correlation(1:ind))

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
     >       'cpmd_default: error writing to rtdb', 0, RTDB_ERR)

      return
      end

