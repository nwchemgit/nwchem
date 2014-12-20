c
c $Id$
c
      subroutine cpsd_input(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 34)
 
      character*30 dirs(num_dirs)
      character*255 test, id

      data dirs / 'cell_name:',
     >            'cell_name',
     >            'geometry_optimize',
     >            'formatted_filename:',
     >            'formatted_filename',
     >            'input_wavefunction_filename:',
     >            'input_wavefunction_filename',
     >            'output_wavefunction_filename:',
     >            'output_wavefunction_filename',
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
     >            'fractional_coordinates',
     >            'mulliken',
     >            'mult',
     >            'dft',
     >            'odft',
     >            'end'/

       character*50 cell_name
       character*50 input_wavefunction_filename
       character*50 output_wavefunction_filename
       character*50 exchange_correlation
       logical geometry_optimize,frac_coord,mulliken
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut
       double precision tolerances(3),ecut,wcut,fe
       integer ispin,multiplicity

*     **** initializations ****
c      call cpsd_input_default(rtdb)
      npsp = 0



 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cpsd_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cpsd_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cpsd_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,100, 200, 300,300, 400,400, 500,500, 600,600, 
     >       700,700, 800,800, 900,900, 1000,1000, 1100,1100,
     >       1200,1200, 1300,1300, 1400,1400,1400,
     >       1500,1600,1700,1800,1900,
     >      9999) ind
      call errquit(
     >      'cpsd_input: unimplemented directive', ind, INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cpsd_input: failed to read cell_name', 0, INPUT_ERR)
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:cell_name',1,cell_name(1:ind))
      if (.not.value)
     > call errquit('cpsd_input: writing cell_name', 0, RTDB_ERR)
      goto 10


c
c geometry_optimize
c
 200  geometry_optimize = .true.
      value = rtdb_put(rtdb,'cpsd:geometry_optimize',mt_log,1,
     >                       geometry_optimize)
      if (.not.value)
     > call errquit('cpsd_input: writing geometry_optimize', 0,
     &       RTDB_ERR)
      goto 10


c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cpsd_input: failed to read psp_filename', 0, INPUT_ERR)
      npsp = npsp + 1
      id = 'cpsd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cpsd_input: rtdb_cput failed', 0, RTDB_ERR)
      value = rtdb_put(rtdb,'cpsd:npsp',     mt_int,1,npsp)
      if (.not.value)
     > call errquit('cpsd_input: writing npsp', 0, RTDB_ERR)
      goto 10


c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cpsd_input: failed to read input_wavefunction', 0,
     &       INPUT_ERR)

      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpsd_input: writing input vectors', 0, RTDB_ERR)
      goto 10

c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cpsd_input: failed to read output_wavefunction', 0,
     &       INPUT_ERR)

      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value)
     > call errquit('cpsd_input: writing output vectors', 0, RTDB_ERR)
      goto 10


c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'wvfnc_init_input: failed to read nup', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value)
     > call errquit('cpsd_input: writing fake_mass', 0, RTDB_ERR)
      goto 10
 
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:time_step',mt_dbl,1,time_step)
      if (.not.value)
     > call errquit('cpsd_input: writing time_step', 0, RTDB_ERR)
      goto 10
 
 
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0, INPUT_ERR)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:loop',mt_int,2,loop)
      if (.not.value)
     > call errquit('cpsd_input: writing loop', 0, RTDB_ERR)
      goto 10
 
 
c
c  tolerances
c 
 900  tolerances(1) = 1.0e-9
      tolerances(2) = 1.0e-9
      tolerances(3) = 1.0e-4
      if (inp_f(fe)) then
         tolerances(1) = fe
         if (inp_f(fe)) then
            tolerances(2) = fe
            if (inp_f(fe)) tolerances(3) = fe
         end if
      end if
 901  value = rtdb_put(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances)
      if (.not.value)
     > call errquit('cpsd_input: writing tolerances', 0, RTDB_ERR)
      goto 10

c
c  energy_cutoff
c 
 1000 if (.not. inp_f(ecut))
     >     call errquit(
     >          'cpsd_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
      if (.not.value)
     > call errquit('cpsd_input: writing ecut', 0, RTDB_ERR)
      goto 10

c
c  wavefunction_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cpsd_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
      if (.not.value)
     > call errquit('cpsd_input: writing wcut', 0, RTDB_ERR)
      goto 10

c
c  ewald_rcut
c 
 1200 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cpsd_input: failed to read rcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:rcut',mt_dbl,1,rcut)
      if (.not.value)
     > call errquit('cpsd_input: writing rcut', 0, RTDB_ERR)
      goto 10
c
c  ewald_ncut
c 
 1300 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cpsd_input: failed to read ncut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cpsd:ncut',mt_int,1,ncut)
      if (.not.value)
     > call errquit('cpsd_input: writing ncut', 0, RTDB_ERR)
      goto 10
c
c  exchange_correlation
c 
 1400  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cpsd_input: failed to read exchange_correlation', 0,
     &       INPUT_ERR)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      if (.not.value)
     > call errquit('cpsd_input: writing xc', 0, RTDB_ERR)
      goto 10

c
c fractional coordinates
c
 1500 frac_coord = .true.
      value = rtdb_put(rtdb,'cpsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      if (.not.value)
     > call errquit('cpsd_input: writing frac', 0, RTDB_ERR)
      goto 10


c
c Mulliken
c
 1600 mulliken = .true.
      value = rtdb_put(rtdb,'cpsd:mulliken',mt_log,1,
     >                       mulliken)
      if (.not.value)
     > call errquit('cpsd_input: writing mulliken', 0, RTDB_ERR)
      goto 10
c
c  multiplicity
c
 1700 if (.not. inp_i(multiplicity))
     >     call errquit(
     >          'cgsd_input: failed to read mult', 0, INPUT_ERR)
      if (multiplicity.lt.1) then
         multiplicity = 1
         ispin        = 1
      else if (multiplicity.gt.1) then
         ispin        = 2
      end if
      value = rtdb_put(rtdb,'cpsd:mult',mt_int,1,multiplicity)
      if (.not.value)
     > call errquit('cpsd_input: writing mult', 0, RTDB_ERR)
      goto 10
c
c  dft
c
 1800 ispin = 1
      value = rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
      if (.not.value)
     > call errquit('cpsd_input: writing ispin', 0, RTDB_ERR)
      goto 10
c
c  odft
c
 1900 ispin = 2
      value = rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
      if (.not.value)
     > call errquit('cpsd_input: writing ispin', 0, RTDB_ERR)
      goto 10



*     ***** add wavefunction to rtdb ****
 9999 continue



      return
      end


c      subroutine cpsd_input_default(rtdb)
c      implicit none
c#include "errquit.fh"
c      integer rtdb
c
c#include "inp.fh"
c#include "bafdecls.fh"
c#include "rtdb.fh"
cc
c
c*     **** local variables ****
c      logical value
c      integer ind             
c
c      character*50 cell_name
c      character*50 input_wavefunction_filename
c      character*50 output_wavefunction_filename
c      character*50 exchange_correlation
c      logical geometry_optimize,frac_coord,mulliken
c      double precision fake_mass,time_step,rcut
c      integer loop(2),npsp,ncut
c      double precision tolerances(3),ecut,wcut
c      integer ispin,multiplicity
c
c
c*     **** don't set defaults if they already exist ****
c      value = rtdb_get(rtdb,'cpsd:ncut',mt_int,1,ncut)
c      if (value) return
c
c*     ***** initializations ****
c      cell_name             = 'cell_default'
c
c      call util_file_prefix('movecs',input_wavefunction_filename)
c      call util_file_prefix('movecs',output_wavefunction_filename)
cc     input_wavefunction_filename  = ' '
cc     output_wavefunction_filename = ' '
c
c      exchange_correlation         = 'vosko'
c      geometry_optimize            = .false.
c      frac_coord                   = .false.
c      mulliken                     = .false.
c      fake_mass = 400000.0d0
c      time_step = 5.8d0
c      loop(1) = 10
c      loop(2) = 1
c      tolerances(1) = 1.0d-9
c      tolerances(2) = 1.0d-9
c      tolerances(3) = 1.0d-4
c      ecut=9000.0d0
c      wcut=9000.0d0
c      rcut = 0.0d0
c      ncut = 0
c      npsp = 0
c      ispin        = 1
c      multiplicity = 1
c
c
c      ind = index(cell_name,' ') - 1
c      value = rtdb_cput(rtdb,'cpsd:cell_name',1,cell_name(1:ind))
c
c      ind = index(input_wavefunction_filename,' ') - 1
c      value = value.and.
c     >        rtdb_cput(rtdb,'cpsd:input_wavefunction_filename',
c     >                    1,input_wavefunction_filename(1:ind))
c
c      ind = index(output_wavefunction_filename,' ') - 1
c      value = value.and.
c     >        rtdb_cput(rtdb,'cpsd:output_wavefunction_filename',
c     >                    1,output_wavefunction_filename(1:ind))
c
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:geometry_optimize',mt_log,1,
c     >                       geometry_optimize)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:fractional_coordinates',mt_log,1,
c     >                       frac_coord)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:mulliken',mt_log,1,
c     >                       mulliken)
c
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:npsp',     mt_int,1,npsp)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:time_step',mt_dbl,1,time_step)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:loop',mt_int,2,loop)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:rcut',mt_dbl,1,rcut)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:ncut',mt_int,1,ncut)
cc      ind = index(exchange_correlation,' ') - 1
cc      value = value.and.
cc     >        rtdb_cput(rtdb,'cpsd:exchange_correlation',
cc     >                    1,exchange_correlation(1:ind))
c
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:mult',mt_int,1,multiplicity)
c      value = value.and.
c     >        rtdb_put(rtdb,'cpsd:ispin',mt_int,1,ispin)
c
c*    ***** Error writing to RTDB *****
c      if (.not.value)
c     >  call errquit(
c     >       'cpsd_default: error writing to rtdb', 0, RTDB_ERR)
c
c
c      return
c      end
