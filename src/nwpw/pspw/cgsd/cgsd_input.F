c
c $Id$
c
      subroutine cgsd_input(rtdb)
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
      parameter (num_dirs = 20)
 
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
     >            'mult',
     >            'dft',
     >            'odft',
     >            'end'/

       character*50 cell_name
       character*50 input_wavefunction_filename
       character*50 output_wavefunction_filename
       character*50 exchange_correlation
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut
       double precision tolerances(3),ecut,wcut,fe
       logical frac_coord,mulliken,allow_translation
       integer ispin,multiplicity
       
*     ***** initializations     ****
*     ***** should be dead code ****
      npsp = 0
      call cgsd_input_default(rtdb)

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cgsd_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cgsd_input: failed to read keyword', 0, INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cgsd_input: unknown directive', 0, INPUT_ERR)


      goto ( 100, 300, 400, 500, 600, 700, 800, 900,
     >       1000,1100,1200,1300,1400,1500,1600,1700,
     >       1800,1900,2000,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cgsd_input: failed to read cell_name', 0, INPUT_ERR)
      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write cell_name',0,
     &       RTDB_ERR)
      goto 10


c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cgsd_input: failed to read psp_filename', 0, INPUT_ERR)
      npsp = npsp + 1
      id = 'cgsd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cgsd_input: rtdb_cput failed', 0, RTDB_ERR)
      value = rtdb_put(rtdb,'cgsd:npsp',     mt_int,1,npsp)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write npsp',0, RTDB_ERR)
      goto 10


c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cgsd_input: failed to read input_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))

*     *********************************
*     **** used by task_save_state ****
*     *********************************
      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.rtdb_cput(rtdb,'pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write input vectors',0,
     &       RTDB_ERR)
      goto 10

c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cgsd_input: failed to read output_wavefunction', 0,
     &       INPUT_ERR)
      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'pspw:output vectors',
     >                    1,output_wavefunction_filename(1:ind))
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write output vectors',0,
     &       RTDB_ERR)
      goto 10


c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'cgsd_input: failed to fake mass', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write fake_mass',0,
     &       RTDB_ERR)
      goto 10
 
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'cgsd_input: failed to time step', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write time_step',0,
     &       RTDB_ERR)
      goto 10
 
 
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'cgsd_input: failed to read loop(1)', 0, INPUT_ERR)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'cgsd_input: failed to read loop(2)', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write loop',0, RTDB_ERR)
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
 901  value = rtdb_put(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write tolerances',0,
     &       RTDB_ERR)
      goto 10


c
c  energy_cutoff
c 
 1000 if (.not. inp_f(ecut))
     >     call errquit(
     >          'cgsd_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write ecut',0, RTDB_ERR)
      goto 10

c
c  wavefunction_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cgsd_input: failed to read ndown', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write wcut',0, RTDB_ERR)
      goto 10

c
c  ewald_rcut
c 
 1200 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cgsd_input: failed to read rcut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write rcut',0, RTDB_ERR)
      goto 10

c
c  ewald_ncut
c 
 1300 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cgsd_input: failed to read ncut', 0, INPUT_ERR)
      value = rtdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write ncut',0, RTDB_ERR)
      goto 10
c
c  exchange_correlation
c 
 1400  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cgsd_input: failed to read exchange_correlation', 0,
     &       INPUT_ERR)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write xc',0, RTDB_ERR)
      goto 10
c
c  fractional_coordinates
c 
 1500 frac_coord = .true.
      value = rtdb_put(rtdb,'cgsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write frac',0, RTDB_ERR)
      goto 10

c
c  Mulliken
c 
 1600 mulliken = .true.
      value = rtdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write mulliken',0, RTDB_ERR)
      goto 10
c
c  allow_translation
c 
 1700 allow_translation = .true.
      value = rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,
     >                       allow_translation)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write translation',0,
     &       RTDB_ERR)
      goto 10
c
c  multiplicity
c 
 1800 if (.not. inp_i(multiplicity))
     >     call errquit(
     >          'cgsd_input: failed to read mult', 0, INPUT_ERR)
      if (multiplicity.lt.1) then
         multiplicity = 1
         ispin        = 1
      else if (multiplicity.gt.1) then
         ispin        = 2
      end if
      value = rtdb_put(rtdb,'cgsd:mult',mt_int,1,multiplicity)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write mult',0, RTDB_ERR)
      goto 10
c
c  dft
c 
 1900 ispin = 1
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write ispin',0, RTDB_ERR)
      goto 10
c
c  odft
c 
 2000 ispin = 2
      value = rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
      if (.not.value) 
     >  call errquit('cgsd_input: failed to write ispin',0, RTDB_ERR)
      goto 10

*     ***** add wavefunction to rtdb ****
 9999 continue

      return
      end


      subroutine cgsd_input_default(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "rtdb.fh"

c
*     ***** local variables *****
      integer ind               
      logical      value
      character*50 cell_name
      character*50 input_wavefunction_filename
      character*50 output_wavefunction_filename
      character*50 exchange_correlation
      double precision fake_mass,time_step,rcut
      integer loop(2),npsp,ncut
      double precision tolerances(3),ecut,wcut
      logical frac_coord,mulliken,allow_translation
      integer ispin,multiplicity
      character*50 rtdb_name
      character*50 boundry
      double precision unita(3,3)
      integer ngrid(3)

       
*     **** don't set defaults if they already exist ****
      value = rtdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut)
      if (value) return


*     ***** initializations ****
      cell_name             = 'cell_default'
      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)
      exchange_correlation         = 'vosko'
      frac_coord                   = .false.
      mulliken                     = .false.
      allow_translation            = .false.
      fake_mass = 400000.0d0
      time_step = 5.8d0
      loop(1) = 10
      loop(2) = 100
      tolerances(1) = 1.0d-7
      tolerances(2) = 1.0d-7
      tolerances(3) = 1.0d-4
      ecut=9000.0d0
      wcut=9000.0d0
      rcut = 0.0d0
      ncut = 0
      npsp = 0
      ispin        = 1
      multiplicity = 1
      ngrid(1) = 32
      ngrid(2) = 32
      ngrid(3) = 32
      unita(1,1) = 0.0d0
      unita(2,1) = 0.0d0
      unita(3,1) = 0.0d0
      unita(1,2) = 0.0d0
      unita(2,2) = 0.0d0
      unita(3,2) = 0.0d0
      unita(1,3) = 0.0d0
      unita(2,3) = 0.0d0
      unita(3,3) = 0.0d0
      boundry = 'periodic'


 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('cgsd_default: nead a cell_name', 0, INPUT_ERR)
      if (input_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cgsd_default: nead an input wavefunction_filename', 0,
     &       INPUT_ERR)
      if (output_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cgsd_default: nead an output wavefunction_filename', 0,
     &       INPUT_ERR)

      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))

      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cgsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

*     *********************************
*     **** used by task_save_state ****
*     *********************************
      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(output_wavefunction_filename,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cgsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))

      value = value.and.
     >        rtdb_put(rtdb,'cgsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:allow_translation',mt_log,1,
     >                       allow_translation)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:npsp',     mt_int,1,npsp)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
      ind = index(exchange_correlation,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      value = value.and.
     >        rtdb_put(rtdb,'cgsd:mult',mt_int,1,multiplicity)
      value = value.and.
     >        rtdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)


*     **** unit cell stuff ****
      ind = index(cell_name,' ') -1
      rtdb_name = cell_name(1:ind)//':boundry'
      ind = index(boundry,' ') - 1
      value = value.and.
     >        rtdb_cput(rtdb,rtdb_name,1,boundry(1:ind))

      ind = index(cell_name,' ') -1
      rtdb_name = cell_name(1:ind)//':unita'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_dbl,9,unita)

      ind = index(cell_name,' ') -1
      rtdb_name = cell_name(1:ind)//':ngrid'
      value = value.and.
     >        rtdb_put(rtdb,rtdb_name,mt_int,3,ngrid)
      

*     ***** Error writing to RTDB *****
      if (.not.value)
     >  call errquit(
     >       'cgsd_default: error writing to rtdb', 0, RTDB_ERR)

      return
      end
