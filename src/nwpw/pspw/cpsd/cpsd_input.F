      subroutine cpsd_input(rtdb)
c
c $Id: cpsd_input.F,v 1.2 2001-10-26 21:05:08 bylaska Exp $
c
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
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
     >            'geometry_optimize',
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
     >            'end'/

       character*30 cell_name
       character*30 input_wavefunction_filename
       character*30 output_wavefunction_filename
       character*30 exchange_correlation
       logical geometry_optimize,frac_coord,mulliken
       double precision fake_mass,time_step,rcut
       integer loop(2),npsp,ncut
       double precision tolerances(3),ecut,wcut,fe
        
       


*     ***** initializations ****
      cell_name             = ' '

      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)
c     input_wavefunction_filename  = ' '
c     output_wavefunction_filename = ' '

      exchange_correlation         = 'vosko'
      geometry_optimize            = .false.
      frac_coord                   = .false.
      mulliken                     = .false.
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
     >           'cpsd_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cpsd_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cpsd_input: unknown directive', 0)


      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900,
     >       1000,1100,1200,1300,1400,1500,1600,
     >      9999) ind
      call errquit(
     >      'psp_formatter_input: unimplemented directive', ind)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cpsd_input: failed to read cell_name', 0)
      goto 10


c
c geometry_optimize
c
 200  geometry_optimize = .true.
      goto 10


c
c  formatted_psp_filename
c 
 300  if (.not. inp_a(test))
     >     call errquit(
     >          'cpsd_input: failed to read psp_filename', 0)
      npsp = npsp + 1
      id = 'cpsd:psp'//CHAR(npsp)
      ind = index(test,' ') - 1
         if (.not. rtdb_cput(rtdb,id,
     >                    1,test(1:ind)))
     >        call errquit(
     >        'cpsd_input: rtdb_cput failed', 0)
      goto 10


c
c  input_wavefunction_filename
c 
 400  if (.not. inp_a(input_wavefunction_filename))
     >     call errquit(
     >          'cpsd_input: failed to read input_wavefunction', 0)
      goto 10

c
c  output_wavefunction_filename
c 
 500  if (.not. inp_a(output_wavefunction_filename))
     >     call errquit(
     >          'cpsd_input: failed to read output_wavefunction', 0)
      goto 10


c
c  fake_mass
c 
 600  if (.not. inp_f(fake_mass))
     >     call errquit(
     >          'wvfnc_init_input: failed to read nup', 0)
      goto 10
 
c
c  time_step
c 
 700  if (.not. inp_f(time_step))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0)
      goto 10
 
 
c
c  loop
c 
 800  if (.not. inp_i(loop(1)))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0)
      if (.not. inp_i(loop(2)))
     >     call errquit(
     >          'wvfnc_init_input: failed to read ndown', 0)
      goto 10
 
 
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
     >          'cpsd_input: failed to read ndown', 0)
      goto 10

c
c  energy_cutoff
c 
 1100 if (.not. inp_f(wcut))
     >     call errquit(
     >          'cpsd_input: failed to read ndown', 0)
      goto 10

c
c  ewald_rcut
c 
 1200 if (.not. inp_f(rcut))
     >     call errquit(
     >          'cpsd_input: failed to read rcut', 0)
      goto 10

c
c  ewald_rcut
c 
 1300 if (.not. inp_i(ncut))
     >     call errquit(
     >          'cpsd_input: failed to read ncut', 0)
      goto 10
c
c  exchange_correlation
c 
 1400  if (.not. inp_a(exchange_correlation))
     >     call errquit(
     >          'cpsd_input: failed to read exchange_correlation', 0)
      goto 10

c
c fractional coordinates
c
 1500 frac_coord = .true.
      goto 10


c
c mullikent
c
 1600 mulliken = .true.
      goto 10


*     ***** add wavefunction to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('cpsd_input: nead a cell_name', 0)
      if (input_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cpsd_input: nead an input wavefunction_filename', 0)
      if (output_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'cpsd_input: nead an output wavefunction_filename', 0)


      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:cell_name',1,cell_name(1:ind))

      ind = index(input_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(output_wavefunction_filename,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))

      value = rtdb_put(rtdb,'cpsd:geometry_optimize',mt_log,1,
     >                       geometry_optimize)
      value = rtdb_put(rtdb,'cpsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      value = rtdb_put(rtdb,'cpsd:mulliken',mt_log,1,
     >                       mulliken)

      value = rtdb_put(rtdb,'cpsd:npsp',     mt_int,1,npsp)
      value = rtdb_put(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass)
      value = rtdb_put(rtdb,'cpsd:time_step',mt_dbl,1,time_step)
      value = rtdb_put(rtdb,'cpsd:loop',mt_int,2,loop)
      value = rtdb_put(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances)
      value = rtdb_put(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
      value = rtdb_put(rtdb,'cpsd:wcut',mt_dbl,1,wcut)
      value = rtdb_put(rtdb,'cpsd:rcut',mt_dbl,1,rcut)
      value = rtdb_put(rtdb,'cpsd:ncut',mt_int,1,ncut)
      ind = index(exchange_correlation,' ') - 1
      value = rtdb_cput(rtdb,'cpsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      return
      end
