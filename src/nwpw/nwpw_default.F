c
c $Id$
c
      subroutine nwpw_default(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

c
*     ***** local variables *****
      integer ind               
      logical      value
      character*50 cell_name,zone_name,structure_name
      character*50 input_wavefunction_filename
      character*50 output_wavefunction_filename
      character*50 input_ewavefunction_filename
      character*50 output_ewavefunction_filename

      character*50 exchange_correlation
      double precision fake_mass,time_step,rcut
      integer loop(2),npsp,ncut
      double precision tolerances(3),ecut,wcut
      logical frac_coord,mulliken,allow_translation
      integer ispin,multiplicity
      character*50 rtdb_name
      character*50 boundry
      double precision unita(3,3)
      integer ngrid(3),nx,ny,nz,ewald_ngrid(3)

       
*     **** don't set defaults if they already exist ****
      value = btdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut)
      if (value) return


*     ***** initializations ****
      cell_name             = 'cell_default'
      zone_name             = 'zone_default'
      structure_name        = 'structure_default'
      call util_file_prefix('movecs',input_wavefunction_filename)
      call util_file_prefix('movecs',output_wavefunction_filename)
      call util_file_prefix('emovecs',input_ewavefunction_filename)
      call util_file_prefix('emovecs',output_ewavefunction_filename)
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
      ncut = 1
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

c ---- commented out default handled in control read ----
c     if (cell_name.eq.' ') 
c    >  call errquit('nwpw_default: nead a cell_name', 0, INPUT_ERR)
c     ind = index(cell_name,' ') - 1
c     value = btdb_cput(rtdb,'cgsd:cell_name',1,cell_name(1:ind))
c     value = btdb_cput(rtdb,'band:cell_name',1,cell_name(1:ind))

      if (input_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'nwpw_default: nead an input wavefunction_filename', 0,
     &       INPUT_ERR)
      if (output_wavefunction_filename.eq.' ') 
     >  call errquit(
     >       'nwpw_default: nead an output wavefunction_filename', 0,
     &       INPUT_ERR)


      ind = index(zone_name,' ') - 1
      value = btdb_cput(rtdb,'band:zone_name',1,zone_name(1:ind))
      ind = index(structure_name,' ') - 1
      value = btdb_cput(rtdb,'band_structure:zone_name',1,
     >                  structure_name(1:ind))

      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'cgsd:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'band:input_wavefunction_filename',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(input_ewavefunction_filename,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'cgsd:input_ewavefunction_filename',
     >                    1,input_ewavefunction_filename(1:ind))

*     *********************************
*     **** used by task_save_state ****
*     *********************************
      ind = index(input_wavefunction_filename,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'pspw:input vectors',
     >                    1,input_wavefunction_filename(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'band:input vectors',
     >                    1,input_wavefunction_filename(1:ind))

      ind = index(output_wavefunction_filename,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'cgsd:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'band:output_wavefunction_filename',
     >                    1,output_wavefunction_filename(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'pspw:output vectors',
     >                    1,output_wavefunction_filename(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'band:output vectors',
     >                    1,output_wavefunction_filename(1:ind))

      ind = index(output_ewavefunction_filename,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'cgsd:output_ewavefunction_filename',
     >                    1,output_ewavefunction_filename(1:ind))

      value = value.and.
     >        btdb_put(rtdb,'cgsd:fractional_coordinates',mt_log,1,
     >                       frac_coord)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:mulliken',mt_log,1,
     >                       mulliken)
      value = value.and.
     >        btdb_put(rtdb,'band:mulliken',mt_log,1,
     >                       mulliken)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:allow_translation',mt_log,1,
     >                       allow_translation)
      value = value.and.
     >        btdb_put(rtdb,'band:allow_translation',mt_log,1,
     >                       allow_translation)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:npsp',     mt_int,1,npsp)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        btdb_put(rtdb,'band:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        btdb_put(rtdb,'band:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:loop',mt_int,2,loop)
      value = value.and.
     >        btdb_put(rtdb,'band:loop',mt_int,2,loop)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      value = value.and.
     >        btdb_put(rtdb,'band:tolerances',mt_dbl,3,tolerances)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        btdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        btdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        btdb_put(rtdb,'band:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:ncut',mt_int,1,ncut)
      value = value.and.
     >        btdb_put(rtdb,'band:ncut',mt_int,1,ncut)
      ind = index(exchange_correlation,' ') - 1
      value = value.and.
     >        btdb_cput(rtdb,'cgsd:exchange_correlation',
     >                    1,exchange_correlation(1:ind))
      value = value.and.
     >        btdb_cput(rtdb,'band:exchange_correlation',
     >                    1,exchange_correlation(1:ind))

      value = value.and.
     >        btdb_put(rtdb,'cgsd:mult',mt_int,1,multiplicity)
      value = value.and.
     >        btdb_put(rtdb,'band:mult',mt_int,1,multiplicity)
      value = value.and.
     >        btdb_put(rtdb,'cgsd:ispin',mt_int,1,ispin)
      value = value.and.
     >        btdb_put(rtdb,'band:ispin',mt_int,1,ispin)


*     **** unit cell stuff ****
c---commented out - defaults handled in control_read ---
c     ind = index(cell_name,' ') -1
c     rtdb_name = cell_name(1:ind)//':boundry'
c     ind = index(boundry,' ') - 1
c     value = value.and.
c    >        btdb_cput(rtdb,rtdb_name,1,boundry(1:ind))

c---commented out - defaults handled in control_read ---
c     ind = index(cell_name,' ') -1
c     rtdb_name = cell_name(1:ind)//':unita'
c     value = value.and.
c    >        btdb_put(rtdb,rtdb_name,mt_dbl,9,unita)

c---commented out - defaults handled in control_read ---
c     ind = index(cell_name,' ') -1
c     rtdb_name = cell_name(1:ind)//':ngrid'
c     value = value.and.
c    >        btdb_put(rtdb,rtdb_name,mt_int,3,ngrid)


*     **** unit brillouin_zone stuff ***
      nx = 1
      ny = 1
      nz = 1
      call monkhorst_pack_set(rtdb,zone_name,nx,ny,nz)
      call monkhorst_pack_set(rtdb,structure_name,nx,ny,nz)
      

*     ***** Error writing to RTDB *****
      if (.not.value)
     >  call errquit(
     >       'nwpw_default: error writing to rtdb', 0, RTDB_ERR)

      return
      end
