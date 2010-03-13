*
* $Id: control.F,v 1.130 2009-03-26 00:18:52 bylaska Exp $
*

*     ***********************************
*     *					*
*     *		control_read		*
*     *					*
*     ***********************************
      logical function control_read(code_in,rtdb)
      implicit none 
      integer code_in
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "control.fh"

      logical value 
      integer ispin0,ne(2),nbrill,np1,np2

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose

*     **** control_rtdb common block ****
      integer trtdb
      common / control_rtdb1 / trtdb

*     **** control_print common block ****
      integer print_level
      common / control_print1 / print_level

      !character*30 cell_name
      character*30 rtdb_unita,rtdb_ngrid,rtdb_boundry
      character*30 exchange_correlation
      integer i,l,ind,np

      integer  control_num_kvectors
      external control_num_kvectors
      logical  control_allow_translation
      external control_allow_translation


      value = rtdb_parallel(.true.)
      code = code_in
      trtdb = rtdb

      call util_print_get_level(print_level)
      !write(*,*) "Print level is ",print_level

*     **** get parallel mappings ****
      if (.not.rtdb_get(rtdb,'nwpw:mapping',mt_int,1,mapping)) then
         mapping = 1
      end if

*     **** set mapping1d ****
      if (.not.rtdb_get(rtdb,'nwpw:mapping1d',mt_int,1,mapping1d)) then
         mapping1d = 1
      end if

*     **** set np_dimensions ****
      call Parallel_np(np)

**** if NOMPI then no processor groups or communicators ****
#ifdef NOMPI
      np_dimensions(1) = np
      np_dimensions(2) = 1
      np_dimensions(3) = 1
#else
      if(.not.rtdb_get(rtdb,'nwpw:np_dimensions',
     >                 mt_int,3,np_dimensions)) then
         np_dimensions(1) = np
         np_dimensions(2) = 1
         np_dimensions(3) = 1
      else
         if (.not.((code.eq.5).or.(code.eq.10))) np_dimensions(3) = 1
         if (np_dimensions(3).lt.1) np_dimensions(3) = 1
         if (np_dimensions(2).lt.1) np_dimensions(2) = 1
         if (np_dimensions(1).lt.1) np_dimensions(1) = 1

*        **** reset np_dimensions(3) if larger than nbrill ****
         nbrill = control_num_kvectors()
         if (np_dimensions(3).gt.nbrill) np_dimensions(3) = nbrill

*        **** reset np_dimensions(3) if it is not a  multiple of np ****
         do while ((mod(np,np_dimensions(3)).ne.0).and.
     >             (np_dimensions(3).gt.1))
            np_dimensions(3) = np_dimensions(3) - 1
         end do
*        **** reset np_dimensions(2) if it is not a  multiple of np2 ****
         np = np/np_dimensions(3)
         do while ((mod(np,np_dimensions(2)).ne.0).and.
     >             (np_dimensions(2).gt.1))
            np_dimensions(2) = np_dimensions(2) - 1
         end do

         !*** temporary restriction until ne parallelized in band ***
         if ((code.eq.5).or.(code.eq.10)) np_dimensions(2) = 1

         np_dimensions(1) = np/np_dimensions(2)

      endif
#endif

*     **** get balance mapping ****
      if (.not.rtdb_get(rtdb,'nwpw:balance',mt_log,1,balance)) then
         balance = .true.
      end if


*     *********************************
*     **** cpsd and band_sd: stuff ****
*     *********************************
      if ((code.eq.1).or.(code.eq.13)) then
      if (.not.rtdb_cget(rtdb,'cpsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      if (.not.rtdb_cget(rtdb,'cpsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cpsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     > call util_file_prefix('movecs',output_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cpsd:exchange_correlation',
     >                  1,exchange_correlation))
     >   exchange_correlation = 'vosko'


#include "control_gga.fh"


      if(.not.rtdb_get(rtdb,'cpsd:geometry_optimize',mt_log,1,move))
     >    move = .false.
      if(.not.rtdb_get(rtdb,'cpsd:fractional_coordinates',
     >          mt_log,1,frac_coord))
     >     frac_coord = .false.
      if(.not.rtdb_get(rtdb,'cpsd:npsp',mt_int,1,npsp))
     >    npsp=0
      if(.not.rtdb_get(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass))
     >     fake_mass = 400000.0d0
      if(.not.rtdb_get(rtdb,'cpsd:time_step',mt_dbl,1,time_step))
     >     time_step = 5.8d0
      if(.not.rtdb_get(rtdb,'cpsd:loop',mt_int,2,loop)) then
       loop(1)=10 
       loop(2)=100
      end if
      if(.not. rtdb_get(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances))
     >then
         tolerances(1) = 1.0d-9
         tolerances(2) = 1.0d-9
         tolerances(3) = 1.0d-4
      end if 
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      if(.not.rtdb_get(rtdb,'cpsd:ecut',mt_dbl,1,ecut))
     >    ecut=9000.0d0
      if(.not.rtdb_get(rtdb,'cpsd:wcut',mt_dbl,1,wcut))
     >  wcut=ecut
      if(.not.rtdb_get(rtdb,'cpsd:rcut',mt_dbl,1,rcut))
     >    rcut = 0.0d0
      if(.not.rtdb_get(rtdb,'cpsd:ncut',mt_int,1,ncut))
     >    ncut = 1
      if(.not.rtdb_get(rtdb,'cpsd:mult',mt_int,1,multiplicity))
     >    multiplicity = 1
      if(.not.rtdb_get(rtdb,'cpsd:ispin',mt_int,1,ispin))
     >    ispin=1

      if (.not.rtdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation))
     >  rotation = .true.

      if (.not.rtdb_get(rtdb,'nwpw:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2


*     ****************************************
*     **** cpmd and band_cpmd code: stuff ****
*     ****************************************
      else if ((code.eq.2) .or. (code.eq.14)) then
      if (.not.rtdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      if(.not.rtdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if(.not.rtdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     >  call util_file_prefix('movecs',output_wavefunction_filename)
      if(.not.rtdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',input_v_wavefunction_filename)
      if(.not.rtdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',output_v_wavefunction_filename)
      if(.not.rtdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)
      ind = index(exchange_correlation,' ') - 1
      if(.not.rtdb_cget(rtdb,'cpmd:exchange_correlation',
     >                  1,exchange_correlation))
     >  exchange_correlation = 'vosko'

#include "control_gga.fh"


      if(.not. rtdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move))
     >    move = .true.
      if(.not. rtdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >   frac_coord = .false.
      if(.not. rtdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp))
     >    npsp = 0
      if(.not. rtdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass))
     >    fake_mass = 800.0d0
      if(.not. rtdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step))
     >    time_step = 5.0d0
      if(.not. rtdb_get(rtdb,'cpmd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not. rtdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)) then
        scaling(1) = 1.0d0
        scaling(2) = 1.0d0
      end if
      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      if(.not. rtdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut))
     >    ecut=9000.0d0
      if(.not. rtdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut))
     >    wcut = ecut
      if(.not. rtdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut))
     >    rcut = 0.0d0
      if(.not. rtdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut))
     >    ncut = 1
      SA = .true.
      if (.not.rtdb_get(rtdb,'cpmd:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.rtdb_get(rtdb,'cpmd:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.rtdb_get(rtdb,'cpmd:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.

      if (.not.rtdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation))
     >  rotation = .true.


*     **** get thermostat information ****
      if (.not.rtdb_get(rtdb,'cpmd:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.rtdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe))
     >   Pe = 1000.0d0
      if (.not.rtdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.rtdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr))
     >   Pr = 1000.0d0
      if (.not.rtdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0

      if (.not.rtdb_get(rtdb,'nwpw:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2


*     ************************************
*     **** cgsd: stuff or  paw: stuff ****
*     ************************************
      else if ((code.eq.3).or.(code.eq.8).or.
     >         (code.eq.11).or.(code.eq.12)) then
      if (.not.rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
c
c     **** Figure input/output MO vectors ****
c
      if (.not. rtdb_cget(rtdb, 'pspw:input vectors', 
     >                    1,input_wavefunction_filename))
     >  input_wavefunction_filename = 'atomic'
      if (.not. rtdb_cget(rtdb, 'pspw:output vectors', 
     >                    1,output_wavefunction_filename))
     >     output_wavefunction_filename = ' '
      if (output_wavefunction_filename.eq.' ')then
         if (input_wavefunction_filename.eq.'atomic')then
           call util_file_prefix('movecs',output_wavefunction_filename)
         else
            output_wavefunction_filename = input_wavefunction_filename
         endif
      endif
      if (input_wavefunction_filename.eq.'atomic')then
         input_wavefunction_filename = output_wavefunction_filename
      end if
      if (.not.rtdb_cget(rtdb,'cgsd:input_ewavefunction_filename',
     >                  1,input_ewavefunction_filename))
     >  call util_file_prefix('emovecs',input_ewavefunction_filename)

      if (.not.rtdb_cget(rtdb,'cgsd:output_ewavefunction_filename',
     >                  1,output_ewavefunction_filename))
     >  call util_file_prefix('emovecs',output_ewavefunction_filename)

      if (code.eq.11) then
         if(.not.rtdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename))
     >    call util_file_prefix('vmovecs',input_v_wavefunction_filename)
         if(.not.rtdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename))
     >   call util_file_prefix('vmovecs',output_v_wavefunction_filename)
      end if



      if (.not.rtdb_cget(rtdb,'cgsd:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'


#include "control_gga.fh"

*     ***** motion filenames ****
      if(.not.rtdb_cget(rtdb,'nwpw:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)


*     **** set Kohn-Sham scf parameters ***
      if (.not. rtdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not.rtdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 3

      if (.not.rtdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   ks_algorithm = 0

*     **** set maxit_orb maxit_orbs ***
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 1


      if (.not.rtdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp))
     >   npsp = 0
      if (.not.rtdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 400000.0d0
      if (.not.rtdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step))
     >  time_step = 5.8d0
      if (.not.rtdb_get(rtdb,'cgsd:loop',mt_int,2,loop)) then
        loop(1) = 10
        loop(2) = 100
      end if
      if (.not.rtdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances))then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if

      if(.not. rtdb_get(rtdb,'nwpw:scaling',mt_dbl,2,scaling)) then
        scaling(1) = 1.0d0
        scaling(2) = 1.0d0
      end if

      if (.not.rtdb_get(rtdb,'cgsd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >    frac_coord = .false.
  
      if (.not.rtdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut))
     >   ecut=9000.0d0
      if (.not.rtdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut))
     >   wcut = ecut
      if (.not.rtdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut))
     >   rcut = 0.0d0
      if (.not.rtdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut))
     >   ncut = 1
      if (.not.rtdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity))
     >   multiplicity = 1
      if (.not.rtdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin))
     >   ispin = 1


*     **** BO parameterss ***
      if (.not.rtdb_get(rtdb,'nwpw:bo_steps',mt_int,2,bo_steps)) then
         bo_steps(1) = 10
         bo_steps(2) = 100
      end if
      if (.not.rtdb_get(rtdb,'nwpw:bo_time_step',mt_dbl,1,bo_time_step))
     >   bo_time_step = time_step
      if (.not.rtdb_get(rtdb,'nwpw:bo_algorithm',mt_int,1,bo_algorithm))
     >   bo_algorithm = 0
      if (.not.rtdb_get(rtdb,'nwpw:bo_fake_mass',mt_dbl,1,bo_fake_mass))
     >   bo_fake_mass = 500.0d0

      if (.not.rtdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation))
     >  rotation = .true.

      SA = .true.
      if (.not.rtdb_get(rtdb,'nwpw:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.rtdb_get(rtdb,'nwpw:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.rtdb_get(rtdb,'nwpw:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.



*     **** get thermostat information ****
      if (.not.rtdb_get(rtdb,'nwpw:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.rtdb_get(rtdb,'nwpw:Pe',mt_dbl,1,Pe))
     >   Pe = 1000.0d0
      if (.not.rtdb_get(rtdb,'nwpw:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.rtdb_get(rtdb,'nwpw:Pr',mt_dbl,1,Pr))
     >   Pr = 1000.0d0
      if (.not.rtdb_get(rtdb,'nwpw:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0


*     ***************************
*     **** pspw_dplot: stuff ****
*     ***************************
      else if (code.eq.4) then
         value = .true.
         if (.not.rtdb_cget(rtdb,'pspw_dplot:wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >     call util_file_prefix('movecs',input_wavefunction_filename)
         call psi_get_header(i,ngrid,unita,ispin0,ne)
         if (i.eq.3) boundry = 'periodic'
         if (i.eq.4) boundry = 'aperiodic'

*        **** dummy variables ****
         move       = .false.
         frac_coord = .false.
         gga = 0
         fake_mass = 400000.0d0
         time_step = 5.8d0
         loop(1) = 0
         loop(2) = 0
         tolerances(1) = 1.0d-9
         tolerances(2) = 1.0d-9
         tolerances(3) = 1.0d-4
         ecut=9000.0d0
         wcut=9000.0d0
         rcut = 0.0d0
         ncut = 1
         npsp = 0

         control_read = value
         return

*     *********************
*     **** band: stuff ****
*     *********************
      else if (code.eq.5) then
     
      if (.not.rtdb_cget(rtdb,'band:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

c     **** Figure input/output MO vectors ****
c
      if (.not. rtdb_cget(rtdb, 'pspw:input vectors', 
     >                    1,input_wavefunction_filename))
     >  input_wavefunction_filename = 'atomic'

      if (.not. rtdb_cget(rtdb, 'pspw:output vectors', 
     >                    1,output_wavefunction_filename))
     >     output_wavefunction_filename = ' '
      if (output_wavefunction_filename.eq.' ')then
         if (input_wavefunction_filename.eq.'atomic')then
           call util_file_prefix('movecs',output_wavefunction_filename)
         else
            output_wavefunction_filename = input_wavefunction_filename
         endif
      endif
      if (input_wavefunction_filename.eq.'atomic')then
         input_wavefunction_filename = output_wavefunction_filename
      end if
      if (.not.rtdb_cget(rtdb,'cgsd:input_ewavefunction_filename',
     >                  1,input_ewavefunction_filename))
     >  call util_file_prefix('emovecs',input_ewavefunction_filename)

      if (.not.rtdb_cget(rtdb,'cgsd:output_ewavefunction_filename',
     >                  1,output_ewavefunction_filename))
     >  call util_file_prefix('emovecs',output_ewavefunction_filename)


      if (.not.rtdb_cget(rtdb,'band:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

      if (inp_compare(.false.,exchange_correlation,'vosko')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
         gga = 10
      else if (inp_compare(.false.,exchange_correlation,'blyp')) then
         gga = 11
      else if (inp_compare(.false.,exchange_correlation,'revpbe')) then
         gga = 12


*     **** Hartree-Fock ***
      else if (inp_compare(.false.,exchange_correlation,'hf')) then
        gga = 200
        value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'hartree-fock')) then
        gga = 200
        value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)


*     **** hybrid functionals ****
      else if (inp_compare(.false.,exchange_correlation,'pbe0')) then
        gga = 110
        value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'band:HFX_parameter',mt_dbl,1,0.25d0)

      else if (inp_compare(.false.,exchange_correlation,'revpbe0')) then
        gga = 112
        value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'band:HFX_parameter',mt_dbl,1,0.25d0)

      else if (inp_compare(.false.,exchange_correlation,'bnl')) then
        gga = 113
        value = value.and.rtdb_put(rtdb,'band:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'band:HFX_parameter',mt_dbl,1,1.00d0)
        value = value.and.
     >          rtdb_put(rtdb,'band:HFX_screening_type',mt_int,1,2)


      else
         gga = 0
      end if

      if(.not.rtdb_get(rtdb,'band:geometry_optimize',mt_log,1,move))
     >    move = .false.
      if (.not.rtdb_get(rtdb,'band:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 400000.0d0
      if (.not.rtdb_get(rtdb,'band:time_step',mt_dbl,1,time_step))
     >   time_step = 5.8d0
      if (.not.rtdb_get(rtdb,'band:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not.rtdb_get(rtdb,'band:tolerances',mt_dbl,3,tolerances)) then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      if (.not.rtdb_get(rtdb,'band:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.rtdb_get(rtdb,'band:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.rtdb_get(rtdb,'band:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.rtdb_get(rtdb,'band:ncut',mt_int,1,ncut))
     >  ncut = 1
      if (.not.rtdb_get(rtdb,'band:mult',mt_int,1,multiplicity))
     >   multiplicity = 1
      if (.not.rtdb_get(rtdb,'band:ispin',mt_int,1,ispin))
     >   ispin = 1

      if (.not.rtdb_get(rtdb,'band:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2
      


*     **** set Kohn-Sham scf parameters ***
      if (.not. rtdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not.rtdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 3
      if (.not.rtdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   ks_algorithm = 0

*     **** set maxit_orb maxit_orbs ***
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 1


*     ***********************
*     **** paw_sd: stuff ****
*     ***********************
      else if (code.eq.6) then

      if (.not.rtdb_get(rtdb,'cgsd:geometry_optimize',mt_log,1,move))
     >   move = .false.

      if (.not.rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      if (.not.rtdb_cget(rtdb,'cgsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cgsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     > call util_file_prefix('movecs',output_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cgsd:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

      if (inp_compare(.false.,exchange_correlation,'vosko')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
         gga = 10
      else if (inp_compare(.false.,exchange_correlation,'blyp')) then
         gga = 11
      else if (inp_compare(.false.,exchange_correlation,'revpbe')) then
         gga = 12
      else
         gga = 0
      end if

      if (.not.rtdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp))
     >   npsp = 0
      if (.not.rtdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass))
     >   fake_mass = 400000.0d0
      if (.not.rtdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step))
     >   time_step = 5.8d0
      if (.not.rtdb_get(rtdb,'cgsd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not.rtdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)) then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0
      if (.not.rtdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.rtdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.rtdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.rtdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut))
     >  ncut = 1
      if (.not.rtdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity))
     >  multiplicity = 1
      if (.not.rtdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin))
     >  ispin = 1


*     *************************
*     **** paw_cpmd: stuff ****
*     *************************
      else if (code.eq.7) then

      if (.not.rtdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

      if (.not.rtdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     >  call util_file_prefix('movecs',output_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',input_v_wavefunction_filename)
      if (.not.rtdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',output_v_wavefunction_filename)


*     ***** motion filenames ****
      if(.not.rtdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)


      if (.not.rtdb_cget(rtdb,'cpmd:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

      if (inp_compare(.false.,exchange_correlation,'vosko')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
         gga = 0
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
         gga = 10
      else if (inp_compare(.false.,exchange_correlation,'blyp')) then
         gga = 11
      else if (inp_compare(.false.,exchange_correlation,'revpbe')) then
         gga = 12
      else
         gga = 0
      end if


      if (.not.rtdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move))
     >   move = .false.
      if (.not.rtdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >   frac_coord = .false.
      if (.not.rtdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp))
     >  npsp = 0
      if (.not.rtdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 800.0d0
      if (.not.rtdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step))
     >  time_step = 5.0d0
      if (.not.rtdb_get(rtdb,'cpmd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if (.not.rtdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)) then
         scaling(1) = 1.0d0
         scaling(2) = 1.0d0
      end if
      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      if (.not.rtdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.rtdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.rtdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.rtdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut))
     >   ncut = 1

      SA = .true.
      if (.not.rtdb_get(rtdb,'cpmd:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.rtdb_get(rtdb,'cpmd:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.rtdb_get(rtdb,'cpmd:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.

      if (.not.rtdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation))
     >  rotation = .true.

*     **** get thermostat information ****
      if (.not.rtdb_get(rtdb,'cpmd:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.rtdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe))
     >   Pe = 1000.0d0
      if (.not.rtdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.rtdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr))
     >   Pr = 1000.0d0
      if (.not.rtdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0


*     **********************
*     **** pspw_wannier ****
*     **********************
      else if (code.eq.9) then

c     **** Figure input/output MO vectors ****
c
      if (.not.rtdb_cget(rtdb,'wannier:input vectors',
     >                  1,input_wavefunction_filename)) then
      if (.not. rtdb_cget(rtdb, 'pspw:input vectors',
     >                    1,input_wavefunction_filename))
     >  input_wavefunction_filename = 'atomic'
      end if

      if (.not.rtdb_cget(rtdb,'wannier:output vectors',
     >                  1,output_wavefunction_filename)) then
      if (.not. rtdb_cget(rtdb, 'pspw:output vectors',
     >                    1,output_wavefunction_filename))
     >     output_wavefunction_filename = ' '

      if (output_wavefunction_filename.eq.' ')then
         if (input_wavefunction_filename.eq.'atomic')then
            call util_file_prefix('movecs',output_wavefunction_filename)
         else
            output_wavefunction_filename = input_wavefunction_filename
         endif
      endif
      endif

      if (input_wavefunction_filename.eq.'atomic')then
         input_wavefunction_filename = output_wavefunction_filename
      end if

         call psi_get_header(i,ngrid,unita,ispin0,ne)
         if (i.eq.3) boundry = 'periodic'
         if (i.eq.4) boundry = 'aperiodic'

*        **** dummy variables ****
         move       = .false.
         frac_coord = .false.
         gga = 0
         fake_mass = 400000.0d0
         time_step = 5.8d0
         loop(1) = 0
         loop(2) = 0
         tolerances(1) = 1.0d-9
         tolerances(2) = 1.0d-9
         tolerances(3) = 1.0d-4
         if(.not.rtdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut)) ecut=9000.0d0
         if(.not.rtdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut)) wcut=ecut
         rcut = 0.0d0
         ncut = 0
         npsp = 0

         control_read = value
         return

*     ***************************
*     **** band_dplot: stuff ****
*     ***************************
      else if (code.eq.10) then
         value = .true.
         if (.not.rtdb_cget(rtdb,'band_dplot:wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >     call util_file_prefix('movecs',input_wavefunction_filename)
         call psi_get_header(i,ngrid,unita,ispin0,ne,nbrill)
         if (i.eq.5) boundry = 'periodic'

*        **** dummy variables ****
         move       = .false.
         frac_coord = .false.
         gga = 0
         fake_mass = 400000.0d0
         time_step = 5.8d0
         loop(1) = 0
         loop(2) = 0
         tolerances(1) = 1.0d-9
         tolerances(2) = 1.0d-9
         tolerances(3) = 1.0d-4
         ecut=9000.0d0
         wcut=9000.0d0
         rcut = 0.0d0
         ncut = 0
         npsp = 0

         control_read = value
         return

*     ***************************
*     **** unknown code type ****
*     ***************************
      else
         value = .false.
         write(*,*) "control_read: unknown code type:",code
         control_read = value
         return
      end if

*     *****************************
*     ***** symmetry variables ****
*     *****************************
      if (.not.rtdb_get(rtdb,'nwpw:symmetry',mt_int,1,symm_number))
     >   symm_number = 0

*     **********************
*     ***** cell: stuff ****
*     **********************
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      rtdb_ngrid = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'

           
*     **** define unita and boundary ****
      if (.not.rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)) then
        call dcopy(9,0.0d0,0,unita,1)
      end if

      if (.not.rtdb_cget(rtdb,rtdb_boundry,1,boundry)) then
         boundry = 'periodic'
      end if
      call check_unita_for_default(rtdb,unita,rtdb_unita,cell_name)


*     **** define ngrid based on unita and ecut ****
      if (.not.rtdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
         if ((ecut.gt.5000.0d0).and.(wcut.gt.5000.0d0)) then
            call control_ecut_wcut_default(rtdb,ecut,wcut)
         else if ((ecut.gt.5000.0d0).and.(wcut.lt.5000.0d0)) then
            ecut = 2.0d0*wcut
         end if
            
         call control_ngrid_default(rtdb,unita,ecut,mapping,ngrid)
         !ngrid(1) = 32
         !ngrid(2) = 32
         !ngrid(3) = 32
      end if


*     **** set fractional (smearing) parameters ****
      if (.not.
     > rtdb_get(rtdb,'nwpw:fractional_orbitals',mt_int,2,frac_ne)) then
         frac_ne(1) = 0
         frac_ne(2) = 0
      end if
      fractional = (frac_ne(1).gt.0).or.(frac_ne(2).gt.0)
      if (.not.rtdb_get(rtdb,'nwpw:fractional_temperature',
     >                  mt_dbl,1,frac_temperature)) then
         frac_temperature = 0.0d0
      end if
      if (.not.rtdb_get(rtdb,'nwpw:fractional_smeartype',
     >                  mt_int,1,frac_smeartype)) then
         frac_smeartype = 0
      end if

*     **** set attenuation parameter ****
      if (.not.rtdb_get(rtdb,'nwpw:attenuation',mt_dbl,1,attenuation))
     >  attenuation = 0.5d0

*     **** set out of time variables ****
      est_step_time   = -1
      est_finish_time = -1
      call current_second(cpu1_time)

*     **** set gram_schmidt ***
      gram_schmidt = .false.
      if (.not.rtdb_get(rtdb,
     >      'nwpw:gram-schmidt',mt_log,1,gram_schmidt))
     >  gram_schmidt = .false.

*     **** set translation ***
      translation = control_allow_translation()
      
*     **** set two component pseudopotential
      two_comp_ppot = .false.      


*     ****  ewald ngrid ************
cccc set to some reasonable default !
      if (.not.rtdb_get(rtdb,'nwpw:ewald_ngrid',mt_int,3,ewald_grid)) 
     > then
     	ewald_grid(1)=ngrid(1)
        ewald_grid(2)=ngrid(2)
        ewald_grid(3)=ngrid(3)
      end if

*     ****  reset mapping if slab and it doesnot fit ****
      if (mapping.eq.1) then
         if ((ngrid(2).ne.ngrid(3)).or.
     >       (ngrid(3).lt.np_dimensions(1))) then
            mapping = 2
         end if
      end if

      control_read = value
      return
      end
*     ***********************************
*     *  control_ewald_ngrid()
*     ***********************************
      integer function control_ewald_ngrid(i)
      implicit none
#include "control.fh"
      integer i
      control_ewald_ngrid=ewald_grid(i)
      return
      end
*     ***********************************
*     *  control_ewald_set_ngrid
*     ***********************************
      subroutine ewald_set_ngrid(enx,eny,enz)
      implicit none
      integer enx,eny,enz

#include "control.fh"

      ewald_grid(1)=enx
      ewald_grid(2)=eny
      ewald_grid(3)=enz
      return
      end

*     ***********************************
*     *                                 *
*     *    control_ngrid_default        *
*     *                                 *
*     ***********************************
      subroutine control_ngrid_default(rtdb,unita,ecut,mapping,ngrid)
      implicit none
      integer rtdb
      real*8 unita(3,3),ecut
      integer mapping
      integer ngrid(3)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** local variables ****
      real*8 unitg(3,3),omega
      real*8 gx,gy,gz
      real*8 xh,yh,zh
      integer  control_set_ngrid
      external control_set_ngrid

      call get_cube(unita,unitg,omega)

      gx = unitg(1,1)
      gy = unitg(2,1)
      gz = unitg(3,1)
      xh = dsqrt(2.0d0*ecut/(gx*gx + gy*gy + gz*gz))+0.5d0

      gx = unitg(1,2)
      gy = unitg(2,2)
      gz = unitg(3,2)
      yh = dsqrt(2.0d0*ecut/(gx*gx + gy*gy + gz*gz))+0.5d0

      gx = unitg(1,3)
      gy = unitg(2,3)
      gz = unitg(3,3)
      zh = dsqrt(2.0d0*ecut/(gx*gx + gy*gy + gz*gz))+0.5d0

      if (mapping.eq.2) then
c        ngrid(1) = control_set_ngrid(2.0d0*xh,.true.)
c        ngrid(2) = control_set_ngrid(2.0d0*yh,.false.)
c        ngrid(3) = control_set_ngrid(2.0d0*zh,.false.)
        ngrid(1) = control_set_ngrid(2.0d0*xh,.true.)
        ngrid(2) = control_set_ngrid(2.0d0*yh,.true.)
        ngrid(3) = control_set_ngrid(2.0d0*zh,.true.)
      else
        ngrid(1) = control_set_ngrid(2.0d0*xh,.true.)
        ngrid(2) = control_set_ngrid(2.0d0*yh,.true.)
        ngrid(3) = control_set_ngrid(2.0d0*zh,.true.)
        if (ngrid(2).gt.ngrid(3)) then
          ngrid(3) = ngrid(2)
        else
          ngrid(2) = ngrid(3)
        end if
      end if


c*     *** write unita to rtdb  - should happen only once during a simulation ***
c      if (.not.rtdb_put(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
c        call errquit('cannot write ngrid to rtdb',0,0)
c      end if

      return
      end

*     ***********************************
*     *                                 *
*     *    control_set_ngrid            *
*     *                                 *
*     ***********************************
*
*     return n so that it is a multiple of 2,3,5,7 
*
      integer function control_set_ngrid(x,mult2)
      implicit none
      real*8 x
      logical mult2

*     **** local variables ****
      integer nx,ntest
      integer nf2,nf3,nf5,nf7

      integer  factor_count2
      external factor_count2

*     **** find prime factors of nx ***
      nx = (x+0.5d0)  !*** crude rounding
      if ((mult2).and.(mod(nx,2).ne.0)) nx = nx+1

      nf2 = factor_count2(nx,2)
      nf3 = factor_count2(nx,3)
      nf5 = factor_count2(nx,5)
      nf7 = factor_count2(nx,7)
      ntest = (2**nf2) * (3**nf3) * (5**nf5) * (7**nf7)
      do while  (nx .ne. ntest)
        nx = nx + 1
        if (mult2) nx = nx + 1
        nf2 = factor_count2(nx,2)
        nf3 = factor_count2(nx,3)
        nf5 = factor_count2(nx,5)
        nf7 = factor_count2(nx,7)
        ntest = (2**nf2) * (3**nf3) * (5**nf5) * (7**nf7)
      end do

      control_set_ngrid = nx
      return 
      end

      integer function factor_count2(n,m)
      implicit none
      integer n,m
      integer f,nn

      f  = 0
      nn = n

      do while (mod(nn,m).eq.0)
        nn = nn/m
        f  = f + 1
      end do

      factor_count2 = f
      return
      end





*     ***********************************
*     *					*
*     *    check_unita_for_default	*
*     *					*
*     ***********************************
      subroutine check_unita_for_default(rtdb,unita,rtdb_unita,
     >                                   cell_name)
      implicit none
      integer rtdb
      real*8 unita(3,3)
      character*30 rtdb_unita
      character*30 cell_name
      


#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,box_orient
      integer geom,box_type,l
      real*8  box_delta
      character*30 rtdb_name

      value = (unita(1,1) .eq. 0.0d0).and.
     >        (unita(2,1) .eq. 0.0d0).and.
     >        (unita(3,1) .eq. 0.0d0).and.
     >        (unita(1,2) .eq. 0.0d0).and.
     >        (unita(2,2) .eq. 0.0d0).and.
     >        (unita(3,2) .eq. 0.0d0).and.
     >        (unita(1,3) .eq. 0.0d0).and.
     >        (unita(2,3) .eq. 0.0d0).and.
     >        (unita(3,3) .eq. 0.0d0)

      if (value) then
         value = geom_create(geom,'geometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
         value = value.and.geom_amatrix_get(geom,unita)
         if (.not. value) call errquit('cannot load geometry',0,
     &       GEOM_ERR)

         value = (unita(1,1) .eq. 1.0d0).and.
     >           (unita(2,1) .eq. 0.0d0).and.
     >           (unita(3,1) .eq. 0.0d0).and.
     >           (unita(1,2) .eq. 0.0d0).and.
     >           (unita(2,2) .eq. 1.0d0).and.
     >           (unita(3,2) .eq. 0.0d0).and.
     >           (unita(1,3) .eq. 0.0d0).and.
     >           (unita(2,3) .eq. 0.0d0).and.
     >           (unita(3,3) .eq. 1.0d0)
         if (value) then
           l = index(cell_name,' ') - 1
           rtdb_name = cell_name(1:l)//':box_delta'
           if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,1,box_delta))
     >       box_delta = 5.0d0
           rtdb_name = cell_name(1:l)//':box_type'
           if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,box_type)) 
     >       box_type = 0
           rtdb_name = cell_name(1:l)//':box_orient'
           if (.not.rtdb_get(rtdb,rtdb_name,mt_log,1,
     >                       box_orient)) 
     >       box_orient = .false.
           call control_find_box(geom,box_type,box_orient,box_delta,
     >                           unita)

*          *** write unita to rtdb  - should happen only once during a simulation ***
           if (.not.rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita)) then
             call errquit('cannot write unita to rtdb',0,0)
           end if

           !unita(1,1) = 20.0d0
           !unita(2,1) =  0.0d0
           !unita(3,1) =  0.0d0
           !unita(1,2) =  0.0d0
           !unita(2,2) = 20.0d0
           !unita(3,2) =  0.0d0
           !unita(1,3) =  0.0d0
           !unita(2,3) =  0.0d0
           !unita(3,3) = 20.0d0
         end if
         value = geom_destroy(geom)
         if (.not. value) call errquit('cannot destroy geom',0,
     &       GEOM_ERR)

      end if

      return
      end 

*     ***********************************
*     *					*
*     *	      control_find_box		*
*     *					*
*     ***********************************
      subroutine control_find_box(geom,box_type,orient,delta,unita)
      implicit none
      integer geom
      integer box_type
      logical orient
      real*8 delta
      real*8 unita(3,3)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     *** local variables ***
      integer lwork
      parameter (lwork=9)
      double precision work(lwork)

      integer ii,nion,ierr
      double precision q,rxyz(3),mtensor(3,3),frac(3),L1,L2,L3
      double precision x,y,z,meig(3)
      character*16     t
      logical value

*     **** external functions ****
      logical   control_notqmmmq
      external  control_notqmmmq

      if (.not.geom_ncent(geom,nion)) then
        call errquit('cannot load nion from geom',0,GEOM_ERR)
      end if

*     *** find principle axes wrt origin (0,0,0) ***
       mtensor(1,1) = 1.0d0
       mtensor(2,1) = 0.0d0
       mtensor(3,1) = 0.0d0
       mtensor(1,2) = 0.0d0
       mtensor(2,2) = 1.0d0
       mtensor(3,2) = 0.0d0
       mtensor(1,3) = 0.0d0
       mtensor(2,3) = 0.0d0
       mtensor(3,3) = 1.0d0
       if (orient) then
         call dcopy(9,0.0d0,0,mtensor,1)
         do ii=1,nion
           value = geom_cent_get(geom,ii,t,rxyz,q)
           if (.not.control_notqmmmq(t)) then
            x =rxyz(1)
            y =rxyz(2)
            z =rxyz(3)
            mtensor(1,1)=mtensor(1,1)+(y*y+z*z)
            mtensor(2,1)=mtensor(2,1)- x*y
            mtensor(1,2)=mtensor(2,1)
            mtensor(3,1)=mtensor(3,1)- x*z
            mtensor(1,3)=mtensor(3,1)
            mtensor(2,2)=mtensor(2,2)+(x*x+z*z)
            mtensor(3,2)=mtensor(3,2)- y*z
            mtensor(2,3)=mtensor(3,2)
            mtensor(3,3)=mtensor(3,3)+(x*x+y*y)
           end if
         end do
 
c        **** longest dimension is along a1 ****
         call DSYEV('V','U',3,mtensor,3,meig,work,lwork,ierr)

c         !*** reorder eigenvectors - make longest dimesion along a3 ****
c         x = mtensor(1,1)
c         y = mtensor(2,1)
c         z = mtensor(3,1)
c         mtensor(1,1) = mtensor(1,3)
c         mtensor(2,1) = mtensor(2,3)
c         mtensor(3,1) = mtensor(3,3)
c         mtensor(1,3) = x
c         mtensor(2,3) = y
c         mtensor(3,3) = z
       end if

*     *****************
*     *** cubic box ***
*     *****************
      if (box_type.eq.0) then

*     **** define L1 ***
      L1 = 0.0d0
      do ii=1,nion
        if (.not.geom_cent_get(geom,ii,t,rxyz,q)) then
          call errquit('cannot load center from geom',0,GEOM_ERR)
        end if
        if (.not.control_notqmmmq(t)) then
           frac(1) = mtensor(1,1)*rxyz(1)
     >             + mtensor(2,1)*rxyz(2)
     >             + mtensor(3,1)*rxyz(3)
           frac(2) = mtensor(1,2)*rxyz(1)
     >             + mtensor(2,2)*rxyz(2)
     >             + mtensor(3,2)*rxyz(3)
           frac(3) = mtensor(1,3)*rxyz(1)
     >             + mtensor(2,3)*rxyz(2)
     >             + mtensor(3,3)*rxyz(3)
           if ((frac(1)+delta).gt.( 0.5*L1)) L1 = 2.0d0*(frac(1)+delta)
           if ((frac(2)+delta).gt.( 0.5*L1)) L1 = 2.0d0*(frac(2)+delta)
           if ((frac(3)+delta).gt.( 0.5*L1)) L1 = 2.0d0*(frac(3)+delta)
           if ((frac(1)-delta).lt.(-0.5*L1)) L1 = 2.0d0*(frac(1)+delta)
           if ((frac(2)-delta).lt.(-0.5*L1)) L1 = 2.0d0*(frac(2)+delta)
           if ((frac(3)-delta).lt.(-0.5*L1)) L1 = 2.0d0*(frac(3)+delta)
        end if
      end do

*     **** define unit cell ****
      unita(1,1) = L1*mtensor(1,1)
      unita(2,1) = L1*mtensor(2,1)
      unita(3,1) = L1*mtensor(3,1)
      unita(1,2) = L1*mtensor(1,2)
      unita(2,2) = L1*mtensor(2,2)
      unita(3,2) = L1*mtensor(3,2)
      unita(1,3) = L1*mtensor(1,3)
      unita(2,3) = L1*mtensor(2,3)
      unita(3,3) = L1*mtensor(3,3)

*     ************************
*     *** orthorhombic box ***
*     ************************
      else if (box_type.eq.1) then

*     **** define L1, L2, and L3 ****
      L1 = 0.0d0
      L2 = 0.0d0
      L3 = 0.0d0
      do ii=1,nion
        if (.not.geom_cent_get(geom,ii,t,rxyz,q)) then
          call errquit('cannot load center from geom',0,GEOM_ERR)
        end if
        if (.not.control_notqmmmq(t)) then
           frac(1) = mtensor(1,1)*rxyz(1)
     >             + mtensor(2,1)*rxyz(2)
     >             + mtensor(3,1)*rxyz(3)
           frac(2) = mtensor(1,2)*rxyz(1)
     >             + mtensor(2,2)*rxyz(2)
     >             + mtensor(3,2)*rxyz(3)
           frac(3) = mtensor(1,3)*rxyz(1)
     >             + mtensor(2,3)*rxyz(2)
     >             + mtensor(3,3)*rxyz(3)
           if ((frac(1)+delta).gt.( 0.5*L1)) L1 = 2.0d0*(frac(1)+delta)
           if ((frac(2)+delta).gt.( 0.5*L2)) L2 = 2.0d0*(frac(2)+delta)
           if ((frac(3)+delta).gt.( 0.5*L3)) L3 = 2.0d0*(frac(3)+delta)
           if ((frac(1)-delta).lt.(-0.5*L1)) L1 = 2.0d0*(frac(1)+delta)
           if ((frac(2)-delta).lt.(-0.5*L2)) L2 = 2.0d0*(frac(2)+delta)
           if ((frac(3)-delta).lt.(-0.5*L3)) L3 = 2.0d0*(frac(3)+delta)
        end if
      end do

*     **** define unit cell ****
      unita(1,1) = L1*mtensor(1,1)
      unita(2,1) = L1*mtensor(2,1)
      unita(3,1) = L1*mtensor(3,1)
      unita(1,2) = L2*mtensor(1,2)
      unita(2,2) = L2*mtensor(2,2)
      unita(3,2) = L2*mtensor(3,2)
      unita(1,3) = L3*mtensor(1,3)
      unita(2,3) = L3*mtensor(2,3)
      unita(3,3) = L3*mtensor(3,3)

*     **** unknown box type ****
      else
        call errquit('invalid box_type',0,0)
      end if
      

      return
      end

*     ***************************
*     *                         *
*     *   control_notqmmmq      *
*     *                         *
*     ***************************
      logical function control_notqmmmq(string)
      implicit none
      character*16 string

      logical qmmmq

      qmmmq = .false.
      if (index(string,'^').gt.0)   qmmmq = .true.
      if (index(string,'x').eq.1)   qmmmq = .true.
      if (index(string,'X').eq.1)   qmmmq = .true.
      if (index(string,'bq').eq.1)  qmmmq = .true.
      if (index(string,'Bq').eq.1)  qmmmq = .true.
      if (index(string,'bQ').eq.1)  qmmmq = .true.
      if (index(string,'BQ').eq.1)  qmmmq = .true.


      control_notqmmmq = qmmmq
      return
      end



*     ***********************************
*     *					*
*     *		control_move 		*
*     *					*
*     ***********************************
      logical function control_move()
      implicit none 

#include "control.fh"

      control_move = move
      return
      end


*     ***********************************
*     *                                 *
*     *         control_rotation        *
*     *                                 *
*     ***********************************
      logical function control_rotation()
      implicit none

#include "control.fh"

      control_rotation = rotation
      return
      end



*     ***********************************
*     *                                 *
*     *       control_translation       *
*     *                                 *
*     ***********************************
      logical function control_translation()
      implicit none

#include "control.fh"

      control_translation = translation
      return
      end




*     ***********************************
*     *					*
*     *	     control_out_of_time 	*
*     *					*
*     ***********************************

*  This function is used to estimate if there is
* enough time to perform another iteration.  The
* routine control_read intializes this routine.  To
* determine if there is enough time left to do another
* iteration this routine uses estimates for the amount
* of time to finish the simulation (est_finish_time) and 
* the amount of time to perform another iteration step 
* (est_step_time).  Where
*
* est_finish_time = 2*(time elapsed from call to control_read
*                      to the first call to control_out_of_time)
*
* est_step_time = (time elapsed between successive calls to 
*                  control_out_of_time)
*
*  Uses: control_blktime common block located in control.fh
*        util_test_time_remaining
*  created: 5-8-2002

      logical function control_out_of_time()
      implicit none 

#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      integer required_time

*     **** external functions ****
      logical  util_test_time_remaining
      external util_test_time_remaining

      call current_second(cpu2_time)

*     **** This is the first time this routine has been called ****
      if (est_finish_time.eq.-1) then 
         est_finish_time = 2*int(cpu2_time-cpu1_time) ! crude estimate
         value           = .false.

*     **** This routine has been called two or more times ****
      else 
         est_step_time = int(cpu2_time-cpu1_time)+1 ! no statistical info used
         required_time = est_step_time + est_finish_time
         value = .not.util_test_time_remaining(rtdb,required_time)
      end if
         
      cpu1_time = cpu2_time

      control_out_of_time = value
      return
      end


*     ***********************************
*     *					*
*     *		control_frac_coord	*
*     *					*
*     ***********************************
      logical function control_frac_coord()
      implicit none 

#include "control.fh"

      control_frac_coord = frac_coord
      return
      end





*     ***********************************
*     *					*
*     *		control_code 		*
*     *					*
*     ***********************************
      integer function control_code()
      implicit none 

#include "control.fh"

      control_code = code
      return
      end



*     ***********************************
*     *					*
*     *		control_ngrid		*
*     *					*
*     ***********************************
      integer function control_ngrid(ijk)
      implicit none 
      integer ijk

#include "control.fh"

      control_ngrid = ngrid(ijk)
      return
      end




*     ***********************************
*     *					*
*     *		control_it_in		*
*     *					*
*     ***********************************
      integer function control_it_in()
      implicit none 

#include "control.fh"

      control_it_in = loop(1)
      return
      end


*     ***********************************
*     *					*
*     *		control_it_out		*
*     *					*
*     ***********************************
      integer function control_it_out()
      implicit none 

#include "control.fh"

      control_it_out = loop(2)
      return
      end

*     ***********************************
*     *					*
*     *		control_bo_steps_in     *
*     *					*
*     ***********************************
      integer function control_bo_steps_in()
      implicit none 

#include "control.fh"

      control_bo_steps_in = bo_steps(1)
      return
      end

*     ***********************************
*     *					*
*     *		control_bo_steps_out    *
*     *					*
*     ***********************************
      integer function control_bo_steps_out()
      implicit none 

#include "control.fh"

      control_bo_steps_out = bo_steps(2)
      return
      end

*     ***********************************
*     *					*
*     *	     control_bo_algorithm       *
*     *					*
*     ***********************************
      integer function control_bo_algorithm()
      implicit none 

#include "control.fh"

      control_bo_algorithm = bo_algorithm
      return
      end



*     ***********************************
*     *					*
*     *		control_time_step	*
*     *					*
*     ***********************************
      real*8 function control_time_step()
      implicit none 

#include "control.fh"

      control_time_step = time_step
      return
      end

*     ***********************************
*     *					*
*     *		control_bo_time_step	*
*     *					*
*     ***********************************
      real*8 function control_bo_time_step()
      implicit none 

#include "control.fh"

      control_bo_time_step = bo_time_step

      return
      end

*     ***********************************
*     *					*
*     *	    control_ion_time_step	*
*     *					*
*     ***********************************
*     Used by ion.F,  the ion_time_step is
*     set to time_step if Car-Parrinello
*     and it  is set to bo_time_step if 
*     Born-Oppenheimer.
      real*8 function control_ion_time_step()
      implicit none 

#include "control.fh"

      !*** BO dynamics ****
      if (code.eq.11) then
         control_ion_time_step = bo_time_step

      !*** CP dynamics ****
      else
         control_ion_time_step = time_step
      end if

      return
      end



*     ***********************************
*     *					*
*     *		control_fake_mass	*
*     *					*
*     ***********************************
      real*8 function control_fake_mass()
      implicit none 

#include "control.fh"

      control_fake_mass = fake_mass
      return
      end


*     ***********************************
*     *                                 *
*     *         control_bo_fake_mass    *
*     *                                 *
*     ***********************************
      real*8 function control_bo_fake_mass()
      implicit none

#include "control.fh"

      control_bo_fake_mass = bo_fake_mass
      return
      end

*     ***********************************
*     *					*
*     *		control_ks_alpha	*
*     *					*
*     ***********************************
      real*8 function control_ks_alpha()
      implicit none 

#include "control.fh"

      control_ks_alpha = ks_alpha
      return
      end


*     ***********************************
*     *                                 *
*     *         control_scf_algorithm   *
*     *                                 *
*     ***********************************
      integer function control_scf_algorithm()
      implicit none

#include "control.fh"

      control_scf_algorithm = scf_algorithm
      return
      end

*     ***********************************
*     *                                 *
*     *         control_ks_algorithm   *
*     *                                 *
*     ***********************************
      integer function control_ks_algorithm()
      implicit none

#include "control.fh"

      control_ks_algorithm = ks_algorithm
      return
      end


*     ***********************************
*     *					*
*     *		control_ks_maxit_orb 	*
*     *					*
*     ***********************************
      integer function control_ks_maxit_orb()
      implicit none 

#include "control.fh"

      control_ks_maxit_orb = maxit_orb
      return
      end

*     ***********************************
*     *					*
*     *	     control_ks_maxit_orbs 	*
*     *					*
*     ***********************************
      integer function control_ks_maxit_orbs()
      implicit none 

#include "control.fh"

      control_ks_maxit_orbs = maxit_orbs
      return
      end



*     ***********************************
*     *					*
*     *		control_tole		*
*     *					*
*     ***********************************
      real*8 function control_tole()
      implicit none 

#include "control.fh"

      control_tole = tolerances(1)
      return
      end


*     ***********************************
*     *					*
*     *		control_tolc		*
*     *					*
*     ***********************************
      real*8 function control_tolc()
      implicit none 

#include "control.fh"

      control_tolc = tolerances(2)
      return
      end


*     ***********************************
*     *					*
*     *		control_tolr		*
*     *					*
*     ***********************************
      real*8 function control_tolr()
      implicit none 

#include "control.fh"

      control_tolr = tolerances(3)
      return
      end

*     ***********************************
*     *					*
*     *		control_rte		*
*     *					*
*     ***********************************
      real*8 function control_rte()
      implicit none 

#include "control.fh"

      control_rte = scaling(1)
      return
      end

*     ***********************************
*     *					*
*     *		control_rti		*
*     *					*
*     ***********************************
      real*8 function control_rti()
      implicit none 

#include "control.fh"

      control_rti = scaling(2)
      return
      end




*     ***********************************
*     *					*
*     *		control_unita		*
*     *					*
*     ***********************************
      real*8 function control_unita(i,j)
      implicit none 
      integer i,j

#include "control.fh"

      control_unita = unita(i,j)
      return
      end


*     ***********************************
*     *					*
*     *		control_boundry		*
*     *					*
*     ***********************************
      character*12 function control_boundry()
      implicit none 

#include "control.fh"

      control_boundry = boundry
      return
      end


c*     ***********************************
c*     *					*
c*     *		control_pspnames	*
c*     *					*
c*     ***********************************
c      character*20  function control_pspnames(i)
c      implicit none 
c      integer i
c
c#include "control.fh"
c
c      control_pspnames = pspnames(i)
c      return
c      end
c
c
c*     ***********************************
c*     *							 		*
c*     *		control_pspstressnames		*
c*     *									*
c*     ***********************************
c      character*20  function control_pspstressnames(i)
c      implicit none 
c      integer i
c
c      integer ind
c      character*20 pspname
c      character*20 control_pspnames
c      external     control_pspnames
c
c      pspname = control_pspnames(i)
c      ind = index(pspname,' ') -1
c      pspname = pspname(1:ind)//'2'
c
c      control_pspstressnames = pspname
c      return
c      end

*     ***********************************
*     *					*
*     *		control_npsp		*
*     *					*
*     ***********************************
      integer  function control_npsp()
      implicit none 

#include "control.fh"

      control_npsp = npsp
      return
      end



*     ***********************************
*     *					*
*     *		control_ecut		*
*     *					*
*     ***********************************
      real*8 function control_ecut()
      implicit none 

#include "control.fh"
      control_ecut = ecut
      return
      end



*     ***********************************
*     *					*
*     *		control_wcut		*
*     *					*
*     ***********************************
      real*8 function control_wcut()
      implicit none 

#include "control.fh"

      control_wcut = wcut
      return
      end


*     ***********************************
*     *					*
*     *		control_rcut		*
*     *					*
*     ***********************************
      real*8 function control_rcut()
      implicit none 

#include "control.fh"

      control_rcut = rcut
      return
      end

*     ***********************************
*     *					*
*     *		control_ncut		*
*     *					*
*     ***********************************
      integer function control_ncut()
      implicit none 

#include "control.fh"

      control_ncut = ncut
      return
      end




*     ***********************************
*     *					*
*     *		control_output_psi	*
*     *					*
*     ***********************************
      character*30 function control_output_psi()
      implicit none

#include "control.fh"

      control_output_psi = output_wavefunction_filename
      return 
      end

*     ***********************************
*     *					*
*     *		control_output_epsi	*
*     *					*
*     ***********************************
      character*30 function control_output_epsi()
      implicit none

#include "control.fh"

      control_output_epsi = output_ewavefunction_filename
      return 
      end


*     ***********************************
*     *					*
*     *		control_input_psi	*
*     *					*
*     ***********************************
      character*30 function control_input_psi()
      implicit none

#include "control.fh"

      control_input_psi = input_wavefunction_filename
      return 
      end



*     ***********************************
*     *					*
*     *		control_input_epsi	*
*     *					*
*     ***********************************
      character*30 function control_input_epsi()
      implicit none

#include "control.fh"

      control_input_epsi = input_ewavefunction_filename
      return 
      end


*     ***********************************
*     *					*
*     *		control_output_v_psi	*
*     *					*
*     ***********************************
      character*30 function control_output_v_psi()
      implicit none

#include "control.fh"

      control_output_v_psi = output_v_wavefunction_filename
      return 
      end


*     ***********************************
*     *					*
*     *		control_input_v_psi	*
*     *					*
*     ***********************************
      character*30 function control_input_v_psi()
      implicit none

#include "control.fh"

      control_input_v_psi = input_v_wavefunction_filename
      return 
      end



*     ***********************************
*     *					*
*     *		control_xyz		*
*     *					*
*     ***********************************
      character*30 function control_xyz()
      implicit none


#include "control.fh"

      control_xyz = xyz_filename
      return 
      end

*     ***********************************
*     *                                 *
*     *         control_cell_name       *
*     *                                 *
*     ***********************************
      character*30 function control_cell_name()
      implicit none


#include "control.fh"

      control_cell_name = cell_name
      return
      end



*     ***********************************
*     *					*
*     *		control_gga		*
*     *					*
*     ***********************************
      integer function control_gga()
      implicit none

#include "control.fh"

      control_gga = gga
      return 
      end


*     ***********************************
*     *					*
*     *		control_multiplicity	*
*     *					*
*     ***********************************
      integer function control_multiplicity()
      implicit none

#include "control.fh"

      control_multiplicity = multiplicity
      return 
      end

*     ***********************************
*     *					*
*     *	    control_multiplicity_set	*
*     *					*
*     ***********************************
      subroutine control_multiplicity_set(new_multiplicity)
      implicit none
      integer new_multiplicity

#include "control.fh"

      multiplicity = new_multiplicity
      return 
      end

*     *****************************************
*     *                                       *
*     *   control_check_charge_multiplicity   *
*     *                                       *
*     *****************************************
      logical function control_check_charge_multiplicity()
      implicit none

#include "control.fh"

*    *** local variables ***
      logical check
      real*8  icharge,tcharge,t
      integer mult,x,x_wf,nextra_orbs
      integer ispin_wf,ne_wf(2),x_f

*     ***** local functions ****
      logical  psi_filefind
      external psi_filefind
      real*8   control_TotalCharge
      external control_TotalCharge
      real*8   ion_TotalCharge_qm
      external ion_TotalCharge_qm
      !integer  control_frac_occ_extra_orbitals
      !external control_frac_occ_extra_orbitals
      integer  control_fractional_orbitals
      external control_fractional_orbitals



*     **** check wavefunction file ****
      if (psi_filefind()) then

*        **** get mult and e-charge from wavefunction file ****
         call psi_get_ne(ispin_wf,ne_wf)

         !nextra_orbs = control_frac_occ_extra_orbitals()
         if (ispin_wf.eq.1) then
            ne_wf(1) = ne_wf(1) - control_fractional_orbitals(1)
         else
            ne_wf(1) = ne_wf(1) - control_fractional_orbitals(1)
            ne_wf(2) = ne_wf(2) - control_fractional_orbitals(2)
         end if

         x_wf = ne_wf(1)+ne_wf(2)
         mult = ne_wf(1)-ne_wf(2) + 1
         if (ispin_wf.eq.1) then
            x_wf = 2*x_wf
            mult = 1
         end if

*        **** get mult and e-charge from control ****
         tcharge = control_TotalCharge()
         icharge = ion_TotalCharge_qm()
         t = icharge - tcharge       !** total number of electrons **
         x = int(NINT(t))


*        **** reassign spin to agree with total number of electrons ****
         if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
            ispin = 2
         end if

*        **** reassign multiplicity to agree with total number of electrons ****
*        *** odd number of electrons and mult odd ***
         if ((mod(x,2).ne.0) .and.(mod(multiplicity,2).ne.0)) then 
            multiplicity = multiplicity - 1
            do while (multiplicity.gt.(x+1))
               multiplicity = multiplicity - 2
            end do
            if (multiplicity.lt.1) multiplicity = 2
         end if
*        *** even number of electrons and mult even ***
         if ((mod(x,2).eq.0) .and.(mod(multiplicity,2).eq.0)) then 
            multiplicity = multiplicity - 1
            do while (multiplicity.gt.(x+1))
               multiplicity = multiplicity - 2
            end do
            if (multiplicity.lt.1) multiplicity = 1
         end if


*        **** compare multiplicity, charge, and ispin ****
         check = ((mult.eq.multiplicity).and.
     >            (x_wf.eq.x)           .and.
     >            (ispin.eq.ispin_wf))

         if (ispin_wf.eq.3) then
            check=((ne_wf(1).eq.ne_wf(2)).and.
     >             (x.eq.ne_wf(1) ) )
         end if
*     **** no wavefunction file ***
      else
         check = .false.
      end if
       
      control_check_charge_multiplicity = check
      return
      end



*     *****************************************
*     *                                       *
*     *   control_check_number_virtuals       *
*     *                                       *
*     *****************************************
      logical function control_check_number_virtuals()
      implicit none

#include "control.fh"

*    *** local variables ***
      logical check
      integer ispin_wf,ne_wf(2),ne(2)

*     ***** local functions ****
      logical  epsi_filefind
      external epsi_filefind
      integer  control_excited_ne
      external control_excited_ne

*     **** check wavefunction file ****
      if (epsi_filefind()) then

*        **** get mult and e-charge from wavefunction file ****
         call psi_get_ne_excited(ispin_wf,ne_wf)
         ne(1) = 0
         ne(2) = 0
         ne(1) = control_excited_ne(1)
         if (ispin.eq.2) ne(2) = control_excited_ne(2)

         check = ((ne(1).eq.ne_wf(1)).and.
     >            (ne(2).eq.ne_wf(2)).and.
     >            (ispin.eq.ispin_wf))

*     **** no wavefunction file ***
      else
         check = .false.
      end if

      control_check_number_virtuals = check
      return
      end



*     ***********************************
*     *					*
*     *		control_ispin  		*
*     *					*
*     ***********************************
      integer function control_ispin()
      implicit none

#include "control.fh"

      control_ispin = ispin
      return 
      end

*     ***********************************
*     *					*
*     *		control_ispin_set	*
*     *					*
*     ***********************************
      subroutine control_ispin_set(new_ispin)
      implicit none
      integer new_ispin

#include "control.fh"

      ispin = new_ispin
      return 
      end


*     *******************************************
*     *						*
*     *		control_gradient_iterations	*
*     *						*
*     *******************************************
      subroutine control_gradient_iterations()
      implicit none

#include "control.fh"

      loop(1) = 1
      loop(2) = 1

      return 
      end

*     ***********************************
*     *					*
*     *		control_version		*
*     *					*
*     ***********************************
      integer function control_version()
      implicit none

#include "inp.fh"
#include "control.fh"

*     **** local variables ****
      integer l,version

      l =index(boundry,' ') - 1

      version = 3
      if (inp_compare(.false.,boundry(1:l),'periodic'))  version=3
      if (inp_compare(.false.,boundry(1:l),'aperiodic')) version=4

      control_version = version
      return 
      end


*     ************************
*     *                	     *
*     *     control_Nose     *
*     *                      *
*     ************************
      logical function control_Nose()
      implicit none

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose


      control_Nose = nose
      return
      end


*     ****************************
*     *                	 	 *
*     *     control_Nose_Pe      *
*     *                 	 *
*     ****************************
      real*8 function control_Nose_Pe()
      implicit none

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose


      control_Nose_Pe = Pe
      return
      end 


*     ****************************
*     *                	 		 *
*     *     control_Nose_Te      *
*     *                 	 	 *
*     ****************************
      real*8 function control_Nose_Te()
      implicit none

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose


      control_Nose_Te = Te
      return
      end 


*     ****************************
*     *                	 	 *
*     *     control_Nose_Pr      *
*     *                 	 *
*     ****************************
      real*8 function control_Nose_Pr()
      implicit none

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose


      control_Nose_Pr = Pr
      return
      end 


*     ****************************
*     *                	 	 *
*     *     control_Nose_Tr      *
*     *                 	 *
*     ****************************
      real*8 function control_Nose_Tr()
      implicit none

*     **** control_nose common block ****
      logical nose
      real*8 Pe,Te,Pr,Tr
      common / control_nblock / Pe,Te,Pr,Tr,nose


      control_Nose_Tr = Tr
      return
      end 


*     ****************************
*     *                	 	 *
*     *     control_Mulliken     *
*     *                 	 *
*     ****************************
      logical function control_Mulliken()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "rtdb.fh"

#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     ***** local variables ****
      logical value

      value = .false.
      if (code.eq.1) then
        if (.not.rtdb_get(rtdb,'cpsd:mulliken',mt_log,1,value))
     >      value = .false.
      end if
      if ((code.eq.2).or.(code.eq.7)) then
        if (.not.rtdb_get(rtdb,'cpmd:mulliken',mt_log,1,value))
     >     value = .false.
      end if


      if ((code.eq.3).or.(code.eq.11)) then
        if (.not.rtdb_get(rtdb,'cgsd:mulliken',mt_log,1,value))
     >     value = .false.
      end if


      control_Mulliken = value
      return
      end




*     *****************************
*     *                	 	  *
*     * control_allow_translation *
*     *                 	  *
*     *****************************
      logical function control_allow_translation()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"

      logical value
      character*30 operation

      logical  control_qmmm
      external control_qmmm

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.rtdb_get(rtdb,'cgsd:allow_translation',
     >                  mt_log,1,value))
     >  value = .true.

*      *** read the current operation ****
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '

*     *** allow translation of operation == freq||hessian ***
      if (inp_compare(.false.,'freq',operation)) value = .true.
      if (inp_compare(.false.,'hessian',operation)) value = .true.

*     *** allow translation for QMMM ***
      if (control_qmmm()) value = .true.


      control_allow_translation = value
      return
      end




*     *****************************
*     *                           *
*     *        control_qmmm       *
*     *                           *
*     *****************************
      logical function control_qmmm()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"

      logical task_qmmm

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if( .not. rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     >  task_qmmm = .false.

      control_qmmm = task_qmmm
      return
      end


*     ****************************
*     *                	 	 *
*     *  control_num_kvectors    *
*     *                 	 *
*     ****************************
      integer function control_num_kvectors()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*30 zone_name
      character*30 rtdb_name
      integer num_kvectors,l

      if (.not.rtdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >   num_kvectors = 1

      control_num_kvectors = num_kvectors
      return
      end

*     ****************************
*     *                	 	 *
*     *      control_ksvector	 *
*     *                 	 *
*     ****************************
      subroutine control_ksvector(i,ks)
      implicit none
      integer i
      real*8 ks(4)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*30 zone_name
      character*30 rtdb_name
      integer num_kvectors,l
      integer kvs(2)

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      num_kvectors = control_num_kvectors()

      if (.not.MA_push_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >  call errquit('control_ksvector: out of stack', 0,MA_ERR)

      if (.not.rtdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'


      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                  (4*num_kvectors),dbl_mb(kvs(1)))) then
         call dcopy(4*num_kvectors,0.0d0,0,dbl_mb(kvs(1)),1)
      end if

      ks(1) = dbl_mb(kvs(1)+4*(i-1))
      ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
      ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
      ks(4) = dbl_mb(kvs(1)+4*(i-1)+3)

      if (.not.MA_pop_stack(kvs(2)))
     >  call errquit('control_ksvector: failed to free stack',0,MA_ERR)

      return
      end

*     ****************************
*     *                	 	 *
*     *      control_kvector	 *
*     *                 	 *
*     ****************************
      subroutine control_kvector(i,kv)
      implicit none
      integer i
      real*8  kv(3)

*     **** local variables ****
      real*8 ks(4)

*     **** external functions ****
      real*8   lattice_unitg
      external lattice_unitg

      call control_ksvector(i,ks)

      kv(1) = ks(1)*lattice_unitg(1,1)
     >      + ks(2)*lattice_unitg(1,2)
     >      + ks(3)*lattice_unitg(1,3)
      kv(2) = ks(1)*lattice_unitg(2,1)
     >      + ks(2)*lattice_unitg(2,2)
     >      + ks(3)*lattice_unitg(2,3)
      kv(3) = ks(1)*lattice_unitg(3,1)
     >      + ks(2)*lattice_unitg(3,2)
     >      + ks(3)*lattice_unitg(3,3)

      return
      end

*     **********************************
*     *                                *
*     *  control_monkhorst_pack_grid   *
*     *                                *
*     **********************************
      subroutine control_monkhorst_pack_grid(grid)
      implicit none
      integer grid(3)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*30 zone_name
      character*30 rtdb_name
      integer l

      if (.not.rtdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':monkhorst-pack'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,3,grid)) then
         grid(1) = 0
         grid(2) = 0
         grid(3) = 0
      end if

      return
      end


*     **********************************
*     *                                *
*     *      control_ksvector_index    *
*     *                                *
*     **********************************
      integer function control_ksvector_index(ks)
      implicit none
      real*8 ks(*)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*30 zone_name
      character*30 rtdb_name
      integer num_kvectors,i,i1,i2
      integer kvs(2)
      real*8  d1,d2

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      num_kvectors = control_num_kvectors()

      if (.not.MA_push_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >  call errquit('control_ksvector: out of stack', 0,MA_ERR)

      if (.not.rtdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      i = index(zone_name,' ') -1
      rtdb_name = zone_name(1:i)//':kvectors'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                  (4*num_kvectors),dbl_mb(kvs(1)))) then
         call dcopy(4*num_kvectors,0.0d0,0,dbl_mb(kvs(1)),1)
      end if

      i1 = -1
      i2 = -1
      do i=1,num_kvectors
         d1 = dsqrt((dbl_mb(kvs(1)+4*(i-1))   - ks(1))**2
     >            + (dbl_mb(kvs(1)+4*(i-1)+1) - ks(2))**2
     >            + (dbl_mb(kvs(1)+4*(i-1)+2) - ks(3))**2)
         if (d1.lt.1.0e-6) i1 = i
         d2 = dsqrt((dbl_mb(kvs(1)+4*(i-1))   + ks(1))**2
     >            + (dbl_mb(kvs(1)+4*(i-1)+1) + ks(2))**2
     >            + (dbl_mb(kvs(1)+4*(i-1)+2) + ks(3))**2)
         if (d2.lt.1.0e-6) i2 = i
      end do
    
      if (.not.MA_pop_stack(kvs(2)))
     >  call errquit('control_ksvector: failed to free stack',0,MA_ERR)

      if (i1.eq.-1) i1 = i2
      control_ksvector_index = i1
      return
      end




*     *****************************
*     *                	 	  *
*     *    control_TotalCharge	  *
*     *                 	  *
*     *****************************
      real*8 function control_TotalCharge()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      double precision charge

      charge = 0.0d0 
      if (.not.rtdb_get(rtdb,'charge',mt_dbl,1,charge)) then
         charge = 0.0d0 
      end if

      control_TotalCharge = charge
      return
      end



*     **************************************
*     *                                    *
*     *   control_frac_occ_extra_orbitals  *
*     *                                    *
*     **************************************
      integer function control_frac_occ_extra_orbitals()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer norbs

      norbs = 0
      if (.not.
     > rtdb_get(rtdb,'nwpw:frac_occ:extra_orbitals',
     >           mt_int,1,norbs)) then
         norbs = 0
      end if

      control_frac_occ_extra_orbitals = norbs
      return
      end




*     *****************************
*     *                	 	  *
*     *       control_rtdb	  *
*     *                 	  *
*     *****************************
      integer function control_rtdb()
      implicit none

*     **** control_rtdb common block ****
      integer trtdb
      common / control_rtdb1 / trtdb

      control_rtdb = trtdb
      return
      end

*     *****************************
*     *                	 	  *
*     *    	control_minimizer *
*     *                 	  *
*     *****************************
      integer function control_minimizer()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer minimizer

      if (.not.rtdb_get(rtdb,'nwpw:minimizer',mt_int,1,minimizer)) then
         minimizer = 1   ! make the default Grassmann cg
      end if

      control_minimizer = minimizer
      return
      end




*     *****************************
*     *                	          *
*     *    control_lmbfgs_size    *
*     *                           *
*     *****************************
      integer function control_lmbfgs_size()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer lmbfgs_size

      if (.not.rtdb_get(rtdb,'nwpw:lmbfgs_size',mt_int,1,lmbfgs_size)) 
     >  then
         lmbfgs_size = 1
      end if

      control_lmbfgs_size = lmbfgs_size
      return
      end

*     *****************************
*     *                	 		  *
*     *    	control_precondition  *
*     *                 	 	  *
*     *****************************
      logical function control_precondition()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical precondition

      if (.not.rtdb_get(rtdb,'nwpw:precondition',
     >                  mt_log,1,precondition)) 
     >  then
         precondition = .false.
      end if

      control_precondition = precondition
      return
      end

*     *****************************
*     *                	 		  *
*     *    	control_lmbfgs_ondisk *
*     *                 	 	  *
*     *****************************
      logical function control_lmbfgs_ondisk()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical ondisk

      if (.not.rtdb_get(rtdb,'nwpw:lmbfgs_ondisk',
     >                  mt_log,1,ondisk)) 
     >  then
         ondisk = .false.
      end if

      control_lmbfgs_ondisk = ondisk
      return
      end

*     *****************************
*     *                	          *		 
*     *   control_pspparameters   *
*     *                 	  *
*     *****************************
      subroutine control_pspparameters(atom,lmax,locp,rlocal)
      implicit none
      character*4 atom
      integer     lmax,locp
      real*8      rlocal

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer l,k
      character*5  element
      character*20 rtdb_name

      element = '     '
      element = atom
      l = index(element,' ') - 1
      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':lmax'
      k = index(rtdb_name,' ') - 1
      if (.not.rtdb_get(rtdb,rtdb_name(1:k),mt_int,1,lmax)) 
     >  lmax = -1

      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':locp'
      k = index(rtdb_name,' ') - 1
      if (.not.rtdb_get(rtdb,rtdb_name(1:k),mt_int,1,locp)) 
     >  locp = -1

      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':rlocal'
      k = index(rtdb_name,' ') - 1
      if (.not.rtdb_get(rtdb,rtdb_name(1:k),mt_dbl,1,rlocal)) 
     >  rlocal = 1.0d0
      
      return
      end


*     ***************************
*     *                	   	*
*     *        control_Ep	*
*     *                  	*
*     ***************************
      real*8 function control_Ep()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 Ep

      if (.not.rtdb_get(rtdb,'nwpw:Eprecondition',
     >                  mt_dbl,1,Ep)) 
     >  then
         Ep = 0.5d0
      end if

      control_Ep = Ep
      return
      end


*     ***************************
*     *                	   	*
*     *        control_SA	*
*     *                  	*
*     ***************************
      logical function control_SA()
      implicit none

#include "control.fh"

      control_SA=SA
      return
      end

*     ***************************
*     *                	   	*
*     *     control_SA_decay	*
*     *                  	*
*     ***************************
      real*8 function control_SA_decay(choice)
      implicit none
      integer choice

#include "control.fh"

      control_SA_decay = sa_decay(choice)
      return
      end


*     ***************************
*     *                	   	*
*     *        control_Fei	*
*     *                  	*
*     ***************************
      logical function control_Fei()
      implicit none

#include "control.fh"

      control_Fei=fei
      return
      end


*     ***************************
*     *                	   	*
*     *    control_Fei_quench	*
*     *                  	*
*     ***************************
      logical function control_Fei_quench()
      implicit none

#include "control.fh"

      control_Fei_quench=fei_quench
      return
      end


*     ***************************
*     *                	   	*
*     *    control_gram_schmidt *
*     *                  	*
*     ***************************
      logical function control_gram_schmidt()
      implicit none

#include "control.fh"

      control_gram_schmidt=gram_schmidt
      return
      end


*     ***************************
*     *                         *
*     *     control_balance     *
*     *                         *
*     ***************************
      logical function control_balance()
      implicit none

#include "control.fh"

      control_balance=balance
      return
      end



*     *****************************
*     *                	 	  *
*     *    control_excited_ne     *
*     *                 	  *
*     *****************************
      integer function control_excited_ne(ii)
      implicit none
      integer ii

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer ne(2)

      if (.not.rtdb_get(rtdb,'nwpw:excited_ne',mt_int,2,ne)) then
         ne(1) = 0
         ne(2) = 0
      end if

      control_excited_ne = ne(ii)
      return
      end



*     *************************************
*     *                                   *
*     *    control_fractional_orbitals    *
*     *                                   *
*     *************************************
      integer function control_fractional_orbitals(ii)
      implicit none
      integer ii
#include "control.fh"
      control_fractional_orbitals = frac_ne(ii)
      return
      end
*     *************************************
*     *                                   *
*     *  control_fractional_kT            *
*     *                                   *
*     *************************************
      real*8 function control_fractional_kT()
      implicit none
#include "control.fh"
*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)
      control_fractional_kT = kb*frac_temperature
      return
      end
*     *************************************
*     *                                   *
*     *  control_fractional_temperature   *
*     *                                   *
*     *************************************
      real*8 function control_fractional_temperature()
      implicit none
#include "control.fh"
      control_fractional_temperature = frac_temperature
      return
      end
*     *************************************
*     *                                   *
*     *  control_fractional_smeartype     *
*     *                                   *
*     *************************************
      integer function control_fractional_smeartype()
      implicit none
#include "control.fh"
      control_fractional_smeartype = frac_smeartype
      return
      end
*     *************************************
*     *                                   *
*     *         control_fractional        *
*     *                                   *
*     *************************************
      logical function control_fractional()
      implicit none
#include "control.fh"
      control_fractional = fractional
      return
      end



*     *****************************
*     *                           *
*     *    control_mapping        *
*     *                           *
*     *****************************
      integer function control_mapping()
      implicit none

#include "control.fh"

      control_mapping = mapping
      return
      end


*     *****************************
*     *                           *
*     *    control_np_dimensions  *
*     *                           *
*     *****************************
      integer function control_np_dimensions(i)
      implicit none
      integer i

#include "control.fh"

      control_np_dimensions = np_dimensions(i)
      return
      end

*     *****************************
*     *                           *
*     *    control_np_orbital     *
*     *                           *
*     *****************************
      integer function control_np_orbital()
      implicit none

#include "control.fh"

      control_np_orbital = np_dimensions(2)
      return
      end



*     *****************************
*     *                           *
*     *    control_mapping1d      *
*     *                           *
*     *****************************
      integer function control_mapping1d()
      implicit none

#include "control.fh"

      control_mapping1d = mapping1d
      return
      end



*     *****************************
*     *                	 	  *
*     *    control_oep            *
*     *                 	  *
*     *****************************
      logical function control_oep()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical oep

      if (.not.rtdb_get(rtdb,'nwpw:oep',mt_log,1,oep)) then
         oep = .false.
      end if

      control_oep = oep
      return
      end



*     *****************************
*     *                	 	  *
*     *    control_new_vpsi       *
*     *                 	  *
*     *****************************
      logical function control_new_vpsi()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical new_vpsi

      if (.not.rtdb_get(rtdb,'nwpw:new_vmovecs',mt_log,1,new_vpsi)) then
         new_vpsi = .false.
      end if

      control_new_vpsi = new_vpsi
      return
      end


*     *****************************
*     *                	 	  *
*     *    control_COM_shift      *
*     *                 	  *
*     *****************************
      logical function control_COM_shift()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical com_shift

      if (.not.rtdb_get(rtdb,'nwpw:com_shift',mt_log,1,com_shift)) then
         com_shift = .true.
      end if

      control_COM_shift = com_shift
      return
      end



*     *****************************
*     *                           *
*     *    control_DOS            *
*     *                           *
*     *****************************
      logical function control_DOS()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical dos
      real*8  alpha

      dos = .false.
      if (rtdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) dos = .true.

      control_DOS = dos
      return
      end


*     *****************************
*     *                           *
*     *    control_psi_tmp        *
*     *                           *
*     *****************************
      logical function control_psi_tmp()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical psitmp

      psitmp = .false.
      if (.not.rtdb_get(rtdb,'nwpw:psi_tmp',mt_log,1,psitmp))
     >   psitmp = .false.

      control_psi_tmp = psitmp
      return
      end

*     *****************************
*     *                           *
*     *   control_mulliken_kawai  *
*     *                           *
*     *****************************
      logical function control_mulliken_kawai()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.rtdb_get(rtdb,'nwpw:mulliken_kawai',mt_log,1,value))
     >   value = .false.

      control_mulliken_kawai = value
      return
      end



*     *****************************
*     *                           *
*     *   control_zero_forces     *
*     *                           *
*     *****************************
      logical function control_zero_forces()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.rtdb_get(rtdb,'nwpw:zero_forces',mt_log,1,value))
     >   value = .false.

      control_zero_forces = value
      return
      end



*     ********************************
*     *                              *
*     *  control_dos_grid_structure  *
*     *                              *
*     ********************************
      subroutine control_dos_grid_structure(grid)
      implicit none
      integer grid(3)

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.rtdb_get(rtdb,'band:dos-grid',mt_int,3,grid)) then
        grid(1) = 1
        grid(2) = 1
        grid(3) = 1
      end if

      return
      end




*     ***********************************
*     *					*
*     *	 control_reset_band_structure	*
*     *					*
*     ***********************************
      subroutine control_reset_band_structure()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not. rtdb_cget(rtdb, 'pspw:input bvectors',
     >                    1,input_wavefunction_filename))
     >  input_wavefunction_filename = 'atomic'

      if (.not. rtdb_cget(rtdb, 'pspw:output bvectors',
     >                    1,output_wavefunction_filename))
     >     output_wavefunction_filename = ' '
      if (output_wavefunction_filename.eq.' ')then
         if (input_wavefunction_filename.eq.'atomic')then
           call util_file_prefix('bmovecs',output_wavefunction_filename)
         else
            output_wavefunction_filename = input_wavefunction_filename
         endif
      endif
      if (input_wavefunction_filename.eq.'atomic')then
         input_wavefunction_filename = output_wavefunction_filename
      end if

     
      return 
      end

*     **************************************
*     *                	 	           *
*     *  control_num_kvectors_structure    *
*     *                 	           *
*     **************************************
      integer function control_num_kvectors_structure()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*30 zone_name
      character*30 rtdb_name
      integer num_kvectors,l

      value = rtdb_cget(rtdb,'band_structure:zone_name',1,zone_name)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors)

      if (.not. value)
     >  call errquit('control_num_kvectors_structure: failed', 
     >               0, RTDB_ERR)

      control_num_kvectors_structure = num_kvectors
      return
      end
  

*     ************************************
*     *                	 	         *
*     *      control_ksvector_structure	 *
*     *                 	         *
*     ************************************
      subroutine control_ksvector_structure(i,ks)
      implicit none
      integer i
      real*8 ks(4)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*30 zone_name
      character*30 rtdb_name
      integer num_kvectors,l
      integer kvs(2)

*     **** external functions ****
      integer  control_num_kvectors_structure
      external control_num_kvectors_structure

      num_kvectors = control_num_kvectors_structure()
      value = MA_push_get(mt_dbl,(4*num_kvectors),
     >        'kvs',kvs(2),kvs(1))
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0,
     &       MA_ERR)

      value = value.and.
     >        rtdb_cget(rtdb,'band_structure:zone_name',1,zone_name)
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0,
     &       RTDB_ERR)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':kvectors'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1)))

      if (.not. value)
     >  call errquit('control_ksvector: failed to get kvs', 0,
     &       RTDB_ERR)

      ks(1) = dbl_mb(kvs(1)+4*(i-1))
      ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
      ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
      ks(4) = dbl_mb(kvs(1)+4*(i-1)+3)

      value = value.and.MA_pop_stack(kvs(2))

      if (.not. value)
     >  call errquit('control_ksvector: failed to free stack', 0,
     &       MA_ERR)
      return
      end


*     ************************************
*     *                	 	         *
*     *    control_kvector_structure	 *
*     *                                  * 
*     ************************************
      subroutine control_kvector_structure(i,kv)
      implicit none
      integer i
      real*8  kv(3)

*     **** local variables ****
      real*8 ks(4)

*     **** external functions ****
      real*8   lattice_unitg
      external lattice_unitg

      call control_ksvector_structure(i,ks)

      kv(1) = ks(1)*lattice_unitg(1,1)
     >      + ks(2)*lattice_unitg(1,2)
     >      + ks(3)*lattice_unitg(1,3)
      kv(2) = ks(1)*lattice_unitg(2,1)
     >      + ks(2)*lattice_unitg(2,2)
     >      + ks(3)*lattice_unitg(2,3)
      kv(3) = ks(1)*lattice_unitg(3,1)
     >      + ks(2)*lattice_unitg(3,2)
     >      + ks(3)*lattice_unitg(3,3)

      return
      end



*     *****************************
*     *                           *
*     *    control_print          *
*     *                           *
*     *****************************
      logical function control_print(level)
      implicit none
      integer level


*     **** control_print common block ****
      integer print_level
      common / control_print1 / print_level


      logical value

      if (level.le.print_level) then
         value = .true.
      else
         value = .false.
      end if

      control_print = value
      return
      end


*     *****************************
*     *                           *
*     *    control_reduce_print   *
*     *                           *
*     *****************************
      subroutine control_reduce_print()
      implicit none

*     **** control_print common block ****
      integer print_level
      common / control_print1 / print_level

      print_level = print_level -10
      return
      end

*     *****************************
*     *                           *
*     *    control_up_print       *
*     *                           *
*     *****************************
      subroutine control_up_print()
      implicit none

*     **** control_print common block ****
      integer print_level
      common / control_print1 / print_level

      print_level = print_level + 10
      return
      end




*     ************************************
*     *                                  *
*     *  control_optimize_cell_strategy  *
*     *                                  *
*     ************************************

      integer function control_optimize_cell_strategy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_strategy

      if (.not.rtdb_get(rtdb,'cell_optimize:optimize_strategy',
     >                 mt_int,1,optimize_strategy))
     >  optimize_strategy = 0

      control_optimize_cell_strategy = optimize_strategy
      return
      end




*     *************************************
*     *                                   *
*     *  control_optimize_lattice_vectors *
*     *                                   *
*     *************************************

      integer function control_optimize_lattice_vectors(u,v)
      implicit none
      integer u,v

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_lattice_vectors(3,3)

      if (.not.rtdb_get(rtdb,'cell_optimize:optimize_lattice_vectors',
     >                 mt_int,9,optimize_lattice_vectors)) then

        optimize_lattice_vectors(1,1) =1
        optimize_lattice_vectors(2,1) =1
        optimize_lattice_vectors(3,1) =1
        optimize_lattice_vectors(1,2) =1
        optimize_lattice_vectors(2,2) =1
        optimize_lattice_vectors(3,2) =1
        optimize_lattice_vectors(1,3) =1
        optimize_lattice_vectors(2,3) =1
        optimize_lattice_vectors(3,3) =1
      end if

      control_optimize_lattice_vectors = optimize_lattice_vectors(u,v)
      return
      end

*     *************************************
*     *                                   *
*     *  control_optimize_lattice         *
*     *                                   *
*     *************************************

      integer function control_optimize_lattice(i)
      implicit none
      integer i

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_lattice(6)

      if (.not.rtdb_get(rtdb,'cell_optimize:optimize_lattice',
     >                 mt_int,6,optimize_lattice)) then

        optimize_lattice(1) =1
        optimize_lattice(2) =1
        optimize_lattice(3) =1
        optimize_lattice(4) =1
        optimize_lattice(5) =1
        optimize_lattice(6) =1
      end if

      control_optimize_lattice = optimize_lattice(i)
      return
      end 



*     ************************************
*     *                                  *
*     *      control_lmax_multipole      *
*     *                                  *
*     ************************************

      integer function control_lmax_multipole()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer lmax

      if (.not.rtdb_get(rtdb,'nwpw:lmax_multipole',mt_int,1,lmax))
     >  lmax = -1

      control_lmax_multipole = lmax
      return
      end



*     ************************************
*     *                                  *
*     *      control_pfft3_qsize          *
*     *                                  *
*     ************************************

      integer function control_pfft3_qsize()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer qmax

      if (.not.rtdb_get(rtdb,'nwpw:pfft3_qsize',mt_int,1,qmax))
     >  qmax = 4

      control_pfft3_qsize = qmax
      return
      end


*     ************************************
*     *                                  *
*     *      control_nprj_mult           *
*     *                                  *
*     ************************************

      integer function control_nprj_mult()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer qmax

      if (.not.rtdb_get(rtdb,'nwpw:nprj_mult',mt_int,1,qmax))
     >  qmax = 1

      control_nprj_mult = qmax
      return
      end



*     ************************************
*     *                                  *
*     *         control_symmetry         *
*     *                                  *
*     ************************************

      integer function control_symmetry()
      implicit none

#include "control.fh"

      control_symmetry = symm_number
      return
      end

*     ************************************
*     *                                  *
*     *         control_spin_orbit       *
*     *                                  *
*     ************************************

      logical function control_spin_orbit()
      implicit none

#include "control.fh"

      control_spin_orbit = spin_orbit
      return
      end


*     *****************************
*     *                           *
*     *      control_pressure     *
*     *                           *
*     *****************************
      logical function control_pressure()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.rtdb_get(rtdb,'cpmd:pressure',mt_log,1,value))
     >   value = .false.

      control_pressure = value
      return
      end


*     ***********************************
*     *                                 *
*     *      control_init_velocities    *
*     *                                 *
*     ***********************************
      logical function control_init_velocities()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value
      real*8  temp

      value = .true.
      if (.not.rtdb_get(rtdb,'cpmd:init_velocities_temperature',
     >     mt_dbl,1,temp))
     >   value = .false.

      control_init_velocities = value
      return
      end


*     ***********************************
*     *                                 *
*     *      control_kbpp_ray           *
*     *                                 *
*     ***********************************
      logical function control_kbpp_ray()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.rtdb_get(rtdb,'nwpw:kbpp_ray',
     >     mt_log,1,value))
     >   value = .false.

      control_kbpp_ray = value
      return
      end



*     ***********************************
*     *                                 *
*     *      control_kbpp_filter        *
*     *                                 *
*     ***********************************
      logical function control_kbpp_filter()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.rtdb_get(rtdb,'nwpw:kbpp_filter',
     >     mt_log,1,value))
     >   value = .false.

      control_kbpp_filter = value
      return
      end


*     ***********************************
*     *                                 *
*     *    control_brillioun_ondisk     *
*     *                                 *
*     ***********************************
      logical function control_brillioun_ondisk()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.rtdb_get(rtdb,'nwpw:brillioun_ondisk',
     >     mt_log,1,value))
     >   value = .false.

      control_brillioun_ondisk = value
      return
      end


*     ***********************************
*     *                                 *
*     *    control_mparallelized        *
*     *                                 *
*     ***********************************
      logical function control_mparallelized()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.rtdb_get(rtdb,'nwpw:mparallelized',mt_log,1,value))
     >   value = .false.

      control_mparallelized = value
      return
      end

*     ***********************************
*     *                                 *
*     *        control_bo_cpmd          *
*     *                                 *
*     ***********************************
      logical function control_bo_cpmd()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.rtdb_get(rtdb,'nwpw:bo_cpmd',mt_log,1,value))
     >   value = .false.

      control_bo_cpmd = value
      return
      end




*     **********************************************
*     *                                            *
*     *      control_init_velocities_temperature   *
*     *                                            *
*     **********************************************
      real*8 function control_init_velocities_temperature()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.rtdb_get(rtdb,'cpmd:init_velocities_temperature',
     >                  mt_dbl,1,value))
     >   value = 300.0d0

      control_init_velocities_temperature = value
      return
      end



*     ***********************************
*     *                                 *
*     *      control_wannier_timestep   *
*     *                                 *
*     ***********************************
      real*8 function control_wannier_timestep()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.rtdb_get(rtdb,'wannier:time_step',mt_dbl,1,value))
     >   value = 2.7e-2

      control_wannier_timestep = value
      return
      end



*     ***********************************
*     *                                 *
*     *      control_wannier_maxiter    *
*     *                                 *
*     ***********************************
      integer function control_wannier_maxiter()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer value

      if (.not.rtdb_get(rtdb,'wannier:maxiter',mt_int,1,value))
     >   value = 500

      control_wannier_maxiter = value
      return
      end


*     **********************************************
*     *                                            *
*     *      control_attenuation                   *
*     *                                            *
*     **********************************************
      real*8 function control_attenuation()
      implicit none
#include "control.fh"
      control_attenuation = attenuation
      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_two_component_pseudopotential()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this called to signal that a two_component ppot is in use
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "control.fh"
      two_comp_ppot=.true.
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
      logical function two_component_pseudopotential()
#include "control.fh"
      two_component_pseudopotential=two_comp_ppot
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       



*     ***********************************
*     *                                 *
*     *    control_ecut_wcut_default    *
*     *                                 *
*     ***********************************
      subroutine control_ecut_wcut_default(rtdb,ecut,wcut)
      implicit none
      integer rtdb
      real*8 ecut,wcut

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value
      integer geom,ii,nion,l,h
      real*8 wcut_max,wcut_atom,q,rion(3)
      character*4 aname
      character*16 t,rtdbname
      character*255 sdir_name


*     **** external functions ****
      real*8      nwpw_libgcutoff
      external    nwpw_libgcutoff
      character*4 ion_aname_geom
      external    ion_aname_geom
      logical     parsepointcharge
      external    parsepointcharge

      !*** parse psp information on rtdb ***
      if (.not.rtdb_get(rtdb,'nwpw:psp:cutoff',mt_dbl,1,wcut_max)) then 
         call Parallel_taskid(taskid)
         wcut_max = -1.0d0
         value = geom_create(geom,'geometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'geometry')
         value = value.and.geom_ncent(geom,nion)
         if (.not. value) 
     >   call errquit('control_ecut_wcut_default:cannot load geometry',
     >              0,GEOM_ERR)

         do ii=1,nion
           if (.not.geom_cent_get(geom,ii,t,rion,q))
     >     call errquit('control_ecut_wcut_default:error reading ions',
     >                  0,GEOM_ERR)
           if (.not.parsepointcharge(t)) then
              aname = ion_aname_geom(geom,ii)
              l = index(aname,' ') - 1
              if (l.le.0) l = 4
              rtdbname = aname(1:l)//':cutoff'
              l = l+7
              if (.not.rtdb_get(rtdb,rtdbname(1:l),mt_dbl,1,wcut_atom)) 
     >         then

                 !*** define wcut_atom from psplibrary ***
                 value = rtdb_parallel(.false.)
                 if (taskid.eq.MASTER) then
                    call util_directory_name(sdir_name,.true.,0)
                    h = index(sdir_name,' ') - 1
                    open(unit=99,file=sdir_name(1:h)//'/junk.inp',
     >                   status='unknown')
                   close(unit=99,status='delete')

                   call nwpw_libgeninp(1,aname,
     >                  sdir_name(1:h)//'/junk.inp')
                   wcut_atom=nwpw_libgcutoff(aname)

                   if(.not.rtdb_put(rtdb,rtdbname(1:l),
     >                              mt_dbl,1,wcut_atom))
     >             call errquit(
     >           'control_ecut_wcut_default:cannot write wcut_atom',0,0)

                 end if
                 value = rtdb_parallel(.true.)
                 call Parallel_Brdcst_value(MASTER,wcut_atom)

              end if

              !*** reset wcut_max ***
              if (wcut_atom.gt.wcut_max) then
                 wcut_max = wcut_atom
                 value = rtdb_parallel(.false.)
                 if (taskid.eq.MASTER) then
                 if (.not.rtdb_put(rtdb,'nwpw:psp:cutoff',
     >                             mt_dbl,1,wcut_max)) 
     >           call errquit(
     >          'control_ecut_wcut_default:cannot write wcut_max',0,0)
                 end if
                 value = rtdb_parallel(.true.)
              end if
           end if
         end do
         if (.not. geom_destroy(geom))
     >   call errquit('control_ecut_wcut_default:cannot destroy geom',
     >                0,GEOM_ERR)
      end if

      wcut = wcut_max
      ecut = 2.0d0*wcut

      return
      end

