*
* $Id$
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
#include "bafdecls.fh"
#include "btdb.fh"
#include "control.fh"
#include "nwpwxc.fh"
#include "util.fh"

      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,np_default,value2,value3
      integer ispin0,ne(2),nbrill,np1,np2
      integer nmult(5)
      real*8  error,thresh

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
      character*50 rtdb_unita,rtdb_unitaf,rtdb_ngrid,rtdb_boundry
      character*50 exchange_correlation,rtdb_ngrid_small
      integer i,l,np
!      integer ind ! unused

      integer  control_num_kvectors
      external control_num_kvectors
      logical  control_allow_translation
      external control_allow_translation
      integer  Parallel_threadid
      external Parallel_threadid


      value = .true.

      call Parallel_taskid(taskid)
      call nwpw_timing_start(50)
      value2 = btdb_parallel(.true.)
      code = code_in
      trtdb = rtdb

      call util_print_get_level(print_level)
      !write(*,*) "Print level is ",print_level


      DO 10, i = 1,100
*     **** get parallel mappings ****
      if (.not.btdb_get(rtdb,'nwpw:mapping',mt_int,1,mapping)) then
         mapping = 2
      end if
10    continue

*     **** set mapping1d ****
      if (.not.btdb_get(rtdb,'nwpw:mapping1d',mt_int,1,mapping1d)) then
         mapping1d = 1
      end if

*     **** set np_dimensions ****
      call Parallel_np(np)

**** if NOMPI then no processor groups or communicators ****
#ifdef NOMPI
      np_dimensions(1) = np
      np_dimensions(2) = 1
      np_dimensions(3) = 1
      np_default = .true.
#else
      if(.not.btdb_get(rtdb,'nwpw:np_dimensions',
     >                 mt_int,3,np_dimensions)) then

         np_dimensions(1) = np
         np_dimensions(2) = 1
         np_dimensions(3) = 1
         np_default = .true.
      else
         np_default = .false.
         if (.not.((code.eq.5) .or.
     >             (code.eq.10).or.
     >             (code.eq.13).or.
     >             (code.eq.14))) np_dimensions(3) = 1
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
         if ((code.eq.5) .or.
     >       (code.eq.10).or.
     >       (code.eq.13).or.
     >       (code.eq.14)) np_dimensions(2) = 1

         np_dimensions(1) = np/(np_dimensions(2)*np_dimensions(3))

      endif
#endif

*     **** get balance mapping ****
      if (.not.btdb_get(rtdb,'nwpw:balance',mt_log,1,balance)) then
         balance = .true.
      end if

*     **** get parallel io ****
      if (.not.btdb_get(rtdb,'nwpw:parallel_io',mt_log,1,pio)) then
         pio = .false.
      end if
      !pio = pio.and.(np_dimensions(2).gt.1)

*     **** get fast_erf ****
      if (.not.btdb_get(rtdb,'nwpw:fast_erf',mt_log,1,fast_erf)) then
         fast_erf = .false.
      end if

*     **** get fmm ****
      if (.not.btdb_get(rtdb,'nwpw:fmm',mt_log,1,fmm)) then
         fmm = .false.
      end if
      if (.not.btdb_get(rtdb,'nwpw:fmm_lmax',mt_int,1,fmm_lmax)) then
         fmm_lmax = 10
      end if
      if (.not.btdb_get(rtdb,'nwpw:fmm_lr',mt_int,1,fmm_lr)) then
         fmm_lr = 1
      end if

*     **** get periodic_dipole ****
      if (.not.btdb_get(rtdb,'nwpw:periodic_dipole',
     >                  mt_log,1,periodic_dipole)) then
         periodic_dipole = .false.
      end if

*     **** get smooth_cuoff ****
      smooth_cutoff = .true.
      if (.not.btdb_get(rtdb,'nwpw:smooth_cutoff',
     >                  mt_dbl,2,smooth_cutoff_values)) then
         smooth_cutoff           = .false.
         smooth_cutoff_values(1) = 2.0d0
         smooth_cutoff_values(2) = 4.0d0
      end if

*     **** get hess_model mapping ****
      if (.not.btdb_get(rtdb,'nwpw:hess_model',mt_log,
     >                  1,hess_model)) then
         hess_model = .false.
      end if



*     *********************************
*     **** cpsd and band_sd: stuff ****
*     *********************************
      if ((code.eq.1).or.(code.eq.13)) then
      if (.not.btdb_cget(rtdb,'cpsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

      if (.not.btdb_cget(rtdb,'cpsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)) then
         if (.not.btdb_cget(rtdb,'pspw:input vectors',
     >                      1,input_wavefunction_filename)) then
            call util_file_prefix('movecs',
     >                            input_wavefunction_filename)
         end if
      end if

      if (.not.btdb_cget(rtdb,'cpsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)) then
         if (.not.btdb_cget(rtdb,'pspw:output vectors',
     >                      1,output_wavefunction_filename)) then
            call util_file_prefix('movecs',
     >                            output_wavefunction_filename)
         end if
      end if


      if (.not.btdb_cget(rtdb,'cpsd:exchange_correlation',
     >                  1,exchange_correlation))
     >   exchange_correlation = 'vosko'

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c        call nwpwxc_print()
         has_disp = nwpwxc_has_disp()
      endif
!$OMP end single

#include "control_gga.fh"

      if(.not.btdb_get(rtdb,'cpsd:geometry_optimize',mt_log,1,move))
     >    move = .false.
      if(.not.btdb_get(rtdb,'cpsd:fractional_coordinates',
     >          mt_log,1,frac_coord))
     >     frac_coord = .false.
      if(.not.btdb_get(rtdb,'cpsd:npsp',mt_int,1,npsp))
     >    npsp=0
      if(.not.btdb_get(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass))
     >     fake_mass = 400000.0d0
      if(.not.btdb_get(rtdb,'cpsd:time_step',mt_dbl,1,time_step))
     >     time_step = 5.8d0
      if(.not.btdb_get(rtdb,'cpsd:loop',mt_int,2,loop)) then
       loop(1)=10 
       loop(2)=100
      end if
      if(.not.btdb_get(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances))
     >then
         tolerances(1) = 1.0d-9
         tolerances(2) = 1.0d-9
         tolerances(3) = 1.0d-4
      end if 
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      if(.not.btdb_get(rtdb,'cpsd:ecut',mt_dbl,1,ecut))
     >    ecut=9000.0d0
      if(.not.btdb_get(rtdb,'cpsd:wcut',mt_dbl,1,wcut))
     >  wcut=ecut
      if(.not.btdb_get(rtdb,'cpsd:rcut',mt_dbl,1,rcut))
     >    rcut = 0.0d0
      if(.not.btdb_get(rtdb,'cpsd:ncut',mt_int,1,ncut))
     >    ncut = 1
      if(.not.btdb_get(rtdb,'cpsd:mult',mt_int,1,multiplicity))
     >    multiplicity = 1
      if(.not.btdb_get(rtdb,'cpsd:ispin',mt_int,1,ispin))
     >    ispin=1

      value3 = .false.
      if (.not.btdb_get(rtdb,'nwpw:dof_rotation',mt_log,1,dof_rotation))
     >    then
        dof_rotation = .false.
        value3= .true.
      end if

      if (.not.btdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation)) then
         rotation = .true.
      else
         if (value3) dof_rotation = rotation
      end if


      if (.not.btdb_get(rtdb,'nwpw:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2

      if (.not.btdb_get(rtdb,'cpsd:fei',mt_log,1,fei))
     >  fei = .false.
      if(.not.btdb_cget(rtdb,'cpsd:xyz_filename',1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)


*     ****************************************
*     **** cpmd and band_cpmd code: stuff ****
*     ****************************************
      else if ((code.eq.2) .or. (code.eq.14)) then
      if (.not.btdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

      if(.not.btdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)) then
         if(.not.btdb_cget(rtdb,'pspw:input vectors',
     >                     1,input_wavefunction_filename)) then
            call util_file_prefix('movecs',
     >                            input_wavefunction_filename)
         end if
      end if

      if(.not.btdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)) then
         if(.not.btdb_cget(rtdb,'pspw:output vectors',
     >                     1,output_wavefunction_filename)) then
            call util_file_prefix('movecs',
     >                            output_wavefunction_filename)
         end if
      end if

      if(.not.btdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename)) then
         if(.not.btdb_cget(rtdb,'pspw:input vvectors',
     >                     1,input_v_wavefunction_filename)) then
            call util_file_prefix('vmovecs',
     >                            input_v_wavefunction_filename)
         end if
      end if

      if(.not.btdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename)) then
         if(.not.btdb_cget(rtdb,'pspw:output vvectors',
     >                     1,output_v_wavefunction_filename)) then
            call util_file_prefix('vmovecs',
     >                            output_v_wavefunction_filename)
         end if
      end if

      if(.not.btdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)
!      ind = index(exchange_correlation,' ') - 1
      if(.not.btdb_cget(rtdb,'cpmd:exchange_correlation',
     >                  1,exchange_correlation))
     >  exchange_correlation = 'vosko'

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c         call nwpwxc_print()
          has_disp = nwpwxc_has_disp()
      endif
!$OMP end single

#include "control_gga.fh"


      if(.not. btdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move))
     >    move = .true.
      if(.not. btdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >   frac_coord = .false.
      if(.not. btdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp))
     >    npsp = 0
      if(.not. btdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass))
     >    fake_mass = 800.0d0
      if(.not. btdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step))
     >    time_step = 5.0d0
      if(.not. btdb_get(rtdb,'cpmd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not. btdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)) then
        scaling(1) = 1.0d0
        scaling(2) = 1.0d0
      end if

      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      if(.not. btdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut))
     >    ecut=9000.0d0
      if(.not. btdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut))
     >    wcut = ecut
      if(.not. btdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut))
     >    rcut = 0.0d0
      if(.not. btdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut))
     >    ncut = 1
      SA = .true.
      if (.not.btdb_get(rtdb,'cpmd:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.btdb_get(rtdb,'nwpw:dipole_motion',mt_log,
     >                  1,dipole_motion))
     >  dipole_motion = .false.

      if (.not.btdb_get(rtdb,'cpmd:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.btdb_get(rtdb,'cpmd:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.

      value3 = .false.
      if (.not.btdb_get(rtdb,'nwpw:dof_rotation',mt_log,1,dof_rotation))
     >    then
        dof_rotation = .false.
        value3= .true.
      end if

      if (.not.btdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation)) then
         rotation = .true.
      else
         if (value3) dof_rotation = rotation
      end if


*     **** get thermostat information ****
      if (.not.btdb_get(rtdb,'cpmd:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.btdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe))
     >   Pe = 1200.0d0
      if (.not.btdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.btdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr))
     >   Pr = 1200.0d0
      if (.not.btdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0

      if (.not.btdb_get(rtdb,'nwpw:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2


*     ************************************
*     **** cgsd: stuff or  paw: stuff ****
*     ************************************
      else if ((code.eq.3).or.(code.eq.8).or.
     >         (code.eq.11).or.(code.eq.12).or.(code.eq.15)) then
      if (.not.btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
c
c     **** Figure input/output MO vectors ****
c
      if (.not. btdb_cget(rtdb, 'pspw:input vectors', 
     >                    1,input_wavefunction_filename)) then
         input_wavefunction_filename = 'atomic'
      end if

      if (.not. btdb_cget(rtdb, 'pspw:output vectors', 
     >                    1,output_wavefunction_filename)) then
         output_wavefunction_filename = ' '

         if (output_wavefunction_filename.eq.' ')then
            if (input_wavefunction_filename.eq.'atomic')then
              call util_file_prefix('movecs',
     >                              output_wavefunction_filename)
            else
               output_wavefunction_filename=input_wavefunction_filename
            endif
         endif
         if (input_wavefunction_filename.eq.'atomic')then
            input_wavefunction_filename = output_wavefunction_filename
         end if
      end if

      !**** determing input_ewavefunction_filename ****
      if (.not.btdb_cget(rtdb, 'pspw:input evectors', 
     >                   1,input_ewavefunction_filename)) then
         if (.not.btdb_cget(rtdb,'cgsd:input_ewavefunction_filename',
     >                     1,input_ewavefunction_filename))
     >     call util_file_prefix('emovecs',input_ewavefunction_filename)
      end if


      !**** determing output_ewavefunction_filename ****
      if (.not.btdb_cget(rtdb,'pspw:output evectors', 
     >                   1,output_ewavefunction_filename)) then
         if (.not.btdb_cget(rtdb,'cgsd:output_ewavefunction_filename',
     >                     1,output_ewavefunction_filename)) then
            output_ewavefunction_filename=input_ewavefunction_filename
         end if
      end if


      if (code.eq.11) then
     
         if(.not.btdb_cget(rtdb,'pspw:input vvectors',
     >                     1,input_v_wavefunction_filename)) then
            if(.not.btdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                        1,input_v_wavefunction_filename)) then
               call util_file_prefix('vmovecs',
     >                               input_v_wavefunction_filename)
            end if
         end if

         if(.not.btdb_cget(rtdb,'pspw:output vvectors',
     >                     1,output_v_wavefunction_filename)) then
           if(.not.btdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                        1,output_v_wavefunction_filename)) then

               call util_file_prefix('vmovecs',
     >                               output_v_wavefunction_filename)
           end if
         end if
      end if


      if (.not.btdb_cget(rtdb,'cgsd:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c        call nwpwxc_print()
         has_disp = nwpwxc_has_disp()
      endif
!$OMP end single

#include "control_gga.fh"

*     ***** motion filenames ****
      if(.not.btdb_cget(rtdb,'nwpw:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)


*     **** set Kohn-Sham scf parameters ***
      if (.not. btdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not. btdb_get(rtdb,'nwpw:fractional_alpha',
     >                   mt_dbl,1,fractional_alpha))
     >   fractional_alpha = 1.00d0
      if (.not.btdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 3

      if (.not.btdb_get(rtdb,'nwpw:diis_histories',
     >                  mt_int,1,diis_histories))
     >   diis_histories = 15

      if (.not.btdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   ks_algorithm = 0
      if (.not.btdb_get(rtdb,'nwpw:kerker_g0',
     >                  mt_dbl,1,kerker_g0))
     >   kerker_g0 = 0.0d0
      if (.not.btdb_get(rtdb,'nwpw:precondition',
     >                  mt_log,1,precondition))
     >   precondition = .false.

*     **** set maxit_orb maxit_orbs ***
      if (.not.btdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.btdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 1


      if (.not.btdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp))
     >   npsp = 0
      if (.not.btdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 400000.0d0
      if (.not.btdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step))
     >  time_step = 5.8d0
      if (.not.btdb_get(rtdb,'cgsd:loop',mt_int,2,loop)) then
        loop(1) = 10
        loop(2) = 100
      end if
      if (.not.btdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances))then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if

      if(.not. btdb_get(rtdb,'nwpw:scaling',mt_dbl,2,scaling)) then
        scaling(1) = 1.0d0
        scaling(2) = 1.0d0
      end if

      if (.not.btdb_get(rtdb,'cgsd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >    frac_coord = .false.
  
      if (.not.btdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut))
     >   ecut=9000.0d0
      if (.not.btdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut))
     >   wcut = ecut
      if (.not.btdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut))
     >   rcut = 0.0d0
      if (.not.btdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut))
     >   ncut = 1
      if (.not.btdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity))
     >   multiplicity = 1
      if (.not.btdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin))
     >   ispin = 1


*     **** BO parameterss ***
      if (.not.btdb_get(rtdb,'nwpw:bo_steps',mt_int,2,bo_steps)) then
         bo_steps(1) = 10
         bo_steps(2) = 100
      end if
      if (.not.btdb_get(rtdb,'nwpw:bo_time_step',mt_dbl,1,bo_time_step))
     >   bo_time_step = time_step
      if (.not.btdb_get(rtdb,'nwpw:bo_algorithm',mt_int,1,bo_algorithm))
     >   bo_algorithm = 0
      if (.not.btdb_get(rtdb,'nwpw:bo_fake_mass',mt_dbl,1,bo_fake_mass))
     >   bo_fake_mass = 500.0d0

      value3 = .false.
      if (.not.btdb_get(rtdb,'nwpw:dof_rotation',mt_log,1,dof_rotation))
     >    then
        dof_rotation = .false.
        value3= .true.
      end if

      if (.not.btdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation))  then
         rotation = .true.
      else
         if (value3) dof_rotation = rotation
      end if


      SA = .true.
      if (.not.btdb_get(rtdb,'nwpw:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.btdb_get(rtdb,'nwpw:dipole_motion',mt_log,
     >                  1,dipole_motion))
     >  dipole_motion = .false.

      if (.not.btdb_get(rtdb,'nwpw:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.btdb_get(rtdb,'nwpw:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.



*     **** get thermostat information ****
      if (.not.btdb_get(rtdb,'nwpw:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.btdb_get(rtdb,'nwpw:Pe',mt_dbl,1,Pe))
     >   Pe = 1200.0d0
      if (.not.btdb_get(rtdb,'nwpw:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.btdb_get(rtdb,'nwpw:Pr',mt_dbl,1,Pr))
     >   Pr = 1200.0d0
      if (.not.btdb_get(rtdb,'nwpw:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0


*     ***************************
*     **** pspw_dplot: stuff ****
*     ***************************
      else if (code.eq.4) then
         if (.not.btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
           cell_name = 'cell_default'
         end if
         value = .true.
         if (.not.btdb_cget(rtdb,'pspw_dplot:wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >     call util_file_prefix('movecs',input_wavefunction_filename)
         call psi_get_header(i,ngrid,unita,ispin0,ne)
         call dcopy(9,unita,1,unita_frozen,1)
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
         if(.not.btdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut)) ecut=9000.0d0
         if(.not.btdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut)) wcut=ecut
         rcut = 0.0d0
         ncut = 1
         npsp = 0

c         control_read = value
c         return

*     *********************
*     **** band: stuff ****
*     *********************
      else if (code.eq.5) then
     
      if (.not.btdb_cget(rtdb,'band:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

c     **** Figure input/output MO vectors ****
c
      if (.not.btdb_cget(rtdb,'pspw:input vectors', 
     >                    1,input_wavefunction_filename)) then
         input_wavefunction_filename = 'atomic'
      end if

      if (.not.btdb_cget(rtdb, 'pspw:output vectors', 
     >                    1,output_wavefunction_filename)) then
         output_wavefunction_filename = ' '

         if (output_wavefunction_filename.eq.' ')then
            if (input_wavefunction_filename.eq.'atomic') then
              call util_file_prefix('movecs',
     >                              output_wavefunction_filename)
            else
               output_wavefunction_filename=input_wavefunction_filename
            endif
         end if
      end if

      if (input_wavefunction_filename.eq.'atomic')then
         input_wavefunction_filename = output_wavefunction_filename
      end if

      if (.not.btdb_cget(rtdb,'cgsd:input_ewavefunction_filename',
     >                  1,input_ewavefunction_filename))
     >  call util_file_prefix('emovecs',input_ewavefunction_filename)

      if (.not.btdb_cget(rtdb,'cgsd:output_ewavefunction_filename',
     >                  1,output_ewavefunction_filename))
     >  call util_file_prefix('emovecs',output_ewavefunction_filename)


      if (.not.btdb_cget(rtdb,'band:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c        call nwpwxc_print()
         has_disp = nwpwxc_has_disp()
      endif   
!$OMP end single

#include "control_gga.fh"


      if(.not.btdb_get(rtdb,'band:geometry_optimize',mt_log,1,move))
     >    move = .false.
      if (.not.btdb_get(rtdb,'band:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 400000.0d0
      if (.not.btdb_get(rtdb,'band:time_step',mt_dbl,1,time_step))
     >   time_step = 5.8d0
      if (.not.btdb_get(rtdb,'band:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not.btdb_get(rtdb,'band:tolerances',mt_dbl,3,tolerances)) then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      if (.not.btdb_get(rtdb,'band:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.btdb_get(rtdb,'band:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.btdb_get(rtdb,'band:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.btdb_get(rtdb,'band:ncut',mt_int,1,ncut))
     >  ncut = 1
      if (.not.btdb_get(rtdb,'band:mult',mt_int,1,multiplicity))
     >   multiplicity = 1
      if (.not.btdb_get(rtdb,'band:ispin',mt_int,1,ispin))
     >   ispin = 1

      if (.not.btdb_get(rtdb,'band:spin_orbit',mt_log,1,spin_orbit))
     > spin_orbit=.false.
      if (spin_orbit) ispin=2
      


*     **** set Kohn-Sham scf parameters ***
      if (.not. btdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not. btdb_get(rtdb,'nwpw:fractional_alpha',
     >                   mt_dbl,1,fractional_alpha))
     >   fractional_alpha = 1.00d0
      if (.not.btdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 3

      if (.not.btdb_get(rtdb,'nwpw:diis_histories',
     >                  mt_int,1,diis_histories))
     >   diis_histories = 15

      if (.not.btdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   ks_algorithm = 0
      if (.not.btdb_get(rtdb,'nwpw:kerker_g0',
     >                  mt_dbl,1,kerker_g0))
     >   kerker_g0 = 0.0d0
      if (.not.btdb_get(rtdb,'nwpw:precondition',
     >                  mt_log,1,precondition))
     >   precondition = .false.


*     **** set maxit_orb maxit_orbs ***
      if (.not.btdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.btdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 1


*     ***********************
*     **** paw_sd: stuff ****
*     ***********************
      else if (code.eq.6) then

      if (.not.btdb_get(rtdb,'cgsd:geometry_optimize',mt_log,1,move))
     >   move = .false.

      if (.not.btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      if (.not.btdb_cget(rtdb,'cgsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if (.not.btdb_cget(rtdb,'cgsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     > call util_file_prefix('movecs',output_wavefunction_filename)
      if (.not.btdb_cget(rtdb,'cgsd:exchange_correlation',
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

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c        call nwpwxc_print()
         has_disp = nwpwxc_has_disp()
      endif   
!$OMP end single

      if (.not.btdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp))
     >   npsp = 0
      if (.not.btdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass))
     >   fake_mass = 400000.0d0
      if (.not.btdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step))
     >   time_step = 5.8d0
      if (.not.btdb_get(rtdb,'cgsd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if(.not.btdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)) then
         tolerances(1) = 1.0d-7
         tolerances(2) = 1.0d-7
         tolerances(3) = 1.0d-4
      end if
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0
      if (.not.btdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.btdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.btdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.btdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut))
     >  ncut = 1
      if (.not.btdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity))
     >  multiplicity = 1
      if (.not.btdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin))
     >  ispin = 1


*     *************************
*     **** paw_cpmd: stuff ****
*     *************************
      else if (code.eq.7) then

      if (.not.btdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

      if (.not.btdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >  call util_file_prefix('movecs',input_wavefunction_filename)
      if (.not.btdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename))
     >  call util_file_prefix('movecs',output_wavefunction_filename)
      if (.not.btdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',input_v_wavefunction_filename)
      if (.not.btdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename))
     >  call util_file_prefix('vmovecs',output_v_wavefunction_filename)


*     ***** motion filenames ****
      if(.not.btdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename))
     >  call util_file_prefix('xyz',xyz_filename)


      if (.not.btdb_cget(rtdb,'cpmd:exchange_correlation',
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

!$OMP single
      if (nwpwxc_rtdb_load(rtdb,"dft")) then
c        call nwpwxc_print()
         has_disp = nwpwxc_has_disp()
      endif   
!$OMP end single


      if (.not.btdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move))
     >   move = .false.
      if (.not.btdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord))
     >   frac_coord = .false.
      if (.not.btdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp))
     >  npsp = 0
      if (.not.btdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass))
     >  fake_mass = 800.0d0
      if (.not.btdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step))
     >  time_step = 5.0d0
      if (.not.btdb_get(rtdb,'cpmd:loop',mt_int,2,loop)) then
         loop(1) = 10
         loop(2) = 100
      end if
      if (.not.btdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)) then
         scaling(1) = 1.0d0
         scaling(2) = 1.0d0
      end if
      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      if (.not.btdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut))
     >  ecut = 9000.0d0
      if (.not.btdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut))
     >  wcut = ecut
      if (.not.btdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut))
     >  rcut = 0.0d0
      if (.not.btdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut))
     >   ncut = 1

      SA = .true.
      if (.not.btdb_get(rtdb,'cpmd:sa_decay',mt_dbl,2,sa_decay)) then
        SA = .false.
        sa_decay(1) = 1.0d0
        sa_decay(2) = 1.0d0
      end if

      if (.not.btdb_get(rtdb,'nwpw:dipole_motion',mt_log,
     >                  1,dipole_motion))
     >  dipole_motion = .false.

      if (.not.btdb_get(rtdb,'cpmd:fei',mt_log,1,fei))
     >  fei = .false.

      if (.not.btdb_get(rtdb,'cpmd:fei_quench',mt_log,1,fei_quench))
     >  fei_quench = .false.

      value3 = .false.
      if (.not.btdb_get(rtdb,'nwpw:dof_rotation',mt_log,1,dof_rotation))
     >    then
        dof_rotation = .false.
        value3= .true.
      end if

      if (.not.btdb_get(rtdb,'nwpw:rotation',mt_log,1,rotation)) then
         rotation = .true.
      else
         if (value3) dof_rotation = rotation
      end if

*     **** get thermostat information ****
      if (.not.btdb_get(rtdb,'cpmd:nose',mt_log,1,nose))
     >   nose = .false.
      if (.not.btdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe))
     >   Pe = 1200.0d0
      if (.not.btdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te))
     >   Te = 298.15d0
      if (.not.btdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr))
     >   Pr = 1200.0d0
      if (.not.btdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr))
     >   Tr = 298.15d0


*     **********************
*     **** pspw_wannier ****
*     **********************
      else if (code.eq.9) then

      if (.not.btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if

c     **** Figure input/output MO vectors ****

      if (.not.btdb_cget(rtdb,'wannier:input vectors',
     >                  1,input_wavefunction_filename)) then
         if (.not.btdb_cget(rtdb, 'pspw:input vectors',
     >                      1,input_wavefunction_filename)) then
            input_wavefunction_filename = 'atomic'
         end if
      end if

      if (.not.btdb_cget(rtdb,'wannier:output vectors',
     >                  1,output_wavefunction_filename)) then
         if (.not.btdb_cget(rtdb,'pspw:output vectors',
     >                      1,output_wavefunction_filename)) then
            output_wavefunction_filename = ' '
         end if
      end if

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

         call psi_get_header(i,ngrid,unita,ispin0,ne)
         call dcopy(9,unita,1,unita_frozen,1)
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
         if(.not.btdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut)) ecut=9000.0d0
         if(.not.btdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut)) wcut=ecut
         rcut = 0.0d0
         ncut = 0
         npsp = 0

c         control_read = value
c         return

*     ***************************
*     **** band_dplot: stuff ****
*     ***************************
      else if (code.eq.10) then

         if (.not.btdb_cget(rtdb,'band:cell_name',1,cell_name)) then
           cell_name = 'cell_default'
         end if

         value = .true.
         if (.not.btdb_cget(rtdb,'band_dplot:wavefunction_filename',
     >                  1,input_wavefunction_filename))
     >     call util_file_prefix('movecs',input_wavefunction_filename)
         call cpsi_get_header(i,ngrid,unita,ispin0,ne,nbrill)
         call dcopy(9,unita,1,unita_frozen,1)
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
         if(.not.btdb_get(rtdb,'band:ecut',mt_dbl,1,ecut)) ecut=9000.0d0
         if(.not.btdb_get(rtdb,'band:wcut',mt_dbl,1,wcut)) wcut=ecut
         rcut = 0.0d0
         ncut = 0
         npsp = 0

c         control_read = value
c         return

*     ********************************
*     **** unknown but dont crash ****
*     ********************************
      else if (code.gt.99) then
         call nwpw_timing_end(50)
         control_read = value
         return

*     ***************************
*     **** unknown code type ****
*     ***************************
      else
         value = .false.
         write(*,*) "control_read: unknown code type:",code
         call nwpw_timing_end(50)
         control_read = value
         return
      end if

*     *****************************
*     ***** symmetry variables ****
*     *****************************
      if (.not.btdb_get(rtdb,'nwpw:symmetry',mt_int,1,symm_number))
     >   symm_number = 0

*     **********************
*     ***** cell: stuff ****
*     **********************
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      rtdb_unitaf = cell_name(1:l)//':unita_frozen'
      rtdb_ngrid = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'
      rtdb_ngrid_small = cell_name(1:l)//':ngrid_small'

           
*     **** define unita and boundary ****
      if (.not.btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)) then
        call dcopy(9,0.0d0,0,unita,1)
      end if

      if (.not.btdb_cget(rtdb,rtdb_boundry,1,boundry)) then
         boundry = 'periodic'
      end if
      call check_unita_for_default(rtdb,unita,rtdb_unita,cell_name)


*     **** define unita_frozen ****
      if (.not.btdb_get(rtdb,'nwpw:frozen_lattice',mt_log,1,frozen)) 
     >   frozen = .true.

      if (frozen) then
         if (.not.btdb_get(rtdb,rtdb_unitaf,mt_dbl,9,unita_frozen)) then
            call dcopy(9,unita,1,unita_frozen,1)
            value2 = btdb_parallel(.false.)
            if (taskid.eq.MASTER) then
               value = value.and.
     >                 btdb_put(rtdb,rtdb_unitaf,mt_dbl,9,unita_frozen)
            end if
            value2 = btdb_parallel(.true.)
         else
            if (.not.btdb_get(rtdb,'nwpw:frozen_lattice:thresh',
     >                        mt_dbl,1,thresh)) thresh = 0.05d0
            error = 0.0d0
            do l=1,3
            do i=1,3
               error =  error  + ((unita_frozen(i,l)-unita(i,l)))**2
            end do
            end do
            error = dsqrt(error)
            if (error.gt.thresh) then
               call dcopy(9,unita,1,unita_frozen,1)
               value2 = btdb_parallel(.false.)
               if (taskid.eq.MASTER) then
                  value = value.and.
     >                 btdb_put(rtdb,rtdb_unitaf,mt_dbl,9,unita_frozen)
               end if
               value2 = btdb_parallel(.true.)
            end if
         end if
      else
         call dcopy(9,unita,1,unita_frozen,1)
      end if


*     **** define ngrid based on unita and ecut ****
      if (.not.btdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
         if ((ecut.gt.5000.0d0).and.(wcut.gt.5000.0d0)) then
            call control_ecut_wcut_default(rtdb,ecut,wcut)
         else if ((ecut.gt.5000.0d0).and.(wcut.lt.5000.0d0)) then
            ecut = 2.0d0*wcut
         end if
            
         !call control_ngrid_default(rtdb,unita,ecut,mapping,ngrid)
         call control_ngrid_default(rtdb,unita_frozen,ecut,
     >                              mapping,ngrid)
         !ngrid(1) = 32
         !ngrid(2) = 32
         !ngrid(3) = 32
      end if


*     **** define ngrid_small ****
      has_ngrid_small = .false.
      ngrid_small(1) = 0
      ngrid_small(2) = 0
      ngrid_small(3) = 0
      if (btdb_get(rtdb,rtdb_ngrid_small,mt_int,3,ngrid_small))
     >   has_ngrid_small = .true.

*     *** set to false if wannier or dplot code ***
      if ((code.eq.4) .or.(code.eq.9).or.(code.eq.10))
     >   has_ngrid_small = .false.

*     **** set single_p;recision hfx ****
      if (.not.btdb_get(rtdb,'pspw:HFX_single_precision',mt_log,1,
     >                  single_precision_on))
     >   single_precision_on = .false.


*     **** set fractional (smearing) parameters ****
      if (.not.
     > btdb_get(rtdb,'nwpw:fractional_orbitals',mt_int,2,frac_ne)) then
         frac_ne(1) = 0
         frac_ne(2) = 0
      end if
      fractional = (frac_ne(1).gt.0).or.(frac_ne(2).gt.0)
      if (.not.btdb_get(rtdb,'nwpw:fractional_temperature',
     >                  mt_dbl,1,frac_temperature)) then
         frac_temperature = 0.0d0
      end if
      if (.not.btdb_get(rtdb,'nwpw:fractional_smeartype',
     >                  mt_int,1,frac_smeartype)) then
         frac_smeartype = 0
      end if
      if (.not. btdb_get(rtdb,'nwpw:fractional_alpha',
     >                   mt_dbl,1,fractional_alpha)) then
         fractional_alpha = 1.00d0
      end if

*     **** set attenuation parameter ****
      if (.not.btdb_get(rtdb,'nwpw:attenuation',mt_dbl,1,attenuation))
     >  attenuation = 0.5d0

*     **** set preconditioning parameters Ep,Sp ****
      if (.not.btdb_get(rtdb,'nwpw:Eprecondition',mt_dbl,1,Ep))
     >   Ep = 20.0d0
      if (.not.btdb_get(rtdb,'nwpw:Sprecondition',mt_dbl,1,Sp))
     >   Sp = 200.0d0

*     **** set out of time variables ****
      est_step_time   = -1
      est_finish_time = -1
      call current_second(cpu1_time)

*     **** set gram_schmidt ***
      gram_schmidt = .false.
      if (.not.btdb_get(rtdb,
     >      'nwpw:gram-schmidt',mt_log,1,gram_schmidt))
     >  gram_schmidt = .false.

*     **** set translation ***
      translation = control_allow_translation()
      dof_translation = translation
      
*     **** set two component pseudopotential
      two_comp_ppot = .false.      


*     ****  ewald ngrid ************
cccc set to some reasonable default !
      if (.not.btdb_get(rtdb,'nwpw:ewald_ngrid',mt_int,3,ewald_grid)) 
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

      if (np_default) then
         if ( (code.eq.1).or.
     >        (code.eq.2).or.
     >        (code.eq.3).or.
     >        (code.eq.9).or.
     >        (code.eq.11)) then
            nmult(1) = 2
            nmult(2) = 3
            nmult(3) = 5
            do i=1,3
            do while ((np_dimensions(1).gt.ngrid(1)).and.
     >                (mod(np_dimensions(1),nmult(i)).eq.0))
               np_dimensions(1) = np_dimensions(1)/nmult(i)
               np_dimensions(2) = np_dimensions(2)*nmult(i)
            end do
            do while ((np_dimensions(1).gt.ngrid(2)).and.
     >                (mod(np_dimensions(1),nmult(i)).eq.0))
               np_dimensions(1) = np_dimensions(1)/nmult(i)
               np_dimensions(2) = np_dimensions(2)*nmult(i)
            end do
            do while ((np_dimensions(1).gt.ngrid(3)).and.
     >                (mod(np_dimensions(1),nmult(i)).eq.0))
               np_dimensions(1) = np_dimensions(1)/nmult(i)
               np_dimensions(2) = np_dimensions(2)*nmult(i)
            end do
            end do
         end if
      end if
      pio = pio.and.((np_dimensions(2).gt.1).or.(np_dimensions(3).gt.1))

*     **** minimizer ****
      if (.not.btdb_get(rtdb,'nwpw:minimizer',mt_int,1,minimizer)) then
         minimizer = 1   ! make the default Grassmann cg
      end if

      call nwpw_timing_end(50)
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

#include "bafdecls.fh"
#include "btdb.fh"
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

      if (mapping.ge.2) then
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
c      if (.not.btdb_put(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
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
      character*50 rtdb_unita
      character*50 cell_name
      


#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,box_orient
      integer geom,box_type,l,isystype
      real*8  box_delta
      character*50 rtdb_name

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
         value = beom_create(geom,'geometry')
         value = value.and.beom_rtdb_load(rtdb,geom,'geometry')
         value = value.and.geom_amatrix_get(geom,unita)
         value = value.and.geom_systype_get(geom,isystype)
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
           if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,box_delta))
     >       box_delta = 5.0d0
           rtdb_name = cell_name(1:l)//':box_type'
           if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,box_type)) 
     >       box_type = 1
           rtdb_name = cell_name(1:l)//':box_orient'
           if (.not.btdb_get(rtdb,rtdb_name,mt_log,1,
     >                       box_orient)) 
     >       box_orient = .true.
           call control_find_box(geom,box_type,box_orient,box_delta,
     >                           unita)

*          *** write unita to rtdb  - should happen only once during a simulation ***
           call Parallel_taskid(taskid)
           value = btdb_parallel(.false.)
           if (taskid.eq.MASTER) then
           if (.not.btdb_put(rtdb,rtdb_unita,mt_dbl,9,unita)) then
             call errquit('cannot write unita to rtdb',0,0)
           end if
           end if
           value = btdb_parallel(.true.)

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
         value = beom_destroy(geom)
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
#include "errquit.fh"

*     *** local variables ***
      integer lwork
      parameter (lwork=9)
      double precision work(lwork)

      integer ii,nion,ierr,i,j
      double precision q,rxyz(3),mtensor(3,3),frac(3),L1,L2,L3
      double precision x,y,z,meig(3)
      double precision Ixx,Iyx,Izx
      double precision Ixy,Iyy,Izy
      double precision Ixz,Iyz,Izz,Lmax,tmass
      double precision txy,txz,tyz,txx_yy,txx_zz,tyy_zz,tx,ty,tz
     
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
         txy = 0.0d0
         txz = 0.0d0
         tyz = 0.0d0
         txx_yy = 0.0d0
         txx_zz = 0.0d0
         tyy_zz = 0.0d0
         tx = 0.0d0
         ty = 0.0d0
         tz = 0.0d0
         do ii=1,nion
           value = geom_cent_get(geom,ii,t,rxyz,q)
           if (.not.control_notqmmmq(t)) then
            x = rxyz(1)
            y = rxyz(2)
            z = rxyz(3)

            txy = txy + x*y
            txz = txz + x*z
            tyz = tyz + y*z

            txx_yy = txx_yy + (x*x + y*y)
            txx_zz = txx_zz + (x*x + z*z)
            tyy_zz = tyy_zz + (y*y + z*z)

            tx  = tx + x
            ty  = ty + y
            tz  = tz + z
           end if
         end do
         tmass = dble(nion)
         Ixx = tyy_zz - ty*ty/tmass - tz*tz/tmass
         Iyy = txx_zz - tx*tx/tmass - tz*tz/tmass
         Izz = txx_yy - tx*tx/tmass - ty*ty/tmass
         Ixy = txy - tx*ty/tmass 
         Iyx = Ixy
         Ixz = txz - tx*tz/tmass
         Izx = Ixz
         Iyz = tyz - ty*tz/tmass
         Izy = Iyz

         mtensor(1,1) = Ixx
         mtensor(2,1) = -Ixy
         mtensor(3,1) = -Ixz

         mtensor(1,2) = -Iyx
         mtensor(2,2) = Iyy
         mtensor(3,2) = -Iyz

         mtensor(1,3) = -Izx
         mtensor(2,3) = -Izy
         mtensor(3,3) = Izz
 
c        **** longest dimension is along a1 ****
         ierr = 0
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
           if ((frac(1)-delta).lt.(-0.5*L1))
     >        L1 = 2.0d0*(dabs(frac(1))+delta)
           if ((frac(2)-delta).lt.(-0.5*L1)) 
     >        L1 = 2.0d0*(dabs(frac(2))+delta)
           if ((frac(3)-delta).lt.(-0.5*L1)) 
     >        L1 = 2.0d0*(dabs(frac(3))+delta)
        end if
      end do
      !*** put threshold on smallest cubic box ***
      if (L1.lt.26.0d0) L1=26.0d0

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
c           write(*,*) "ii,frac=",ii,frac
c           write(*,*) "                  =",
c     >                mtensor(1,1),mtensor(2,1),mtensor(3,1)
c           write(*,*) "                  =",rxyz
c           write(*,*)
           if ((frac(1)+delta).gt.(0.5d0*L1)) L1 = 2.0d0*(frac(1)+delta)
           if ((frac(2)+delta).gt.(0.5d0*L2)) L2 = 2.0d0*(frac(2)+delta)
           if ((frac(3)+delta).gt.(0.5d0*L3)) L3 = 2.0d0*(frac(3)+delta)
           if ((frac(1)-delta).lt.(-0.5d0*L1)) 
     >        L1 = 2.0d0*(dabs(frac(1))+delta)
           if ((frac(2)-delta).lt.(-0.5d0*L2)) 
     >        L2 = 2.0d0*(dabs(frac(2))+delta)
           if ((frac(3)-delta).lt.(-0.5d0*L3)) 
     >        L3 = 2.0d0*(dabs(frac(3))+delta)
        end if
      end do

c      write(*,*) "L1=",L1
c      write(*,*) "L2=",L2
c      write(*,*) "L3=",L3
      !*** put threshold on smallest box ***
      Lmax = L1
      if (Lmax.lt.L2) Lmax = L2
      if (Lmax.lt.L3) Lmax = L3
      if (Lmax.lt.26.0d0) then
         L1 = 19.0d0
         L2 = 19.0d0
         L3 = 19.0d0
         mtensor(1,1) = 1.0d0
         mtensor(2,1) = 1.0d0
         mtensor(3,1) = 0.0d0
         mtensor(1,2) = 1.0d0
         mtensor(2,2) = 0.0d0
         mtensor(3,2) = 1.0d0
         mtensor(1,3) = 0.0d0
         mtensor(2,3) = 1.0d0
         mtensor(3,3) = 1.0d0
      end if

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
*     *      control_dof_rotation       *
*     *                                 *
*     ***********************************
      logical function control_dof_rotation()
      implicit none

#include "control.fh"

      control_dof_rotation = dof_rotation
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
*     *                                 *
*     *     control_dof_translation     *
*     *                                 *
*     ***********************************
      logical function control_dof_translation()
      implicit none

#include "control.fh"

      control_dof_translation = dof_translation
      return
      end





*     ***********************************
*     *                                 *
*     *    control_translate_vector     *
*     *                                 *
*     ***********************************
      subroutine control_translate_vector(rtrans)
      implicit none
      real*8 rtrans(3)

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'nwpw:translate_vector',
     >                  mt_dbl,3,rtrans)) then
         rtrans(1) = 0.0d0
         rtrans(2) = 0.0d0
         rtrans(3) = 0.0d0
      end if
      return
      end

*     ***********************************
*     *                                 *
*     *    control_translate_reorder    *
*     *                                 *
*     ***********************************
      subroutine control_translate_reorder(reorder)
      implicit none
      logical reorder

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'nwpw:translate_reorder',
     >                  mt_log,1,reorder)) then
         reorder = .true.
      end if
      return
      end



*     ***********************************
*     *                                 *
*     *  control_translate_geom_name    *
*     *                                 *
*     ***********************************
      subroutine control_translate_geom_name(geom_name)
      implicit none
      character*(*) geom_name

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb


      if (.not.btdb_cget(rtdb,'nwpw:translate_geom_name',
     >                  1,geom_name)) then
         geom_name = "translated_geometry"
      end if
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
*     *                                 *
*     *         control_ngrid_small     *
*     *                                 *
*     ***********************************
      integer function control_ngrid_small(ijk)
      implicit none
      integer ijk

#include "control.fh"

      control_ngrid_small = ngrid_small(ijk)
      return
      end

*     ***********************************
*     *                                 *
*     *     control_has_ngrid_small     *
*     *                                 *
*     ***********************************
      logical function control_has_ngrid_small()
      implicit none

#include "control.fh"

      control_has_ngrid_small = has_ngrid_small
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
*     *					*
*     *	   control_fractional_alpha	*
*     *					*
*     ***********************************
      real*8 function control_fractional_alpha()
      implicit none 

#include "control.fh"

      control_fractional_alpha = fractional_alpha
      return
      end


*     ***********************************
*     *					*
*     *		control_kerker_g0	*
*     *					*
*     ***********************************
      real*8 function control_kerker_g0()
      implicit none 

#include "control.fh"

      control_kerker_g0 = kerker_g0
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
*     *       control_diis_histories    *
*     *                                 *
*     ***********************************
      integer function control_diis_histories()
      implicit none

#include "control.fh"

      control_diis_histories = diis_histories
      return
      end

*     ***********************************
*     *                                 *
*     *         control_ks_algorithm    *
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
*     *                                 *
*     *         control_set_unita       *
*     *                                 *
*     ***********************************
      subroutine control_set_unita(unita_in)
      implicit none
      real*8 unita_in(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer l
      character*50 rtdb_unita

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb


      call dcopy(9,unita_in,1,unita,1)

*     **** put unita on the rtdb ****
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'

      call Parallel_taskid(taskid)
      value=btdb_parallel(.false.)
      if (taskid.eq.MASTER) then
         if (.not.btdb_put(rtdb,rtdb_unita,mt_dbl,9,unita))
     >   call errquit('control_set_unita:writing unita',0,RTDB_ERR)
      end if
      value=btdb_parallel(.true.)
      return
      end


*     ***********************************
*     *                                 *
*     *         control_get_unita       *
*     *                                 *
*     ***********************************
      subroutine control_get_unita(unita_out)
      implicit none
      real*8 unita_out(*)

#include "control.fh"

      call dcopy(9,unita,1,unita_out,1)
      return
      end



*     ***********************************
*     *                                 *
*     *    control_set_unita_frozen     *
*     *                                 *
*     ***********************************
      subroutine control_set_unita_frozen(unita_in)
      implicit none
      real*8 unita_in(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer l
      character*50 rtdb_unita

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      call dcopy(9,unita_in,1,unita_frozen,1)

*     **** put unita_frozen on the rtdb ****
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita_frozen'

      call Parallel_taskid(taskid)
      value=btdb_parallel(.false.)
      if (taskid.eq.MASTER) then
         if (.not.btdb_put(rtdb,rtdb_unita,mt_dbl,9,unita))
     >   call errquit('control_set_unita:writing unita',0,RTDB_ERR)
      end if
      value=btdb_parallel(.true.)
      return
      end


*     ***********************************
*     *					*
*     *		control_unita_frozen	*
*     *					*
*     ***********************************
      real*8 function control_unita_frozen(i,j)
      implicit none 
      integer i,j

#include "control.fh"

      control_unita_frozen = unita_frozen(i,j)
      return
      end

*     ***********************************
*     *					*
*     *		control_frozen		*
*     *					*
*     ***********************************
      logical function control_frozen()
      implicit none 

#include "control.fh"

      control_frozen = frozen
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
      character*50 function control_output_psi()
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
      character*50 function control_output_epsi()
      implicit none

#include "control.fh"

      control_output_epsi = output_ewavefunction_filename
      return 
      end

*     ***********************************
*     *					*
*     *		control_output_qpsi	*
*     *					*
*     ***********************************
      character*50 function control_output_qpsi()
      implicit none

#include "control.fh"

      control_output_qpsi = output_qwavefunction_filename
      return 
      end


*     ***********************************
*     *					*
*     *		control_input_psi	*
*     *					*
*     ***********************************
      character*50 function control_input_psi()
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
      character*50 function control_input_epsi()
      implicit none

#include "control.fh"

      control_input_epsi = input_ewavefunction_filename
      return 
      end

*     ***********************************
*     *					*
*     *		control_input_qpsi	*
*     *					*
*     ***********************************
      character*50 function control_input_qpsi()
      implicit none

#include "control.fh"

      control_input_qpsi = input_qwavefunction_filename
      return 
      end


*     ***********************************
*     *					*
*     *		control_output_v_psi	*
*     *					*
*     ***********************************
      character*50 function control_output_v_psi()
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
      character*50 function control_input_v_psi()
      implicit none

#include "control.fh"

      control_input_v_psi = input_v_wavefunction_filename
      return 
      end


*     ***********************************
*     *                                 *
*     *   control_use_fractional_rho    *
*     *                                 *
*     ***********************************
      logical function control_use_fractional_rho()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical use_rho

      if (.not.btdb_get(rtdb,"nwpw:use_fractional_rho", 
     >     mt_log,1,use_rho)) then
         use_rho = .true. 
      end if
      control_use_fractional_rho = use_rho
      return
      end


*     ***********************************
*     *					*
*     *		control_input_rho 	*
*     *					*
*     ***********************************

      character*50 function control_input_rho()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      character*50 filename

      if (.not.btdb_cget(rtdb,'nwpw:input rho',1,filename)) then
         call util_file_prefix('rho',filename)
      end if
      control_input_rho = filename
      return 
      end

*     ***********************************
*     *                                 *
*     *         control_output_rho      *
*     *                                 *
*     ***********************************

      character*50 function control_output_rho()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      character*50 filename

      if (.not.btdb_cget(rtdb,'nwpw:output rho',1,filename)) then
         call util_file_prefix('rho',filename)
      end if
      control_output_rho = filename
      return
      end




*     ***********************************
*     *					*
*     *		control_xyz		*
*     *					*
*     ***********************************
      character*50 function control_xyz()
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
      character*50 function control_cell_name()
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
*     *		control_only_lda	*
*     *					*
*     ***********************************
      logical function control_only_lda()
      implicit none

#include "control.fh"

*     **** external functions ****
      logical  nwpwxc_is_on,nwpwxc_is_lda
      external nwpwxc_is_on,nwpwxc_is_lda

      control_only_lda = (.not.nwpwxc_is_on().and.gga.eq.0).or.
     >                   (nwpwxc_is_on().and.nwpwxc_is_lda())
      return 
      end



*     ***********************************
*     *                                 *
*     *        control_is_grimme2       *
*     *                                 *
*     ***********************************
      logical function control_is_grimme2()
      implicit none

#include "control.fh"

      control_is_grimme2 = is_grimme2
      return
      end


*     ***********************************
*     *                                 *
*     *        control_has_disp         *
*     *                                 *
*     ***********************************
      logical function control_has_disp()
      implicit none

#include "control.fh"

      control_has_disp = has_disp
      return
      end


*     ***********************************
*     *                                 *
*     *        control_options_disp     *
*     *                                 *
*     ***********************************
      character*80 function control_options_disp()
      implicit none

#include "control.fh"

      control_options_disp = options_disp
      return
      end

*     ***********************************
*     *                                 *
*     *        control_has_vdw          *
*     *                                 *
*     ***********************************
      logical function control_has_vdw()
      implicit none

#include "control.fh"

      control_has_vdw = has_vdw
      return
      end 

*     ***********************************
*     *                                 *
*     *        control_is_vdw2          *
*     *                                 *
*     ***********************************
      logical function control_is_vdw2()
      implicit none

#include "control.fh"

      control_is_vdw2 = is_vdw2
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
      logical  control_mult_fixed
      external control_mult_fixed



*     **** check wavefunction file ****
      if (psi_filefind()) then

*        **** get mult and e-charge from wavefunction file ****
         call psi_get_ne(ispin_wf,ne_wf)

c         write(*,*) "ne_wf = ",ne_wf
c         x_wf0 = ne_wf(1)+ne_wf(2)
c         mult0 = ne_wf(1)-ne_wf(2) + 1
         !nextra_orbs = control_frac_occ_extra_orbitals()
         if (ispin_wf.eq.1) then
            ne_wf(1) = ne_wf(1) - control_fractional_orbitals(1)
         else
            ne_wf(1) = ne_wf(1) - control_fractional_orbitals(1)
            ne_wf(2) = ne_wf(2) - control_fractional_orbitals(2)
         end if
c         write(*,*) "after ne_wf = ",ne_wf

         x_wf = ne_wf(1)+ne_wf(2)
         mult = ne_wf(1)-ne_wf(2) + 1
         if (ispin_wf.eq.1) then
c            x_wf0 = 2*x_wf0
c            mult0 = 1
            x_wf  = 2*x_wf
            mult = 1
         end if

*        **** get mult and e-charge from control ****
         tcharge = control_TotalCharge()
         icharge = ion_TotalCharge_qm()
         t = icharge - tcharge       !** total number of electrons **
         x = int(NINT(t))

         if (.not.control_mult_fixed()) then

*           **** reassign spin to agree with total number of electrons ****
            if ((mod(x,2).ne.0).and.(ispin.eq.1)) then !** odd number of electrons **
               ispin = 2
            end if

*           **** reassign multiplicity to agree with total number of electrons ****
*           *** odd number of electrons and mult odd ***
            if ((mod(x,2).ne.0) .and.(mod(multiplicity,2).ne.0)) then 
               multiplicity = multiplicity - 1
               do while (multiplicity.gt.(x+1))
                  multiplicity = multiplicity - 2
               end do
               if (multiplicity.lt.1) multiplicity = 2
            end if
*           *** even number of electrons and mult even ***
            if ((mod(x,2).eq.0) .and.(mod(multiplicity,2).eq.0)) then 
               multiplicity = multiplicity - 1
               do while (multiplicity.gt.(x+1))
                  multiplicity = multiplicity - 2
               end do
               if (multiplicity.lt.1) multiplicity = 1
            end if


*           **** compare multiplicity, charge, and ispin ****
c            check = ( (((mult.eq.multiplicity).and.(x_wf.eq.x))  .or.
c     >                ((mult0.eq.multiplicity).and.(x_wf0.eq.x))).and.
c     >               (ispin.eq.ispin_wf))
            check = ( (mult.eq.multiplicity).and.(x_wf.eq.x).and.
     >               (ispin.eq.ispin_wf))

            if (ispin_wf.eq.3) then
               check=((ne_wf(1).eq.ne_wf(2)).and.
     >                (x.eq.ne_wf(1) ) )
            end if

         else
*           **** compare multiplicity, and ispin ****
            check = ((mult.eq.multiplicity).and.(ispin.eq.ispin_wf))
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


*     *****************************************
*     *                                       *
*     *      control_check_q_orbitals         *
*     *                                       *
*     *****************************************
      logical function control_check_q_orbitals()
      implicit none

#include "control.fh"

*    *** local variables ***
      logical check
      integer ispin_wf,ne_wf(2),ne(2)

*     ***** local functions ****
      logical  qpsi_filefind
      external qpsi_filefind
      integer  psi_ne
      external psi_ne

*     **** check wavefunction file ****
      if (qpsi_filefind()) then

*        **** get mult and e-charge from wavefunction file ****
         call psi_get_ne_qpsi(ispin_wf,ne_wf)
         ne(1) = psi_ne(1)
         ne(2) = psi_ne(2)

         check = ((ne(1).eq.ne_wf(1)).and.
     >            (ne(2).eq.ne_wf(2)).and.
     >            (ispin.eq.ispin_wf))

*     **** no wavefunction file ***
      else
         check = .false.
      end if

      control_check_q_orbitals = check
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     ***** local variables ****
      logical value

      value = .false.
      if (code.eq.1) then
        if (.not.btdb_get(rtdb,'cpsd:mulliken',mt_log,1,value))
     >      value = .false.
      end if
      if ((code.eq.2).or.(code.eq.7)) then
        if (.not.btdb_get(rtdb,'cpmd:mulliken',mt_log,1,value))
     >     value = .false.
      end if


      if ((code.eq.3).or.(code.eq.11)) then
        if (.not.btdb_get(rtdb,'cgsd:mulliken',mt_log,1,value))
     >     value = .false.
      end if

      if ((code.eq.5).or.(code.eq.11)) then
        if (.not.btdb_get(rtdb,'band:mulliken',mt_log,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"

      logical value
      character*50 operation

      logical  control_qmmm
      external control_qmmm

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'cgsd:allow_translation',
     >                  mt_log,1,value))
     >  value = .true.

*      *** read the current operation ****
      if (.not. btdb_cget(rtdb, 'task:operation', 1, operation))
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

      logical task_qmmm

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if( .not. btdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     >  task_qmmm = .false.

      control_qmmm = task_qmmm
      return
      end


*     *****************************
*     *                           *
*     *     control_makehmass2    *
*     *                           *
*     *****************************
      logical function control_makehmass2()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

      logical makehmass2

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if( .not. btdb_get(rtdb,'nwpw:makehmass2',mt_log,1,makehmass2))
     >  makehmass2 = .true.

      control_makehmass2 = makehmass2
      return
      end



*     *****************************
*     *                           *
*     *        control_MP2        *
*     *                           *
*     *****************************
      logical function control_MP2()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

      logical mp2

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if( .not. btdb_get(rtdb,'nwpw:MP2',mt_log,1,mp2))
     >  mp2 = .false.

      control_MP2 = mp2
      return
      end

      
*     *****************************
*     *                           *
*     *     control_2qintegrals   *
*     *                           *
*     *****************************
      logical function control_2qintegrals()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

      logical q2

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if( .not. btdb_get(rtdb,'nwpw:2qintegrals',mt_log,1,q2))
     >  q2 = .false.

      control_2qintegrals = q2
      return
      end


*     ****************************
*     *                	 	 *
*     *  control_num_kvectors    *
*     *                 	 *
*     ****************************
      integer function control_num_kvectors()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors,l

      if (.not.btdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >   num_kvectors = 1

      control_num_kvectors = num_kvectors
      return
      end

*     **********************************
*     *                                *
*     *  control_num_kvectors_print    *
*     *                                *
*     **********************************
      integer function control_num_kvectors_print()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
     
*     **** local variables ****
      logical value
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors_print,l
      
      if (.not.btdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'
      
      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors_print'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors_print))
     >   num_kvectors_print = 0
         
      control_num_kvectors_print = num_kvectors_print
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors,l
      integer kvs(2)

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      num_kvectors = control_num_kvectors()

      if (.not.BA_push_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >  call errquit('control_ksvector: out of stack', 0,MA_ERR)

      if (.not.btdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'


      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,
     >                  (4*num_kvectors),dbl_mb(kvs(1)))) then
         call dcopy(4*num_kvectors,0.0d0,0,dbl_mb(kvs(1)),1)
      end if

      ks(1) = dbl_mb(kvs(1)+4*(i-1))
      ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
      ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
      ks(4) = dbl_mb(kvs(1)+4*(i-1)+3)

      if (.not.BA_pop_stack(kvs(2)))
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*50 zone_name
      character*50 rtdb_name
      integer l

      if (.not.btdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':monkhorst-pack'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,3,grid)) then
         grid(1) = 0
         grid(2) = 0
         grid(3) = 0
      end if

      return
      end

*     **********************************
*     *                                *
*     *      control_fft_dos_grid      *
*     *                                *
*     **********************************
      subroutine control_fft_dos_grid(grid)
      implicit none
      integer grid(3)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*50 zone_name
      character*50 rtdb_name
      integer l

      if (.not.btdb_cget(rtdb,'band_fft:zone_name',
     >                   1,zone_name))
     >   zone_name = 'zone_fft_default'

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':dos-grid'
      if(.not.btdb_get(rtdb,rtdb_name,mt_int,3,grid)) then
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors,i,i1,i2
      integer kvs(2)
      real*8  d1,d2

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      num_kvectors = control_num_kvectors()

      if (.not.BA_push_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >  call errquit('control_ksvector: out of stack', 0,MA_ERR)

      if (.not.btdb_cget(rtdb,'band:zone_name',1,zone_name))
     >   zone_name = 'zone_default'

      i = index(zone_name,' ') -1
      rtdb_name = zone_name(1:i)//':kvectors'
      if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,
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
    
      if (.not.BA_pop_stack(kvs(2)))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      double precision charge

      charge = 0.0d0 
      if (.not.btdb_get(rtdb,'charge',mt_dbl,1,charge)) then
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer norbs

      norbs = 0
      if (.not.
     > btdb_get(rtdb,'nwpw:frac_occ:extra_orbitals',
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
*     *   control_minimizer       *
*     *                 	  *
*     *****************************
      integer function control_minimizer()
      implicit none

#include "control.fh"

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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer lmbfgs_size

      if (.not.btdb_get(rtdb,'nwpw:lmbfgs_size',mt_int,1,lmbfgs_size)) 
     >  then
         lmbfgs_size = 1
      end if

      control_lmbfgs_size = lmbfgs_size
      return
      end

*     *****************************
*     *                		  *
*     *    control_precondition   *
*     *                  	  *
*     *****************************
      logical function control_precondition()
      implicit none

#include "control.fh"

c#include "bafdecls.fh"
c#include "btdb.fh"
c
c*     **** control_rtdb common block ****
c      integer rtdb
c      common / control_rtdb1 / rtdb
c
c      logical precondition
c
c      if (.not.btdb_get(rtdb,'nwpw:precondition',
c     >                  mt_log,1,precondition)) 
c     >  then
c         precondition = .false.
c      end if

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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical ondisk

      if (.not.btdb_get(rtdb,'nwpw:lmbfgs_ondisk',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer l,k
      character*5  element
      character*50 rtdb_name

      element = '     '
      element = atom
      l = index(element,' ') - 1
      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':lmax'
      k = index(rtdb_name,' ') - 1
      if (.not.btdb_get(rtdb,rtdb_name(1:k),mt_int,1,lmax)) 
     >  lmax = -1

      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':locp'
      k = index(rtdb_name,' ') - 1
      if (.not.btdb_get(rtdb,rtdb_name(1:k),mt_int,1,locp)) 
     >  locp = -1

      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':rlocal'
      k = index(rtdb_name,' ') - 1
      if (.not.btdb_get(rtdb,rtdb_name(1:k),mt_dbl,1,rlocal)) 
     >  rlocal = 1.0d0
      
      return
      end

*     **********************************
*     *                                *              
*     *   control_mullikenparameters   *
*     *                                *
*     **********************************
      subroutine control_mullikenparameters(atom,rcut,lmbda)
      implicit none
      character*4 atom
      real*8 rcut,lmbda

#include "bafdecls.fh"
#include "btdb.fh"
      
*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
      
      logical value
      integer l,k
      character*5  element
      character*50 rtdb_name
      
      element = '     '
      element = atom     
      l = index(element,' ') - 1
      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':mulliken:rcut'
      k = index(rtdb_name,' ') - 1
      if (.not.btdb_get(rtdb,rtdb_name(1:k),mt_dbl,1,rcut))
     >  rcut = 1.0d0

      rtdb_name = '                   '
      rtdb_name       = element(1:l)//':mulliken:lmbda'
      k = index(rtdb_name,' ') - 1
      if (.not.btdb_get(rtdb,rtdb_name(1:k),mt_dbl,1,lmbda))
     >  lmbda = 0.0d0
           
      return
      end
         


*     ***************************
*     *                	   	*
*     *        control_Ep	*
*     *                  	*
*     ***************************
      real*8 function control_Ep()
      implicit none

#include "control.fh"

      control_Ep = Ep
      return
      end


      
*     ***************************
*     *                         *
*     *        control_Sp       *
*     *                         *
*     ***************************
      real*8 function control_Sp()
      implicit none

#include "control.fh"

      control_Sp = Sp
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

*     **********************************
*     *                                *
*     *        control_dipole_motion   *
*     *                                *
*     **********************************
      logical function control_dipole_motion()
      implicit none

#include "control.fh"

      control_dipole_motion=dipole_motion
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
*     *                         *
*     *   control_gram_schmidt  *
*     *                         *
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer ne(2)

      if (.not.btdb_get(rtdb,'nwpw:excited_ne',mt_int,2,ne)) then
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

*     *************************************
*     *                                   *
*     *         control_ortho             *
*     *                                   *
*     *************************************
      logical function control_ortho()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical ortho

      if (.not.
     >    btdb_get(rtdb,'nwpw:ortho_initialize',mt_log,1,ortho)) then
         ortho = .true.
      end if

      control_ortho = ortho
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
*     *                           *
*     *    control_parallel_io    *
*     *                           *
*     *****************************
      logical function control_parallel_io()
      implicit none

#include "control.fh"

      control_parallel_io = pio
      return
      end



*     *****************************
*     *                	 	  *
*     *    control_oep            *
*     *                 	  *
*     *****************************
      logical function control_oep()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical oep

      if (.not.btdb_get(rtdb,'nwpw:oep',mt_log,1,oep)) then
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical new_vpsi

      if (.not.btdb_get(rtdb,'nwpw:new_vmovecs',mt_log,1,new_vpsi)) then
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical com_shift

      if (.not.btdb_get(rtdb,'nwpw:com_shift',mt_log,1,com_shift)) then
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical dos
      real*8  alpha

      dos = .false.
      if (btdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) dos = .true.

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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical psitmp

      psitmp = .false.
      if (.not.btdb_get(rtdb,'nwpw:psi_tmp',mt_log,1,psitmp))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:mulliken_kawai',mt_log,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:zero_forces',mt_log,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'band:dos-grid',mt_int,3,grid)) then
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "control.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not. btdb_cget(rtdb, 'pspw:input bvectors',
     >                    1,input_wavefunction_filename))
     >  input_wavefunction_filename = 'atomic'

      if (.not. btdb_cget(rtdb, 'pspw:output bvectors',
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors,l

      value = btdb_cget(rtdb,'band_structure:zone_name',1,zone_name)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors'
      value = value.and.
     >        btdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors)

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

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** local variables ****
      logical value
      character*50 zone_name
      character*50 rtdb_name
      integer num_kvectors,l
      integer kvs(2)

*     **** external functions ****
      integer  control_num_kvectors_structure
      external control_num_kvectors_structure

      num_kvectors = control_num_kvectors_structure()
      value = BA_push_get(mt_dbl,(4*num_kvectors),
     >        'kvs',kvs(2),kvs(1))
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0,
     &       MA_ERR)

      value = value.and.
     >        btdb_cget(rtdb,'band_structure:zone_name',1,zone_name)
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0,
     &       RTDB_ERR)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':kvectors'
      value = value.and.
     >        btdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1)))

      if (.not. value)
     >  call errquit('control_ksvector: failed to get kvs', 0,
     &       RTDB_ERR)

      ks(1) = dbl_mb(kvs(1)+4*(i-1))
      ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
      ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
      ks(4) = dbl_mb(kvs(1)+4*(i-1)+3)

      value = value.and.BA_pop_stack(kvs(2))

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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_strategy

      if (.not.btdb_get(rtdb,'cell_optimize:optimize_strategy',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_lattice_vectors(3,3)

      if (.not.btdb_get(rtdb,'cell_optimize:optimize_lattice_vectors',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer optimize_lattice(6)

      if (.not.btdb_get(rtdb,'cell_optimize:optimize_lattice',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer lmax

      if (.not.btdb_get(rtdb,'nwpw:lmax_multipole',mt_int,1,lmax))
     >  lmax = 0

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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer qmax

      if (.not.btdb_get(rtdb,'nwpw:pfft3_qsize',mt_int,1,qmax))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer qmax

      if (.not.btdb_get(rtdb,'nwpw:nprj_mult',mt_int,1,qmax))
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

*     ************************************
*     *                                  *
*     *         control_fast_erf         *
*     *                                  *
*     ************************************

      logical function control_fast_erf()
      implicit none

#include "control.fh"

      control_fast_erf = fast_erf
      return
      end



*     ************************************
*     *                                  *
*     *         control_fmm              *
*     *                                  *
*     ************************************

      logical function control_fmm()
      implicit none

#include "control.fh"

      control_fmm = fmm
      return
      end

*     ************************************
*     *                                  *
*     *      control_periodic_dipole     *
*     *                                  *
*     ************************************
      logical function control_periodic_dipole()
      implicit none

#include "control.fh"

      control_periodic_dipole = periodic_dipole
      return
      end

*     ************************************
*     *                                  *
*     *         control_smooth_cutoff    *
*     *                                  *
*     ************************************
      logical function control_smooth_cutoff()
      implicit none

#include "control.fh"

      control_smooth_cutoff = smooth_cutoff
      return
      end

*     ************************************
*     *                                  *
*     *    control_smooth_cutoff_values  *
*     *                                  *
*     ************************************
      real*8 function control_smooth_cutoff_values(i)
      implicit none
      integer i

#include "control.fh"

      control_smooth_cutoff_values = smooth_cutoff_values(i)
      return
      end




*     ************************************
*     *                                  *
*     *         control_fmm_lmax         *
*     *                                  *
*     ************************************

      integer function control_fmm_lmax()
      implicit none

#include "control.fh"

      control_fmm_lmax = fmm_lmax
      return
      end

*     ************************************
*     *                                  *
*     *         control_fmm_lr           *
*     *                                  *
*     ************************************

      integer function control_fmm_lr()
      implicit none

#include "control.fh"

      control_fmm_lr = fmm_lr
      return
      end





*     *****************************
*     *                           *
*     *      control_pressure     *
*     *                           *
*     *****************************
      logical function control_pressure()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'cpmd:pressure',mt_log,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical ok,tmp

      if (.not.btdb_get(rtdb,'nwpw:init_velocities',mt_log,1,tmp))
     >   tmp = .false.

      if (tmp) ok = rtdb_delete(rtdb,'nwpw:init_velocities')

      control_init_velocities = tmp
      return
      end


*     ***********************************
*     *                                 *
*     *      control_kbpp_ray           *
*     *                                 *
*     ***********************************
      logical function control_kbpp_ray()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:kbpp_ray',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:kbpp_filter',
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:brillioun_ondisk',
     >     mt_log,1,value))
     >   value = .false.

      control_brillioun_ondisk = value
      return
      end



*     **************************************
*     *                                    *
*     *   control_brillioun_use_symmetry   *
*     *                                    *
*     **************************************
      logical function control_brillioun_use_symmetry()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:brillioun_use_symmetry',
     >     mt_log,1,value))
     >   value = .true.

      control_brillioun_use_symmetry = value
      return
      end

*     **************************************
*     *                                    *
*     *      control_rho_use_symmetry      *
*     *                                    *
*     **************************************
      logical function control_rho_use_symmetry()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"
      
*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value
      
      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:rho_use_symmetry',
     >     mt_log,1,value))
     >   value = .false.

      control_rho_use_symmetry = value
      return                
      end




*     ***********************************
*     *                                 *
*     *    control_mparallelized        *
*     *                                 *
*     ***********************************
      logical function control_mparallelized()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:mparallelized',mt_log,1,value))
     >   value = .false.

      control_mparallelized = value
      return
      end


*     ***********************************
*     *                                 *
*     *    control_mreplicate_size      *
*     *                                 *
*     ***********************************
      integer function control_mreplicate_size()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer value

      value = 16
      if (.not.btdb_get(rtdb,'nwpw:mreplicate_size',mt_int,1,value))
     >   value = 16

      control_mreplicate_size = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_bo_cpmd          *
*     *                                 *
*     ***********************************
      logical function control_bo_cpmd()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:bo_cpmd',mt_log,1,value))
     >   value = .false.

      control_bo_cpmd = value
      return
      end

*     ***********************************
*     *                                 *
*     *        control_hfxon_virtual    *
*     *                                 *
*     ***********************************
      logical function control_hfxon_virtual()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:hfxon_virtual',mt_log,1,value))
     >   value = .true.

      control_hfxon_virtual = value
      return
      end


*     ***********************************
*     *                                 *
*     *     control_gradient_virtual    *
*     *                                 *
*     ***********************************
      logical function control_gradient_virtual()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:gradient_virtual',mt_log,1,value))
     >   value = .false.

      control_gradient_virtual = value
      return
      end

*     ***********************************
*     *                                 *
*     *  control_gradient_virtual_fac   *
*     *                                 *
*     ***********************************
      integer function control_gradient_virtual_fac()
      implicit none 
      
#include "bafdecls.fh"
#include "btdb.fh"
            
*     **** control_rtdb common block ****   
      integer rtdb
      common / control_rtdb1 / rtdb
            
      integer fac
              
      fac = 19
      if (.not.btdb_get(rtdb,'nwpw:gradient_virtual_fac',mt_int,1,fac))
     >   fac = 19
            
      control_gradient_virtual_fac = fac
      return   
      end

      
*     ***********************************
*     *                                 *
*     *     control_bound_virtual       *
*     *                                 *
*     ***********************************
      logical function control_bound_virtual()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
            
      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:bound_virtual',mt_log,1,value))
     >   value = .false.

      control_bound_virtual = value
      return
      end   

*     ***********************************
*     *                                 *
*     *       control_2x2_virtual       *
*     *                                 *
*     ***********************************
      logical function control_2x2_virtual()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
            
      logical value
      
      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:2x2_virtual',mt_log,1,value))
     >   value = .false.
      
      control_2x2_virtual = value
      return
      end

*     ***********************************
*     *                                 *
*     *       control_2x2ne_virtual     *
*     *                                 *
*     ***********************************
      logical function control_2x2ne_virtual()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:2x2ne_virtual',mt_log,1,value))
     >   value = .false.

      control_2x2ne_virtual = value
      return
      end



*     ***********************************
*     *                                 *
*     *       control_3x3_virtual       *
*     *                                 *
*     ***********************************
      logical function control_3x3_virtual()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:3x3_virtual',mt_log,1,value))
     >   value = .false.

      control_3x3_virtual = value
      return
      end


*     ***********************************
*     *                                 *
*     *       control_3x3ne_virtual     *
*     *                                 *
*     ***********************************
      logical function control_3x3ne_virtual()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:3x3ne_virtual',mt_log,1,value))
     >   value = .false.

      control_3x3ne_virtual = value
      return
      end


*     ***********************************
*     *                                 *
*     *       control_4x4_virtual       *
*     *                                 *
*     ***********************************
      logical function control_4x4_virtual()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
            
      logical value
      
      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:4x4_virtual',mt_log,1,value))
     >   value = .false.
      
      control_4x4_virtual = value
      return
      end


*     ***********************************
*     *                                 *
*     *       control_CI_virtual        *
*     *                                 *
*     ***********************************
      integer function control_CI_virtual()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
            
      integer ivalue

      logical  control_2x2_virtual,control_3x3_virtual
      external control_2x2_virtual,control_3x3_virtual
      logical  control_4x4_virtual
      external control_4x4_virtual
      logical  control_2x2ne_virtual,control_3x3ne_virtual
      external control_2x2ne_virtual,control_3x3ne_virtual
      
      ivalue = 0
      if (.not.btdb_get(rtdb,'nwpw:CI_virtual',mt_int,1,ivalue)) then
         if (control_2x2_virtual()) then
            ivalue = 2
         else if (control_3x3_virtual()) then
            ivalue = 3
         else if (control_4x4_virtual()) then
            ivalue = 4
         else if (control_2x2ne_virtual()) then
            ivalue = 5
         else if (control_3x3ne_virtual()) then
            ivalue = 6
         end if
      end if
      
      control_CI_virtual = ivalue
      return
      end





*     ***********************************
*     *                                 *
*     *     control_confine_virtual     *
*     *                                 *   
*     ***********************************
      logical function control_confine_virtual()
      implicit none
            
#include "bafdecls.fh"
#include "btdb.fh"
            
*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
               
      logical value
                 
      value = .false.
      if (.not.btdb_get(rtdb,'nwpw:confine_virtual',mt_log,1,value))
     >   value = .false.
                 
      control_confine_virtual = value
      return
      end   
              
              
           





*     **********************************************
*     *                                            *
*     *      control_init_velocities_temperature   *
*     *                                            *
*     **********************************************
      real*8 function control_init_velocities_temperature()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:init_velocities_temperature',
     >                  mt_dbl,1,value))
     >   value = 300.0d0

      control_init_velocities_temperature = value
      return
      end

*     **********************************************
*     *                                            *
*     *      control_init_velocities_seed          *
*     *                                            *
*     **********************************************
      integer function control_init_velocities_seed()
      implicit none
         
#include "bafdecls.fh"
#include "btdb.fh"
         
*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
            
      integer seed
            
      if (.not.btdb_get(rtdb,'nwpw:init_velocities_seed',
     >                  mt_int,1,seed)) 
     >   seed = 494
            
      control_init_velocities_seed = seed
      return
      end      




*     ***********************************
*     *                                 *
*     *      control_wannier_timestep   *
*     *                                 *
*     ***********************************
      real*8 function control_wannier_timestep()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'wannier:time_step',mt_dbl,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer value

      if (.not.btdb_get(rtdb,'wannier:maxiter',mt_int,1,value))
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "beom.fh"
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
      if (.not.btdb_get(rtdb,'nwpw:psp:cutoff',mt_dbl,1,wcut_max)) then 
         call Parallel_taskid(taskid)
         wcut_max = -1.0d0
         value = beom_create(geom,'geometry')
         value = value.and.beom_rtdb_load(rtdb,geom,'geometry')
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
              if (.not.btdb_get(rtdb,rtdbname(1:l),mt_dbl,1,wcut_atom)) 
     >         then

                 !*** define wcut_atom from psplibrary ***
                 value = btdb_parallel(.false.)
                 if (taskid.eq.MASTER) then
                    call util_directory_name(sdir_name,.true.,0)
                    h = index(sdir_name,' ') - 1
                    open(unit=99,file=sdir_name(1:h)//'/junk.inp',
     >                   status='unknown')
                   close(unit=99,status='delete')

                   call nwpw_libgeninp(1,aname,
     >                  sdir_name(1:h)//'/junk.inp')
                   wcut_atom=nwpw_libgcutoff(aname)

                   if(.not.btdb_put(rtdb,rtdbname(1:l),
     >                              mt_dbl,1,wcut_atom))
     >             call errquit(
     >           'control_ecut_wcut_default:cannot write wcut_atom',0,0)

                 end if
                 value = btdb_parallel(.true.)
                 call Parallel_Brdcst_value(MASTER,wcut_atom)

              end if

              !*** reset wcut_max ***
              if (wcut_atom.gt.wcut_max) then
                 wcut_max = wcut_atom
                 value = btdb_parallel(.false.)
                 if (taskid.eq.MASTER) then
                 if (.not.btdb_put(rtdb,'nwpw:psp:cutoff',
     >                             mt_dbl,1,wcut_max)) 
     >           call errquit(
     >          'control_ecut_wcut_default:cannot write wcut_max',0,0)
                 end if
                 value = btdb_parallel(.true.)
              end if
           end if
         end do
         if (.not. beom_destroy(geom))
     >   call errquit('control_ecut_wcut_default:cannot destroy geom',
     >                0,GEOM_ERR)
      end if

      wcut = wcut_max
      ecut = 2.0d0*wcut

      return
      end



*     ***********************************
*     *                                 *
*     *        control_pspspin          *
*     *                                 *
*     ***********************************
      logical function control_pspspin()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:pspspin',mt_log,1,value))
     >   value = .false.

      control_pspspin = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_set_pspspin      *
*     *                                 *
*     ***********************************
      subroutine control_set_pspspin(nion,
     >                               upscale,downscale,
     >                               upl,downl,
     >                               upm,downm,
     >                               upions,downions)
      implicit none
      integer nion
      real*8 upscale(*),downscale(*)
      integer upl(*),downl(*),upm(*),downm(*)
      logical upions(*),downions(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer taskid,MASTER
      parameter (MASTER=0)
      logical iamup
      integer i,ma_type,nactive_atoms,h_actlist,l_actlist
      integer pcount,ip,tl,tm
      real*8  tscale
      character*50 rtdb_name

*     *** external functions ***
      character*7 c_index_name
      external    c_index_name

      call Parallel_taskid(taskid)
      if (.not.btdb_get(rtdb,'nwpw:pspspin_count',mt_int,1,pcount))
     >   pcount = 0

      if (taskid.eq.MASTER) write(*,2300) 
      do ip=1,pcount

         rtdb_name = 'nwpw:pspspin_iamup:'//c_index_name(ip)
         if (.not.btdb_get(rtdb,rtdb_name,mt_log,1,iamup)) iamup=.true.

         if (iamup) then
            rtdb_name = 'nwpw:pspspin_upscale:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,tscale)) 
     >         tscale = 1.0d0
            rtdb_name = 'nwpw:pspspin_upl:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,tl))
     >         tl = -1
            rtdb_name = 'nwpw:pspspin_upm:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,tm))
     >         tm = 99999
            if (taskid.eq.MASTER) then
               if (tm.lt.999) then
                  write(*,2305) tl,tm,tscale
               else
                  write(*,2301) tl,tscale
               end if
            end if

            rtdb_name = 'nwpw:pspspin_upions:'//c_index_name(ip)
            if (rtdb_ma_get(rtdb,rtdb_name, ma_type,
     >                nactive_atoms, h_actlist)) then
               if (.not.BA_get_index(h_actlist,l_actlist))
     >            call errquit(
     >             'control_set_pspspin: ma_get_index failed',0,
     >             MA_ERR)
         
               if (taskid.eq.MASTER) 
     >         write(*,2302) (int_mb(l_actlist+i),i=0,nactive_atoms-1)
         
               do i=1,nactive_atoms
                  upions(int_mb(l_actlist+i-1))  = .true.
                  upl(int_mb(l_actlist+i-1))     = tl
                  upm(int_mb(l_actlist+i-1))     = tm
                  upscale(int_mb(l_actlist+i-1)) = tscale
               end do
               if (.not. BA_free_heap(h_actlist))
     >         call errquit(
     >         'control_set_pspspin:error freeing heap memory',0,MA_ERR)
            end if


         else
            rtdb_name = 'nwpw:pspspin_downscale:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,tscale))
     >         tscale = 1.0d0
            rtdb_name = 'nwpw:pspspin_downl:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,tl))
     >         tl = -1
            rtdb_name = 'nwpw:pspspin_downm:'//c_index_name(ip)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,tm))
     >         tm = 99999
            if (taskid.eq.MASTER) then
               if (tm.lt.999) then
                  write(*,2306) tl,tm,tscale
               else
                  write(*,2303) tl,tscale
               end if
            end if

            rtdb_name = 'nwpw:pspspin_downions:'//c_index_name(ip)
            if (rtdb_ma_get(rtdb,rtdb_name,ma_type,
     >                      nactive_atoms,h_actlist)) then
               if (.not.BA_get_index(h_actlist,l_actlist))
     >            call errquit(
     >             'control_set_pspspin: ma_get_index failed',1,
     >             MA_ERR)

               if (taskid.eq.MASTER) 
     >         write(*,2304) (int_mb(l_actlist+i),i=0,nactive_atoms-1)

               do i=1,nactive_atoms
                  downions(int_mb(l_actlist+i-1))  = .true.
                  downl(int_mb(l_actlist+i-1))     = tl
                  downm(int_mb(l_actlist+i-1))     = tm
                  downscale(int_mb(l_actlist+i-1)) = tscale
               end do
               if (.not.BA_free_heap(h_actlist))
     >         call errquit(
     >         'control_set_pspspin:error freeing heap memory',1,MA_ERR)
            end if
         end if
      end do

 2300 format(/1x,"Antiferromagnetic Pentalty Function Input:")
 2301 format(2x," - pspspin: up    l =",I2,10x," scale =",F8.3)
 2302 format(2x," - pspspin: up    ion indexes  =",10I5)
 2303 format(2x," - pspspin: down  l =",I2,10x," scale =",F8.3)
 2304 format(2x," - pspspin: down  ion indexes  =",10I5)
 2305 format(2x," - pspspin: up    l =",I2," not_m=",I3," scale =",F8.3)
 2306 format(2x," - pspspin: down  l =",I2," not_m=",I3," scale =",F8.3)

      return
      end

*     ***********************************
*     *                                 *
*     *        control_psputerm         *
*     *                                 *
*     ***********************************
      logical function control_psputerm()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:uterm',mt_log,1,value))
     >   value = .false.

      control_psputerm = value
      return
      end

*     ***********************************
*     *                                 *
*     *        control_pspnuterms       *
*     *                                 *
*     ***********************************
      integer function control_pspnuterms()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer nuterms

      nuterms = 0
      if (.not.btdb_get(rtdb,'nwpw:nuterms',mt_int,1,nuterms))
     >   nuterms = 0

      control_pspnuterms = nuterms
      return
      end


*     ***********************************
*     *                                 *
*     *        control_set_psputerm     *
*     *                                 *
*     ***********************************
      subroutine control_set_psputerm(nion,nuterms,l,uscale,jscale,ions)
      implicit none
      integer nion,nuterms,l(nuterms)
      real*8 uscale(nuterms),jscale(nuterms)
      logical ions(nion,nuterms)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer taskid,MASTER
      parameter (MASTER=0)
      integer nu,i,ma_type,nactive_atoms,h_actlist,l_actlist
      character*50 rtdb_name

*     **** external functions ****
      character*7 c_index_name
      external    c_index_name

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) write(*,3300) 

      do nu=1,nuterms

          rtdb_name = 'nwpw:uterm_scale:'//c_index_name(nu)
          if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,uscale(nu)))
     >      uscale(nu) = 0.0d0

          rtdb_name = 'nwpw:jterm_scale:'//c_index_name(nu)
          if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,1,jscale(nu)))
     >      jscale(nu) = 0.0d0

          rtdb_name = 'nwpw:uterm_l:'//c_index_name(nu)
          if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,l(nu)))
     >       l(nu) = -1

         do i=1,nion
            ions(i,nu) = .false.
         end do

         if (taskid.eq.MASTER) write(*,3301) l(nu),uscale(nu),jscale(nu)

          rtdb_name = 'nwpw:uterm_ions:'//c_index_name(nu)
         if (rtdb_ma_get(rtdb, rtdb_name,ma_type,
     >                   nactive_atoms, h_actlist)) then
            if (.not.BA_get_index(h_actlist,l_actlist))
     >         call errquit(
     >          'control_set_psputerm: ma_get_index failed',0,
     >          MA_ERR)
         
            if (taskid.eq.MASTER) then
               write(*,3302) (int_mb(l_actlist+i),i=0,nactive_atoms-1)
               write(*,*)
            end if
         
            do i=1,nactive_atoms
               ions(int_mb(l_actlist+i-1),nu) = .true.
            end do
            if (.not. BA_free_heap(h_actlist))
     >       call errquit(
     >        'control_set_psputerm:error freeing heap memory',0,MA_ERR)
         end if

      end do

 3300 format(/1x,"Hubbard Uterm Function Input:")
 3301 format(2x," - uterm: l =",I2,"  U=",F8.3," J=",F8.3)
 3302 format(2x," - uterm: ion indexes  =",10I5)
      return
      end



*     ***********************************
*     *                                 *
*     *        control_use_grid_cmp     *
*     *                                 *
*     ***********************************
      logical function control_use_grid_cmp()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:use_grid_cmp',mt_log,1,value))
     >   value = .false.

      control_use_grid_cmp = value
      return
      end


*     ***********************************
*     *                                 *
*     *  control_wgc_alphabetalambda    *
*     *                                 *
*     ***********************************
      real*8 function control_wgc_alphabetalambda(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (i.eq.1) then
         if (.not.btdb_get(rtdb,'nwpw:wgc_alpha',mt_dbl,1,value))
     >      value = 5.0d0/6.0d0
      else if (i.eq.2) then
         if (.not.btdb_get(rtdb,'nwpw:wgc_beta',mt_dbl,1,value))
     >      value = 5.0d0/6.0d0
      else 
         if (.not.btdb_get(rtdb,'nwpw:wgc_lambda',mt_dbl,1,value))
     >      value = 1.0d0
      end if

      control_wgc_alphabetalambda = value
      return
      end



*     ***********************************
*     *                                 *
*     *        control_mc_step_size     *
*     *                                 *
*     ***********************************
      real*8 function control_mc_step_size()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_step_size',mt_dbl,1,value))
     >  value = 0.50d0

      control_mc_step_size = value
      return
      end




*     ***********************************
*     *                                 *
*     *     control_mc_volume_step      *
*     *                                 *
*     ***********************************
      real*8 function control_mc_volume_step()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_volume_step',mt_dbl,1,value))
     >  value = 0.10d0

      control_mc_volume_step = value
      return
      end




*     ***********************************
*     *                                 *
*     *        control_mc_aratio        *
*     *                                 *
*     ***********************************
      real*8 function control_mc_aratio()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_aratio',mt_dbl,1,value))
     >  value = 0.234d0

      control_mc_aratio = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_mc_Temperature   *
*     *                                 *
*     ***********************************
      real*8 function control_mc_Temperature()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_temperature',mt_dbl,1,value))
     >  value = 298.15d0

      control_mc_Temperature = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_mc_pressure      *
*     *                                 *
*     ***********************************
      real*8 function control_mc_pressure()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      real*8 autoMbar,autoGPa,autoatm
      parameter (autoMbar=294.214239071d0)
      parameter (autoGPa=autoMbar*100.0d0)
      parameter (autoatm =290.360032539d6)


*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_pressure',mt_dbl,1,value))
     >  value = 1.0d0/autoatm

      control_mc_pressure = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_mc_ddx           *
*     *                                 *
*     ***********************************
      real*8 function control_mc_ddx()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_ddx',mt_dbl,1,value))
     >  value = 0.0d0

      control_mc_ddx = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_mc_ddv           *
*     *                                 *
*     ***********************************
      real*8 function control_mc_ddv()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      real*8 value

      if (.not.btdb_get(rtdb,'nwpw:mc_ddv',mt_dbl,1,value))
     >  value = 0.0d0

      control_mc_ddv = value
      return
      end



*     ***********************************
*     *                                 *
*     *        control_mc_seed          *
*     *                                 *
*     ***********************************
      integer function control_mc_seed()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer ivalue

      if (.not.btdb_get(rtdb,'nwpw:mc_seed',mt_int,1,ivalue))
     >  ivalue =  9484943

      control_mc_seed = ivalue
      return
      end


*     ***********************************
*     *                                 *
*     *      control_mc_algorithm       *
*     *                                 *
*     ***********************************
      integer function control_mc_algorithm()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer ivalue

      if (.not.btdb_get(rtdb,'nwpw:mc_algorithm',mt_int,1,ivalue))
     >  ivalue =  1

      control_mc_algorithm = ivalue
      return
      end




*     ***********************************
*     *                                 *
*     *    control_mc_atom_direction    *
*     *                                 *
*     ***********************************
      subroutine control_mc_atom_direction(mc_atom_direction)
      implicit none
      real*8 mc_atom_direction(3)

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'nwpw:mc_atom_direction',
     >                  mt_dbl,3,mc_atom_direction)) then
         mc_atom_direction(1) = 1.0d0
         mc_atom_direction(2) = 1.0d0
         mc_atom_direction(3) = 1.0d0
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *       control_mc_ngroups        *
*     *                                 *
*     ***********************************
      subroutine control_mc_ngroups(mc_napply,mc_ngroups,mc_group_size)
      implicit none
      integer mc_napply,mc_ngroups,mc_group_size

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_get(rtdb,'nwpw:mc_napply',
     >                  mt_int,1,mc_napply)) then
         mc_napply = 1
      end if
      if (.not.btdb_get(rtdb,'nwpw:mc_ngroups',
     >                  mt_int,1,mc_ngroups)) then
         mc_ngroups = 0
      end if
      if (.not.btdb_get(rtdb,'nwpw:mc_group_size',
     >                  mt_int,1,mc_group_size)) then
         mc_group_size = 0
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *       control_mc_groups        *
*     *                                 *
*     ***********************************
      subroutine control_mc_groups(mc_group_start,
     >                             mc_group_end,
     >                             mc_group)
      implicit none
      integer mc_group_start(*)
      integer mc_group_end(*)
      integer mc_group(*)

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer ng,ngs

      if (.not.btdb_get(rtdb,'nwpw:mc_ngroups',mt_int,1,ng)) ng=1
      if (.not.btdb_get(rtdb,'nwpw:mc_group_size',mt_int,1,ngs)) ngs=1

      if (.not.btdb_get(rtdb,'nwpw:mc_group_start',
     >                  mt_int,ng,mc_group_start)) then
         call icopy(ng,0,1,mc_group_start,1)
      end if
      if (.not.btdb_get(rtdb,'nwpw:mc_group_end',
     >                  mt_int,ng,mc_group_end)) then
         call icopy(ng,0,1,mc_group_end,1)
      end if
      if (.not.btdb_get(rtdb,'nwpw:mc_group',
     >                  mt_int,ngs,mc_group)) then
         call icopy(ngs,0,1,mc_group,1)
      end if

      return
      end


*     ************************************
*     *                                  *
*     *         control_hess_model       *
*     *                                  *
*     ************************************
      logical function control_hess_model()
      implicit none

#include "control.fh"

      control_hess_model = hess_model
      return
      end

*     ************************************
*     *                                  *
*     *      control_hess_filename       *
*     *                                  *
*     ************************************
      subroutine control_hess_filename(filehess)
      implicit none
      character*(*) filehess

#include "control.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_cget(rtdb,'nwpw:hess_model:filename',
     >                   1,filehess)) then
         call util_file_name('hess',.false.,.false.,filehess)
      end if
      return
      end
      


*     ***********************************
*     *                                 *
*     *    control_psp_semicore_small   *
*     *                                 *
*     ***********************************
      logical function control_psp_semicore_small()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:psp:semicore_small',mt_log,1,value))
     >   value = .false.

      control_psp_semicore_small = value
      return
      end

*     ***********************************
*     *                                 *
*     *         control_incell          *
*     *                                 *
*     ***********************************
      logical function control_incell()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:incell',mt_log,1,value))
     >   value = .true.

      control_incell = value
      return
      end


*     ***********************************
*     *                                 *
*     *    control_fetch_confine        *
*     *                                 *
*     ***********************************
      subroutine control_fetch_confine(symbol,q1,rs,re)
      implicit none
      character*4 symbol
      real*8 q1,rs,re

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value
      integer l
      real*8 rrr(3)

      l = index(symbol,' ') - 1
      if (.not.btdb_get(rtdb,'nwpw:confine:'//symbol(1:l),
     >                  mt_dbl,3,rrr)) then

         rrr(1) = 1.0d0
         rrr(2) = 6.0d0
         rrr(3) = 14.0d0
      end if
      q1 = rrr(1)
      rs = rrr(2)
      re = rrr(3)
      return
      end


*     ***********************************
*     *                                 *
*     *       control_runsocket         *
*     *                                 *
*     ***********************************
      logical function control_runsocket()
      implicit none

#include "inp.fh"

*     **** local variables ****
      logical value
      character*40 stype

      value = .false.
      call control_socket_type(stype)
      if (inp_compare(.false.,'ipi_client',stype)) value = .true.
      if (inp_compare(.false.,'unix',stype)) value = .true.

      control_runsocket = value
      return
      end


*     ************************************
*     *                                  *
*     *      control_socket_type         *
*     *                                  *
*     ************************************
      subroutine control_socket_type(stype)
      implicit none
      character*(*) stype

#include "control.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_cget(rtdb,'nwpw:socket_type',
     >                   1,stype)) then
         stype = "no socket"
      end if
      return
      end

*     ************************************
*     *                                  *
*     *      control_socket_ip           *
*     *                                  *
*     ************************************
      subroutine control_socket_ip(socket_ip)
      implicit none
      character*(*) socket_ip
      character*40 stype

#include "control.fh"
#include "btdb.fh"
#include "inp.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.btdb_cget(rtdb,'nwpw:socket_ip',
     >                   1,socket_ip)) then
         call control_socket_type(stype)
         if (inp_compare(.false.,'unix',stype)) then
            socket_ip = "nwchem"
         else
            socket_ip = "127.0.0.1:31415"
         end if
      end if
      return
      end



*     ************************************
*     *                                  *
*     *      control_vfield_filenames    *
*     *                                  *
*     ************************************
      logical function control_vfield_filenames(filenames)
      implicit none
      character*(*) filenames

#include "control.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical found

      found = .false.
      if (btdb_cget(rtdb,'nwpw:vfield_filenames',1,filenames)) then
         found = .true.
      end if

      control_vfield_filenames = found
      return
      end

*     ***********************************
*     *                                 *
*     *        control_mult_fixed       *
*     *                                 *
*     ***********************************
      logical function control_mult_fixed()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:mult_fixed',mt_log,1,value))
     >   value = .false.

      control_mult_fixed = value
      return
      end


*     ***********************************
*     *                                 *
*     *        control_gas_energy       *
*     *                                 *
*     ***********************************
      real*8 function control_gas_energy()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** external functions ****
      logical  nwpw_cosmo_on, nwpw_born_on
      external nwpw_cosmo_on, nwpw_born_on

      real*8 egas

      egas = 0.0d0
      if (nwpw_cosmo_on().or.nwpw_born_on()) then
         if (.not.btdb_get(rtdb,'nwpw:gas_energy',mt_dbl,1,egas))
     >      egas = 0.0d0
      end if

      control_gas_energy = egas
      return
      end
      

*     ***********************************
*     *                                 *
*     *      control_gas_energy_set     *
*     *                                 *
*     ***********************************
      logical function control_gas_energy_set(egas)
      implicit none
      real*8 egas

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

*     **** extertnal functions ****
      logical  nwpw_cosmo_on,nwpw_born_on
      external nwpw_cosmo_on,nwpw_born_on
      logical ok

      ok=.true.
      if (.not.(nwpw_cosmo_on().or.nwpw_born_on())) then
         ok=ok.and.btdb_put(rtdb,'nwpw:gas_energy',mt_dbl,1,egas)
      end if
      control_gas_energy_set=ok
      return
      end


*     ***********************************
*     *                                 *
*     *        control_H1_it_in          *
*     *                                 *
*     ***********************************
      integer function control_H1_it_in()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer it_in

      if (.not.btdb_get(rtdb,'nwpw:H1_it_in',mt_int,1,it_in)) 
     >   it_in = 50

      control_H1_it_in = it_in
      return
      end

*     ***********************************
*     *                                 *
*     *        control_H1_it_out        *
*     *                                 *
*     ***********************************
      integer function control_H1_it_out()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer it_out

      if (.not.btdb_get(rtdb,'nwpw:H1_it_out',mt_int,1,it_out)) 
     >   it_out = 1

      control_H1_it_out = it_out
      return
      end


*     ***********************************
*     *                                 *
*     *      control_H1_it_ortho        *
*     *                                 *
*     ***********************************
      integer function control_H1_it_ortho()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer it_ortho

      if (.not.btdb_get(rtdb,'nwpw:H1_it_ortho',mt_int,1,it_ortho))
     >   it_ortho = -1

      control_H1_it_ortho = it_ortho
      return
      end

*     ***********************************
*     *                                 *
*     *       control_CI_maxit_sweeps   *
*     *                                 *
*     ***********************************
      integer function control_CI_maxit_sweeps()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb
      
      integer it_out
      
      if (.not.btdb_get(rtdb,'nwpw:CI_maxit_sweeps',mt_int,1,it_out))
     >   it_out = 3
      
      control_CI_maxit_sweeps = it_out
      return
      end


*     ***********************************
*     *                                 *
*     *       control_CI_maxit_orb      *
*     *                                 *
*     ***********************************
      integer function control_CI_maxit_orb()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer it_out

      if (.not.btdb_get(rtdb,'nwpw:CI_maxit_orb',mt_int,1,it_out))
     >   it_out = 120

      control_CI_maxit_orb = it_out
      return
      end

*     ***********************************
*     *                                 *
*     *   control_CI_maxit_linesearch   *
*     *                                 *
*     ***********************************
      integer function control_CI_maxit_linesearch()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer it_out

      if (.not.btdb_get(rtdb,'nwpw:CI_maxit_linesearch',
     >                  mt_int,1,it_out))
     >   it_out = 15

      control_CI_maxit_linesearch = it_out
      return
      end


*     ***********************************
*     *                                 *
*     *        control_use_director     *
*     *                                 *
*     ***********************************
      logical function control_use_director()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:use_director',mt_log,1,value))
     >   value = .false.

      control_use_director = value
      return
      end


*     ***********************************
*     *                                 *
*     *      control_unset_excited_ne   *
*     *                                 *
*     ***********************************
      subroutine control_unset_excited_ne()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value
      integer ne(2)

      ne(1) = 0
      ne(2) = 0
      value = btdb_put(rtdb,'nwpw:excited_ne',mt_int,2,ne)

      !value = rtdb_delete(rtdb,'nwpw:excited_ne')

      return
      end


*     ***********************************
*     *                                 *
*     *        control_Efield_on        *
*     *                                 *
*     ***********************************
      logical function control_Efield_on()
      implicit none

#include "bafdecls.fh"
#include "btdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      logical value

      value = .true.
      if (.not.btdb_get(rtdb,'nwpw:efield',mt_log,1,value))
     >   value = .false.

      control_Efield_on = value
      return
      end


*     *********************************************
*     *                                           *
*     *        control_single_precision_on        *
*     *                                           *
*     *********************************************
      logical function control_single_precision_on()
      implicit none

#include "control.fh"

      control_single_precision_on = single_precision_on
      return
      end
