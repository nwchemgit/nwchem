*
* $Id: control.F,v 1.63 2004-12-06 20:03:24 bylaska Exp $
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
      integer ispin0,ne(2)

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
      integer i,l


      value = rtdb_parallel(.true.)
      code = code_in
      trtdb = rtdb

      call util_print_get_level(print_level)
      !write(*,*) "Print level is ",print_level


*     **** cpsd: stuff ****

      if (code.eq.1) then
      if (.not.rtdb_cget(rtdb,'cpsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:exchange_correlation',
     >                  1,exchange_correlation)

      value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.false.)
      value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.false.)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,1.0d0)
      if (inp_compare(.false.,exchange_correlation,'off')) then
        gga = -1
      else if (inp_compare(.false.,exchange_correlation,'vosko')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lsda')) then
        gga = 0
        value = rtdb_put(rtdb,'cpsd:ispin',mt_int,1,2)
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
        gga = 1
      else if (inp_compare(.false.,exchange_correlation,'hf')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'hartree-fock')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,'pbe0')) then
        gga = 3
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,0.25d0)
      else if (inp_compare(.false.,exchange_correlation,'lda-sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'pbe96-sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                    'lda-0.4sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else if (inp_compare(.false.,exchange_correlation,
     >         'pbe96-0.4sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else
        gga = 0
      end if

      value = value.and.
     >        rtdb_get(rtdb,'cpsd:geometry_optimize',mt_log,1,move)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:fractional_coordinates',
     >                 mt_log,1,frac_coord)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:npsp',mt_int,1,npsp)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:tolerances',mt_dbl,3,tolerances)
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:wcut',mt_dbl,1,wcut)

      value = value.and.
     >        rtdb_get(rtdb,'cpsd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:ncut',mt_int,1,ncut)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:mult',mt_int,1,multiplicity)
      value = value.and.
     >        rtdb_get(rtdb,'cpsd:ispin',mt_int,1,ispin)


*     **** cpmd code: stuff ****
      else if (code.eq.2) then
      if (.not.rtdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename)

      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename)

      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:exchange_correlation',
     >                  1,exchange_correlation)
      value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.false.)
      value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.false.)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,1.00d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,1.00d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,1.0d0)
      if (inp_compare(.false.,exchange_correlation,'off')) then
        gga = -1
      else if (inp_compare(.false.,exchange_correlation,'vosko')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lsda')) then
        gga = 0
        value = rtdb_put(rtdb,'cpsd:ispin',mt_int,1,2)
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
        gga = 1
      else if (inp_compare(.false.,exchange_correlation,'hf')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'hartree-fock')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,'pbe0')) then
        gga = 3
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,0.25d0)
      else if (inp_compare(.false.,exchange_correlation,'lda-sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'pbe96-sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                    'lda-0.4sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else if (inp_compare(.false.,exchange_correlation,
     >         'pbe96-0.4sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else
        gga = 0
      end if


      value = value.and.
     >        rtdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)
      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut)

      value = value.and.
     >        rtdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut)

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


c     do i=1,npsp
c       id = 'cpmd:psp'//CHAR(i)
c       value = value.and.
c    >          rtdb_cget(rtdb,id,1,pspnames(i))
c     end do

*     **** get thermostat information ****
      value = value .and. rtdb_get(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value .and. rtdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value .and. rtdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value .and. rtdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value .and. rtdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr)


*     **** cgsd: stuff or  paw: stuff ****
      else if ((code.eq.3).or.(code.eq.8)) then
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
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:input_ewavefunction_filename',
     >                  1,input_ewavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:output_ewavefunction_filename',
     >                  1,output_ewavefunction_filename)

      if (.not.rtdb_cget(rtdb,'cgsd:exchange_correlation',
     >                   1,exchange_correlation))
     >  exchange_correlation = 'vosko'

      value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.false.)
      value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.false.)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,1.0d0)
      value = value.and.
     >        rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,1.0d0)
      if (inp_compare(.false.,exchange_correlation,'off')) then
        gga = -1
      else if (inp_compare(.false.,exchange_correlation,'vosko')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'svwn5')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lda')) then
        gga = 0
      else if (inp_compare(.false.,exchange_correlation,'lsda')) then
        gga = 0
        value = rtdb_put(rtdb,'cpsd:ispin',mt_int,1,2)
      else if (inp_compare(.false.,exchange_correlation,'pbe96')) then
        gga = 1
      else if (inp_compare(.false.,exchange_correlation,'hf')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'hartree-fock')) then
        gga = 2
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,'pbe0')) then
        gga = 3
        value = value.and.rtdb_put(rtdb,'pspw:HFX',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:HFX_parameter',mt_dbl,1,0.25d0)
      else if (inp_compare(.false.,exchange_correlation,'lda-sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                     'pbe96-sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
      else if (inp_compare(.false.,exchange_correlation,
     >                    'lda-0.4sic')) then
        gga = 0
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else if (inp_compare(.false.,exchange_correlation,
     >         'pbe96-0.4sic')) then
        gga = 1
        value = value.and.rtdb_put(rtdb,'pspw:SIC_all',mt_log,1,.true.)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_xc_parameter',mt_dbl,1,0.40d0)
        value = value.and.
     >          rtdb_put(rtdb,'pspw:SIC_h_parameter',mt_dbl,1,0.40d0)
      else
        gga = 0
      end if


*     **** set Kohn-Sham scf parameters ***
      if (.not. rtdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not.rtdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 0

      if (.not.rtdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   scf_algorithm = 0

*     **** set maxit_orb maxit_orbs ***
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 0


      value = value.and.
     >        rtdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      value = value.and.
     >        rtdb_get(rtdb,'cgsd:fractional_coordinates',
     >                 mt_log,1,frac_coord)
  
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut)     
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin)

c     do i=1,npsp
c       id = 'cgsd:psp'//CHAR(i)
c       value = value.and.
c    >          rtdb_cget(rtdb,id,1,pspnames(i))
c     end do

*     **** pspw_dplot: stuff ****
      else if (code.eq.4) then
         value = .true.
         value = rtdb_cget(rtdb,'pspw_dplot:wavefunction_filename',
     >                  1,input_wavefunction_filename)
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
         ncut = 0
         npsp = 0

         control_read = value
         return

*     **** band: stuff ****

      else if (code.eq.5) then
      if (.not.rtdb_cget(rtdb,'band:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
cc      value = value.and.
cc     >        rtdb_cget(rtdb,'band:input_wavefunction_filename',
cc     >                  1,input_wavefunction_filename)
cc      value = value.and.
cc     >        rtdb_cget(rtdb,'band:output_wavefunction_filename',
cc     >                  1,output_wavefunction_filename)
cc
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
         gga = 1
      else
         gga = 0
      end if

      value = value.and.
     >        rtdb_get(rtdb,'band:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'band:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'band:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'band:tolerances',mt_dbl,3,tolerances)
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0

      value = value.and.
     >        rtdb_get(rtdb,'band:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'band:wcut',mt_dbl,1,wcut)
      value = value.and.
     >        rtdb_get(rtdb,'band:rcut',mt_dbl,1,rcut)     
      value = value.and.
     >        rtdb_get(rtdb,'band:ncut',mt_int,1,ncut)
      value = value.and.
     >        rtdb_get(rtdb,'band:mult',mt_int,1,multiplicity)
      value = value.and.
     >        rtdb_get(rtdb,'band:ispin',mt_int,1,ispin)

*     **** set Kohn-Sham scf parameters ***
      if (.not. rtdb_get(rtdb,'nwpw:ks_alpha',mt_dbl,1,ks_alpha))
     >   ks_alpha = 0.25d0
      if (.not.rtdb_get(rtdb,'nwpw:scf_algorithm',
     >                  mt_int,1,scf_algorithm))
     >   scf_algorithm = 0
      if (.not.rtdb_get(rtdb,'nwpw:ks_algorithm',
     >                  mt_int,1,ks_algorithm))
     >   ks_algorithm = 0

*     **** set maxit_orb maxit_orbs ***
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orb',mt_int,1,maxit_orb))
     >  maxit_orb = 5
      if (.not.rtdb_get(rtdb,
     >      'nwpw:ks_maxit_orbs',mt_int,1,maxit_orbs))
     >  maxit_orbs = 0


*     **** paw_sd: stuff ****
      else if (code.eq.6) then

      if (.not.rtdb_get(rtdb,'cgsd:geometry_optimize',mt_log,1,move))
     >   move = .false.

      if (.not.rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)

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
         gga = 1
      else
         gga = 0
      end if

      value = value.and.
     >        rtdb_get(rtdb,'cgsd:npsp',mt_int,1,npsp)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:tolerances',mt_dbl,3,tolerances)
      scaling(1) = 0.0d0
      scaling(2) = 0.0d0
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:wcut',mt_dbl,1,wcut)

      value = value.and.
     >        rtdb_get(rtdb,'cgsd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ncut',mt_int,1,ncut)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:mult',mt_int,1,multiplicity)
      value = value.and.
     >        rtdb_get(rtdb,'cgsd:ispin',mt_int,1,ispin)




*     **** paw_cpmd: stuff ****
      else if (code.eq.7) then

      if (.not.rtdb_cget(rtdb,'cpmd:cell_name',1,cell_name)) then
        cell_name = 'cell_default'
      end if
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:input_v_wavefunction_filename',
     >                  1,input_v_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:output_v_wavefunction_filename',
     >                  1,output_v_wavefunction_filename)

      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:xyz_filename',
     >                  1,xyz_filename)

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
         gga = 1
      else
         gga = 0
      end if


      value = value.and.
     >        rtdb_get(rtdb,'cpmd:geometry_optimize',mt_log,1,move)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:fractional_coordinates',
     >                 mt_log,1,frac_coord)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:npsp',mt_int,1,npsp)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:fake_mass',mt_dbl,1,fake_mass)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:time_step',mt_dbl,1,time_step)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:loop',mt_int,2,loop)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:scaling',mt_dbl,2,scaling)
      tolerances(1) = 0.0d0
      tolerances(2) = 0.0d0
      tolerances(3) = 0.0d0
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:ecut',mt_dbl,1,ecut)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:wcut',mt_dbl,1,wcut)

      value = value.and.
     >        rtdb_get(rtdb,'cpmd:rcut',mt_dbl,1,rcut)
      value = value.and.
     >        rtdb_get(rtdb,'cpmd:ncut',mt_int,1,ncut)

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

*     **** get thermostat information ****
      value = value .and. rtdb_get(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value .and. rtdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value .and. rtdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value .and. rtdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value .and. rtdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr)





*     **** pspw_wannier ****
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

c         value = rtdb_cget(rtdb,'wannier:input_wavefunction_filename',
c     >                  1,input_wavefunction_filename)
c         value = value.and.
c     >           rtdb_cget(rtdb,'wannier:output_wavefunction_filename',
c     >                  1,output_wavefunction_filename)

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
         ncut = 0
         npsp = 0

         control_read = value
         return


*     **** unknown code type ****
   
      else
         value = .false.
         write(*,*) "control_read: unknown code type:",code
         control_read = value
         return
      end if

*     ***** cell: stuff ****
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      rtdb_ngrid = cell_name(1:l)//':ngrid'
      rtdb_boundry = cell_name(1:l)//':boundry'

           
      if (.not.rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)) then
        call dcopy(9,0.0d0,0,unita,1)
      end if

      if (.not.rtdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)) then
        ngrid(1) = 32
        ngrid(2) = 32
        ngrid(3) = 32
      end if

      if (.not.rtdb_cget(rtdb,rtdb_boundry,1,boundry)) then
         boundry = 'periodic'
      end if

      call check_unita_for_default(rtdb,unita)


*     **** set out of time variables ****
      est_step_time   = -1
      est_finish_time = -1
      call current_second(cpu1_time)

*     **** set gram_schmidt ***
      gram_schmidt = .false.
      if (.not.rtdb_get(rtdb,
     >      'nwpw:gram-schmidt',mt_log,1,gram_schmidt))
     >  gram_schmidt = .false.

      
      control_read = value
      return
      end


*     ***********************************
*     *					*
*     *    check_unita_for_default	*
*     *					*
*     ***********************************
      subroutine check_unita_for_default(rtdb,unita)
      implicit none
#include "errquit.fh"
      integer rtdb
      real*8 unita(3,3)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer geom

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
           unita(1,1) = 20.0d0
           unita(2,1) =  0.0d0
           unita(3,1) =  0.0d0
           unita(1,2) =  0.0d0
           unita(2,2) = 20.0d0
           unita(3,2) =  0.0d0
           unita(1,3) =  0.0d0
           unita(2,3) =  0.0d0
           unita(3,3) = 20.0d0
         end if
         value = geom_destroy(geom)
         if (.not. value) call errquit('cannot destroy geom',0,
     &       GEOM_ERR)

      end if

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

#include "control.fh"

*     **** local variables ****
      integer l,version

      l =index(boundry,' ') - 1

      version = 3
      if (boundry(1:l).eq.'periodic')  version=3
      if (boundry(1:l).eq.'aperiodic') version=4

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

      if (code.eq.1) then
        if (.not.rtdb_get(rtdb,'cpsd:mulliken',mt_log,1,value))
     >        call errquit(
     >        'control_mulliken: rtdb_get mulliken failed', 0, RTDB_ERR)
      end if
      if ((code.eq.2).or.(code.eq.7)) then
        if (.not.rtdb_get(rtdb,'cpmd:mulliken',mt_log,1,value))
     >        call errquit(
     >        'control_mulliken: rtdb_get mulliken failed', 0, RTDB_ERR)
      end if


      if (code.eq.3) then
        if (.not.rtdb_get(rtdb,'cgsd:mulliken',mt_log,1,value))
     >        call errquit(
     >        'control_mulliken: rtdb_get mulliken failed', 0, RTDB_ERR)
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

      control_allow_translation = value
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

      value = rtdb_cget(rtdb,'band:zone_name',1,zone_name)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':number_kvectors'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors)

      if (.not. value)
     >  call errquit('control_num_kvectors: failed', 0, RTDB_ERR)

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
#include "errquit.fh"
      integer i
      real*8 ks(4)

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
      integer kvs(2)

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      num_kvectors = control_num_kvectors()
      value = MA_push_get(mt_dbl,(4*num_kvectors),
     >        'kvs',kvs(2),kvs(1))
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0,
     &       MA_ERR)

      value = value.and.
     >        rtdb_cget(rtdb,'band:zone_name',1,zone_name)
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
      character*2 atom
      integer     lmax,locp
      real*8      rlocal

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer l,k
      character*4  element
      character*20 rtdb_name

      element = '    '
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


*     *****************************
*     *                           *
*     *    control_mapping        *
*     *                           *
*     *****************************
      integer function control_mapping()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      integer mapping

      if (.not.rtdb_get(rtdb,'nwpw:mapping',mt_int,1,mapping)) then
         mapping = 1
      end if

      control_mapping = mapping
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

