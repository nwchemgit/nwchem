*
* $Id: control.F,v 1.11 2001-12-14 20:02:46 bylaska Exp $
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


      character*30 cell_name,id
      character*30 rtdb_unita,rtdb_ngrid,rtdb_boundry
      character*30 exchange_correlation
      integer i,l

      value = rtdb_parallel(.true.)
      code = code_in
      trtdb = rtdb

*     **** cpsd: stuff ****

      if (code.eq.1) then
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:cell_name',1,cell_name)
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cpsd:exchange_correlation',
     >                  1,exchange_correlation)
      if (exchange_correlation.eq.'vosko') then
         gga = 0
      else
         gga = 1
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

      do i=1,npsp
        id = 'cpsd:psp'//CHAR(i)
        value = value.and.
     >          rtdb_cget(rtdb,id,1,pspnames(i))
      end do

*     **** cpmd code: stuff ****
      else if (code.eq.2) then
      value = value.and.
     >        rtdb_cget(rtdb,'cpmd:cell_name',1,cell_name)
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
      if (exchange_correlation.eq.'vosko') then
         gga = 0
      else
         gga = 1
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


      do i=1,npsp
        id = 'cpmd:psp'//CHAR(i)
        value = value.and.
     >          rtdb_cget(rtdb,id,1,pspnames(i))
      end do

*     **** get thermostat information ****
      value = value .and. rtdb_get(rtdb,'cpmd:nose',mt_log,1,nose)
      value = value .and. rtdb_get(rtdb,'cpmd:Pe',mt_dbl,1,Pe)
      value = value .and. rtdb_get(rtdb,'cpmd:Te',mt_dbl,1,Te)
      value = value .and. rtdb_get(rtdb,'cpmd:Pr',mt_dbl,1,Pr)
      value = value .and. rtdb_get(rtdb,'cpmd:Tr',mt_dbl,1,Tr)


*     **** cgsd: stuff ****

      else if (code.eq.3) then
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:exchange_correlation',
     >                  1,exchange_correlation)
      if (exchange_correlation.eq.'vosko') then
         gga = 0
      else
         gga = 1
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

      do i=1,npsp
        id = 'cgsd:psp'//CHAR(i)
        value = value.and.
     >          rtdb_cget(rtdb,id,1,pspnames(i))
      end do

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
      value = value.and.
     >        rtdb_cget(rtdb,'band:cell_name',1,cell_name)
      value = value.and.
     >        rtdb_cget(rtdb,'band:input_wavefunction_filename',
     >                  1,input_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'band:output_wavefunction_filename',
     >                  1,output_wavefunction_filename)
      value = value.and.
     >        rtdb_cget(rtdb,'band:exchange_correlation',
     >                  1,exchange_correlation)
      if (exchange_correlation.eq.'vosko') then
         gga = 0
      else
         gga = 1
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

           
      value = value.and.
     >        rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      value = value.and.
     >        rtdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)
      value = value.and.
     >        rtdb_cget(rtdb,rtdb_boundry,1,boundry)

*     **** check unita for defaults ****
      call check_unita_for_default(rtdb,unita)

      control_read = value
      return
      end

*     ***********************************
*     *									*
*     *		check_unita_for_default		*
*     *									*
*     ***********************************
      subroutine check_unita_for_default(rtdb,unita)
      implicit none
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
         if (.not. value) call errquit('cannot load geometry',0)

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
         if (.not. value) call errquit('cannot destroy geom',0)

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


*     ***********************************
*     *					*
*     *		control_pspnames	*
*     *					*
*     ***********************************
      character*20  function control_pspnames(i)
      implicit none 
      integer i

#include "control.fh"

      control_pspnames = pspnames(i)
      return
      end

*     ***********************************
*     *							 		*
*     *		control_pspstressnames		*
*     *									*
*     ***********************************
      character*20  function control_pspstressnames(i)
      implicit none 
      integer i

      integer ind
      character*20 pspname
      character*20 control_pspnames
      external     control_pspnames

      pspname = control_pspnames(i)
      ind = index(pspname,' ') -1
      pspname = pspname(1:ind)//'2'

      control_pspstressnames = pspname
      return
      end

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
*     *									*
*     *		control_multiplicity		*
*     *									*
*     ***********************************
      integer function control_multiplicity()
      implicit none

#include "control.fh"

      control_multiplicity = multiplicity
      return 
      end


*     ***********************************
*     *									*
*     *			control_ispin      		*
*     *									*
*     ***********************************
      integer function control_ispin()
      implicit none

#include "control.fh"

      control_ispin = ispin
      return 
      end


*     *******************************************
*     *											*
*     *		control_gradient_iterations		*
*     *											*
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
*     *                		 *
*     *     control_Nose     *
*     *                 	 *
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
*     *                	 		 *
*     *     control_Nose_Pe      *
*     *                 	 	 *
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
*     *                	 		 *
*     *     control_Nose_Pr      *
*     *                 	 	 *
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
*     *                	 		 *
*     *     control_Nose_Tr      *
*     *                 	 	 *
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
*     *                	 		 *
*     *     control_Mulliken     *
*     *                 	 	 *
*     ****************************
      logical function control_Mulliken()
      implicit none

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
     >        'control_mulliken: rtdb_get mulliken failed', 0)
      end if
      if (code.eq.2) then
        if (.not.rtdb_get(rtdb,'cpmd:mulliken',mt_log,1,value))
     >        call errquit(
     >        'control_mulliken: rtdb_get mulliken failed', 0)
      end if
      if (code.eq.3) then
        if (.not.rtdb_get(rtdb,'cgsd:mulliken',mt_log,1,value))
     >        call errquit(
     >        'control_mulliken: rtdb_get mulliken failed', 0)
      end if


      control_Mulliken = value
      return
      end


*     *****************************
*     *                	 		  *
*     * control_allow_translation *
*     *                 	 	  *
*     *****************************
      logical function control_allow_translation()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"

      logical value

*     **** control_rtdb common block ****
      integer rtdb
      common / control_rtdb1 / rtdb

      if (.not.rtdb_get(rtdb,'cgsd:allow_translation',
     >                  mt_log,1,value))
     >  call errquit(
     > 'control_allow_translation:rtdb_get allow_translation failed',0)

      control_allow_translation = value
      return
      end

*     ****************************
*     *                	 		 *
*     *  control_num_kvectors    *
*     *                 	 	 *
*     ****************************
      integer function control_num_kvectors()
      implicit none

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
     >  call errquit('control_num_kvectors: failed', 0)

      control_num_kvectors = num_kvectors
      return
      end

*     ****************************
*     *                	 		 *
*     *      control_ksvector	 *
*     *                 	 	 *
*     ****************************
      subroutine control_ksvector(i,ks)
      implicit none
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
     >  call errquit('control_ksvector: failed to get zone name', 0)

      value = value.and.
     >        rtdb_cget(rtdb,'band:zone_name',1,zone_name)
      if (.not. value)
     >  call errquit('control_ksvector: failed to get zone name', 0)

      l = index(zone_name,' ') -1
      rtdb_name = zone_name(1:l)//':kvectors'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1)))

      if (.not. value)
     >  call errquit('control_ksvector: failed to get kvs', 0)

      ks(1) = dbl_mb(kvs(1)+4*(i-1))
      ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
      ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
      ks(4) = dbl_mb(kvs(1)+4*(i-1)+3)

      value = value.and.MA_pop_stack(kvs(2))

      if (.not. value)
     >  call errquit('control_ksvector: failed to free stack', 0)
      return
      end

*     ****************************
*     *                	 		 *
*     *      control_kvector	 *
*     *                 	 	 *
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
*     *                	 		  *
*     *    control_TotalCharge	  *
*     *                 	 	  *
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
*     *                	 		  *
*     *    	   control_rtdb		  *
*     *                 	 	  *
*     *****************************
      integer function control_rtdb()
      implicit none

*     **** control_rtdb common block ****
      integer trtdb
      common / control_rtdb1 / trtdb

      control_rtdb = trtdb
      return
      end

