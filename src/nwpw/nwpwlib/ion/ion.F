*
* $Id: ion.F,v 1.8 2001-11-29 02:00:40 bylaska Exp $
*

*     ***************************
*     *							*
*     *		ion_end 			*
*     *							*
*     ***************************
      subroutine ion_end()
      implicit none

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      logical value

      value = MA_free_heap(r2(2))
      value = value.and.MA_free_heap(r1(2))
      value = value.and.MA_free_heap(r0(2))
      value = value.and.MA_free_heap(katm(2))
      value = value.and.MA_free_heap(natm(2))
      value = value.and.MA_free_heap(atom(2))
      value = value.and.MA_free_heap(dti(2))
      if (.not. value) call errquit('ion_end:free heap ',0)

      return
      end

*     ***************************
*     *				*
*     *		ion_write	*
*     *				*
*     ***************************
      subroutine ion_write(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

*     **** ion common block ****
#include "ion.fh"

      integer geom
      common / ION2/ geom

      logical mmexist
      common / pspwqmm/ mmexist

*     **** local variables ****
      integer i
      logical value
      double precision rxyz(3),q
      character*16     t

*     **** external functions ****
      logical  control_frac_coord
      integer  control_code
      external control_frac_coord
      external control_code

*     *** write out CIF file ***
      call CIF_write(rtdb)

      if (control_frac_coord()) call fcoord_to_frac(nion,dbl_mb(r1(1)))     

      do i=1,nion
         value = geom_cent_get(geom,i,t,rxyz,q)
         value = geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
      end do

*     **** if md code then write velocities ****
      if (control_code().eq.2) then
         value = geom_vel_set(geom,dbl_mb(r0(1)))
      end if

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
        value = geom_rtdb_delete(rtdb,'pspwgeometry')
      if ((control_code().eq.1).or.
     >    (control_code().eq.2)) then
        value = geom_rtdb_store(rtdb,geom,'pspwgeometry')
      end if
      value = geom_destroy(geom)

      if (control_frac_coord()) call fcoord_to_real(nion,dbl_mb(r1(1)))     

*     *****************************************************
*     **** update rtdb if md or steepest descent code  ****
*     *****************************************************
      if ((control_code().eq.1).or.
     >    (control_code().eq.2)) then


*       *****************************************
*       **** put together ions,charges       ****
*       ***** and molecular pseudopotentials ****
*       *****************************************
        call combine_pointcharge(rtdb)
        call combine_molpsp(rtdb)

      else
        value = geom_rtdb_delete(rtdb,'chargepspwgeometry')
        value = geom_rtdb_delete(rtdb,'pspwgeometry')
        value = geom_rtdb_delete(rtdb,'qmmmgeometry')
      end if
      return
      end

*     ***************************
*     *							*
*     *		ion_destroy 		*
*     *							*
*     ***************************
      subroutine ion_destroy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer geom
      common / ION2/ geom
      logical value

      value = geom_destroy(geom)
      return
      end
    



*     ***************************
*     *							*
*     *		ion_init			*
*     *							*
*     ***************************
*
*        This routine reads in the ion data structure from the
*     the runtime database
*
*     Entry - rtdb
*     Exit  - 
*
*     Uses - Parallel and MPI routines
*
      logical function ion_init(rtdb)
      implicit none 
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"


***** ion common block ****
#include "ion.fh"

      integer geom
      common / ION2/ geom

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables ***
      integer     i
      double precision q,dt,vgx,vgy,vgz
      character*16     t
      logical value

      integer MASTER,taskid
      parameter (MASTER=0)

*     **** external functions ****
      character*2 ion_aname
      integer     control_code
      external    ion_aname
      external    control_code
      
      logical          control_frac_coord
      double precision control_time_step,ion_amass
      double precision control_rti
      external         control_frac_coord
      external         control_time_step,ion_amass
      external         control_rti

      call Parallel_taskid(taskid)

*     *****************************************
*     **** seperate ions and charges       ****
*     **** from molecular pseudopotentials ****
*     *****************************************
      call seperate_molpsp(rtdb)
      call seperate_pointcharge(rtdb)


*     **************************
*     **** read in ion data ****
*     **************************
      value = geom_create(geom,'pspwgeometry')
      value = geom_rtdb_load(rtdb,geom,'pspwgeometry')

*     **** get the number ions ****
      value = geom_ncent(geom,nion)

*     *** set ke_count, ke_total, and kg_total to zero ***
      ke_count = 0
      ke_total = 0.0d0
      kg_total = 0.0d0

*     ***** allocate ion data structure *****
      value = MA_alloc_get(mt_dbl,(3*nion),'r2',r2(2),r2(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nion),'r0',r0(2),r0(1))
      value = value.and.
     >       MA_alloc_get(mt_dbl,(nion),'dti',dti(2),dti(1))
      if (.not. value) call errquit('out of heap memory',0)

      call dcopy(3*nion,0.0d0,0,dbl_mb(r2(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r1(1)),1)
      call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)
      do i=1,nion
         value = geom_cent_get(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
         if (.not. value) call errquit('error reading ions',0)
      end do

*     **** define atom list: nkatm,katm,natm,atom ****
      call set_katm_list(1,nion,ion_aname,nkatm,katm,natm,atom)

      if (control_frac_coord()) then
         if (taskid.eq.MASTER) then
           write(*,*) ".... converting from fractional coordinates"
         end if 
        call fcoord_to_real(nion,dbl_mb(r1(1)))
      end if         

*     **** make sure ions are in cell ***
      call incell1(nion,dbl_mb(r1(1)))

      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)


*     **** define dti ****
      dt = control_time_step()
      do i=1,nion
         dbl_mb(dti(1)+i-1) = dt*dt/(ion_amass(i))
      end do

c     value = geom_destroy(geom)
c     call ion_shift()

*     **** get velocities if md code ****
      if (control_code().eq.2) then
*         **** get velocities ****
         value = geom_vel_get(geom,dbl_mb(r0(1)))

*        ***** scale initial velocities and       ****
*        ***** determine initial kinetic energies ****
         call center_v_mass(vgx,vgy,vgz)
         mass=0.0d0
         do i=1,nion
            mass=mass + ion_amass(i)
         end do
         ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)

         eki0 = 0.0d0
         do i=1,nion
            dbl_mb(r0(1)+(i-1)+0) = dbl_mb(r0(1)+(i-1)+0)-vgx
            dbl_mb(r0(1)+(i-1)+1) = dbl_mb(r0(1)+(i-1)+1)-vgy
            dbl_mb(r0(1)+(i-1)+2) = dbl_mb(r0(1)+(i-1)+2)-vgz
            eki0 = eki0 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                                 + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                                 + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
         end do
         eki0=0.5d0*eki0

*        **** scale velocities then find kinetic energy ***
         call dscal((3*nion),control_rti(),dbl_mb(r0(1)),1)
         eki1 = 0.0d0
         do i=1,nion
            eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                                 + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                                 + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
         end do
         eki1=0.5d0*eki1




      end if


      ion_init = value
      return
      end




*     ***************************
*     *							*
*     *	     ion_aname			*
*     *							*
*     ***************************
      character*2 function ion_aname(i)
      implicit none
      integer i

#include "stdio.fh"
#include "geom.fh"

      integer geom
      common / ION2/ geom

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' ion_aname  failed ',i)

      if (.not. geom_tag_to_element(t,symbol,name,q))
     >        call errquit(' ion_aname  failed ',i)
      
      ion_aname = symbol
      return
      end

*     ***************************
*     *							*
*     *	     ion_amass			*
*     *							*
*     ***************************
      real*8 function ion_amass(i)
      implicit none
      integer i

#include "stdio.fh"
#include "geom.fh"


      integer geom
      common / ION2/ geom

*     **** local variables ****
      real*8 mass

      if(.not.geom_mass_get(geom,i,mass))
     >  call errquit(' geom_mass_get  failed ',i)

      ion_amass = mass*1822.89d0
      return
      end

*     ***************************
*     *				*
*     *	     ion_ke		*
*     *				*
*     ***************************
      double precision function ion_ke()
      implicit none

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

      ion_ke=eki1
      return 
      end

*     ***************************
*     *				*
*     *	     ion_com_ke		*
*     *				*
*     ***************************
      double precision function ion_com_ke()
      implicit none

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

      ion_com_ke=ekg
      return 
      end


*     ***************************
*     *				*
*     *	     ion_Temperature    *
*     *				*
*     ***************************
      double precision function ion_Temperature()
      implicit none

***** ion common block ****
#include "ion.fh"

*     **** ion3 common block ****
      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)

      double precision ave,T

      T = 0.0d0
      ave = (ke_total/dble(ke_count))
      if (nion.gt.2) then
         T = 2.0d0*ave/(3.0d0*nion-6.0d0)/kb
      else
         T = 2.0d0*ave/kb
      end if

      ion_Temperature = T
      return 
      end

*     ***************************
*     *				*
*     *	  ion_com_Temperature   *
*     *				*
*     ***************************
      double precision function ion_com_Temperature()
      implicit none

***** ion common block ****
#include "ion.fh"

*     **** ion3 common block ****
      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     *** local variables and parameters ****
      double precision kb
      parameter (kb=3.16679d-6)

      double precision ave,T

      ave = (kg_total/dble(ke_count))
      T   = 2.0d0*ave/kb

      ion_com_Temperature = T
      return 
      end


*     ***************************
*     *				*
*     *	     ion_init_ke	*
*     *				*
*     ***************************
      subroutine ion_init_ke(ekg_out,eki0_out,eki1_out)
      implicit none
      double precision ekg_out
      double precision eki0_out
      double precision eki1_out


      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

      ekg_out  = ekg
      eki0_out = eki0
      eki1_out = eki1

      return
      end

*     ***************************
*     *							*
*     *		ion_nion			*
*     *							*
*     ***************************
      integer function ion_nion()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nion = nion
      return
      end

*     ***************************
*     *							*
*     *		ion_TotalCharge		*
*     *							*
*     ***************************
      real*8 function ion_TotalCharge()
      implicit none

*     **** local variables ****
      integer ia
      real*8  icharge

*     **** external functions ****
      real*8   psp_zv
      integer  ion_natm,ion_nkatm
      external psp_zv
      external ion_natm,ion_nkatm

      icharge = 0.0d0
      do ia=1,ion_nkatm()
         icharge = icharge + ion_natm(ia)*psp_zv(ia)
      end do

      ion_TotalCharge = icharge
      return
      end


*     ***************************
*     *							*
*     *		ion_nkatm			*
*     *							*
*     ***************************
      integer function ion_nkatm()
      implicit none

***** ion common block ****
#include "ion.fh"

      ion_nkatm = nkatm
      return
      end


*     ***************************
*     *							*
*     *		ion_katm			*
*     *							*
*     ***************************

      integer function ion_katm(i)
      implicit none
      integer i

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_katm = int_mb(katm(1)+i-1)
      return
      end

*     ***************************
*     *							*
*     *		ion_natm			*
*     *							*
*     ***************************
      integer function ion_natm(i)
      implicit none
      integer i

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_natm = int_mb(natm(1)+i-1)
      return
      end

*     ***************************
*     *							*
*     *		ion_atom			*
*     *							*
*     ***************************
      character*2 function ion_atom(i)
      implicit none
      integer i

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"
      character*2 aname

      call pspw_copy_atom(1,byte_mb(atom(1)+2*(i-1)),aname)
      ion_atom = aname
      return
      end


*     ***************************
*     *				*
*     *		ion_dti		*
*     *				*
*     ***************************

      double precision function ion_dti(i)
      implicit none
      integer i

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"


      ion_dti = dsqrt(dbl_mb(dti(1)+i-1))
      return
      end



*     ***************************
*     *				*
*     *		ion_rion	*
*     *				*
*     ***************************

      double precision function ion_rion(i,j)
      implicit none
      integer i,j

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"


      ion_rion = dbl_mb(r1(1)+(i-1) + (j-1)*3)
      return
      end

*     ***************************
*     *				*
*     *		ion_vion	*
*     *				*
*     ***************************

      double precision function ion_vion(i,j)
      implicit none
      integer i,j

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      ion_vion = dbl_mb(r0(1)+(i-1) + (j-1)*3)
      return
      end


*     ***************************
*     *				*
*     *	   ion_optimize_step   	*
*     *				*
*     ***************************

      subroutine ion_optimize_step(fion)
      implicit none
      double precision fion(3,*)

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

*     **** local variables ****
      integer i
      double precision scale

      do i=1,nion
         scale = dsqrt(dbl_mb(dti(1)+i-1))
         dbl_mb(r2(1)+(i-1)*3)   = dbl_mb(r1(1)+(i-1)*3) 
     >                           + scale*fion(1,i)
         dbl_mb(r2(1)+(i-1)*3+1) = dbl_mb(r1(1)+(i-1)*3+1) 
     >                           + scale*fion(2,i)
         dbl_mb(r2(1)+(i-1)*3+2) = dbl_mb(r1(1)+(i-1)*3+2) 
     >                           + scale*fion(3,i)
      end do

*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))

*     **** make sure ions are in cell ***
      call incell1(nion,dbl_mb(r2(1)))

      return
      end

*     ***************************
*     *				*
*     *	   ion_verlet_step   	*
*     *				*
*     ***************************

      subroutine ion_verlet_step(fion)
      implicit none
      double precision fion(3,*)

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
      integer i
      double precision scale,dt,h
      double precision vgx,vgy,vgz

*     *** external functions ****
      double precision control_time_step,ion_amass
      external         control_time_step,ion_amass

      dt = control_time_step()
      h = 1.0d0/(2.0d0*dt)

*     **** do a verlet step ***
      do i=1,nion
         scale = (dbl_mb(dti(1)+i-1))
         dbl_mb(r2(1)+(i-1)*3)   = 2.0d0*dbl_mb(r1(1)+(i-1)*3) 
     >                           -       dbl_mb(r0(1)+(i-1)*3) 
     >                           +       scale*fion(1,i)
         dbl_mb(r2(1)+(i-1)*3+1) = 2.0d0*dbl_mb(r1(1)+(i-1)*3+1) 
     >                           -       dbl_mb(r0(1)+(i-1)*3+1) 
     >                           +       scale*fion(2,i)
         dbl_mb(r2(1)+(i-1)*3+2) = 2.0d0*dbl_mb(r1(1)+(i-1)*3+2) 
     >                           -       dbl_mb(r0(1)+(i-1)*3+2) 
     >                           +       scale*fion(3,i)
      end do

*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))

*     **** make sure ions are in cell ***
      call incell3(nion,dbl_mb(r2(1)),
     >                  dbl_mb(r1(1)),
     >                  dbl_mb(r0(1)))

*     **** make r0 the velocity - note that the velocity is deleted ****
*     **** after an ion_shift call                                  ****
      do i=1,nion
         dbl_mb(r0(1)+(i-1)*3)   = h*( dbl_mb(r2(1)+(i-1)*3)
     >                               - dbl_mb(r0(1)+(i-1)*3))
         dbl_mb(r0(1)+(i-1)*3+1) = h*( dbl_mb(r2(1)+(i-1)*3+1)
     >                               - dbl_mb(r0(1)+(i-1)*3+1))
         dbl_mb(r0(1)+(i-1)*3+2) = h*( dbl_mb(r2(1)+(i-1)*3+2)
     >                               - dbl_mb(r0(1)+(i-1)*3+2))
      end do

*     **** determine current kinetic energy and add to running average ****
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1


      return
      end


*     ***************************
*     *							*
*     *	   ion_newton_step   	*
*     *							*
*     ***************************

      subroutine ion_newton_step(fion)
      implicit none
      double precision fion(3,*)

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count
      

*     **** local variables ****
      integer i
      double precision scale,dt
      double precision vgx,vgy,vgz

*     *** external functions ****
      double precision control_time_step,ion_amass
      external         control_time_step,ion_amass

      dt = control_time_step()
      do i=1,nion
         scale = 0.5d0*(dbl_mb(dti(1)+i-1))

         dbl_mb(r2(1)+(i-1)*3)   =    dbl_mb(r1(1)+(i-1)*3) 
     >                           + dt*dbl_mb(r0(1)+(i-1)*3) 
     >                           +       scale*fion(1,i)
         dbl_mb(r2(1)+(i-1)*3+1) =    dbl_mb(r1(1)+(i-1)*3+1) 
     >                           + dt*dbl_mb(r0(1)+(i-1)*3+1) 
     >                           +       scale*fion(2,i)
         dbl_mb(r2(1)+(i-1)*3+2) =    dbl_mb(r1(1)+(i-1)*3+2) 
     >                           + dt*dbl_mb(r0(1)+(i-1)*3+2) 
     >                           +       scale*fion(3,i)
      end do

*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))

*     **** make sure ions are in cell ***
      call incell3(nion,dbl_mb(r2(1)),
     >                  dbl_mb(r1(1)),
     >                  dbl_mb(r0(1)))


*     **** determine current kinetic energy and add to running average ***
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1

      return
      end

*     ***************************
*     *							*
*     *	   ion_nose_step  	 	*
*     *							*
*     ***************************

      subroutine ion_nose_step(ssr,fion)
      implicit none
      double precision ssr
      double precision fion(3,*)

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"

      integer          ke_count
      double precision ekg,eki0,eki1
      double precision ke_total,kg_total,mass
      common / ION3/ ekg,eki0,eki1,mass,
     >               ke_total,kg_total,ke_count

*     **** local variables ****
      integer i
      double precision scale,dt,h
      double precision vgx,vgy,vgz,smr

*     *** external functions ****
      double precision control_time_step,ion_amass
      external         control_time_step,ion_amass

      dt = control_time_step()
      h = 1.0d0/(2.0d0*dt)
      smr = 2.0d0*ssr - 1.0d0

*     **** do a Nose-Hoover verlet step ***
      do i=1,nion
         scale = (dbl_mb(dti(1)+i-1))
         dbl_mb(r2(1)+(i-1)*3)   = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3) 
     >                           -       smr*dbl_mb(r0(1)+(i-1)*3) 
     >                           +       ssr*scale*fion(1,i)
         dbl_mb(r2(1)+(i-1)*3+1) = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3+1) 
     >                           -       smr*dbl_mb(r0(1)+(i-1)*3+1) 
     >                           +       ssr*scale*fion(2,i)
         dbl_mb(r2(1)+(i-1)*3+2) = 2.0d0*ssr*dbl_mb(r1(1)+(i-1)*3+2) 
     >                           -       smr*dbl_mb(r0(1)+(i-1)*3+2) 
     >                           +       ssr*scale*fion(3,i)
      end do

*     ***** impose molecular constraints ***
      call molecular_fix(dbl_mb(r2(1)), dbl_mb(r1(1)))

*     **** make sure ions are in cell ***
      call incell3(nion,dbl_mb(r2(1)),
     >                  dbl_mb(r1(1)),
     >                  dbl_mb(r0(1)))

*     **** make r0 the velocity - note that the velocity is deleted ****
*     **** after an ion_shift call                                  ****
      do i=1,nion
         dbl_mb(r0(1)+(i-1)*3)   = h*( dbl_mb(r2(1)+(i-1)*3)
     >                               - dbl_mb(r0(1)+(i-1)*3))
         dbl_mb(r0(1)+(i-1)*3+1) = h*( dbl_mb(r2(1)+(i-1)*3+1)
     >                               - dbl_mb(r0(1)+(i-1)*3+1))
         dbl_mb(r0(1)+(i-1)*3+2) = h*( dbl_mb(r2(1)+(i-1)*3+2)
     >                               - dbl_mb(r0(1)+(i-1)*3+2))
      end do

*     **** determine current kinetic energy and add to running average ****
      eki1 = 0.0d0
      do i=1,nion
         eki1 = eki1 + ion_amass(i)*( dbl_mb(r0(1)+(i-1)*3 + 0)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 1)**2
     >                              + dbl_mb(r0(1)+(i-1)*3 + 2)**2)
      end do
      eki1=0.5d0*eki1
      call center_v_mass(vgx,vgy,vgz)
      ekg = 0.5d0*mass*(vgx**2+vgy**2+vgz**2)
      ke_total = ke_total + eki1
      kg_total = kg_total + ekg
      ke_count = ke_count + 1


      return
      end



*     ***************************
*     *							*
*     *		ion_shift			*
*     *							*
*     ***************************

      subroutine ion_shift()
      implicit none

#include "mafdecls.fh"

***** ion common block ****
#include "ion.fh"


      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r0(1)),1)
      call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)

      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_e            *
*     *                              *
*     ********************************
*
*   This function computes the (free-space) Coulomb energy between
*   ion cores.
*
*   Entry -
*       nion  ---- number of ions
*       katm  ---- element id number
*       r_ion ---- coordinates of ions
*

      real*8 function ion_ion_e()
      implicit none
      

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,e

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

      e = 0.0d0
      do j=1,ion_nion()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)
          e = e + psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r
        end do
      end do

      ion_ion_e = e
      return
      end


*     ********************************
*     *                              *
*     *         ion_ion_f            *
*     *                              *
*     ********************************
*
*   This routine computes the (free-space) Coulomb forces between
*   ion cores.
*
*   Entry -
*       nion  ---- number of ions
*       katm  ---- element id number
*       r_ion ---- coordinates of ions
*
*   Exit -
*       f_ion ---- force vectors
*

      subroutine ion_ion_f(f_ion)
      implicit none
      real*8 f_ion(3,*)

*     *** local variables ***
      integer i,j
      real*8 x,y,z,r,v

*     *** external functions ***
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

*     *** do nothing if only one atom ***
      if (ion_nion().eq.1) return

      do j=1,ion_nion()
        do i=1,j-1
          x = ion_rion(1,i) - ion_rion(1,j)
          y = ion_rion(2,i) - ion_rion(2,j)
          z = ion_rion(3,i) - ion_rion(3,j)
          r = dsqrt(x*x + y*y + z*z)

          v = psp_zv(ion_katm(i))*psp_zv(ion_katm(j))/r**3

          f_ion(1,i) = f_ion(1,i) + (x*v)
          f_ion(2,i) = f_ion(2,i) + (y*v)
          f_ion(3,i) = f_ion(3,i) + (z*v)

          f_ion(1,j) = f_ion(1,j) - (x*v)
          f_ion(2,j) = f_ion(2,j) - (y*v)
          f_ion(3,j) = f_ion(3,j) - (z*v)
        end do
      end do

      return
      end




*     ********************************
*     *                              *
*     *         ion_Print_XYZ         *
*     *                              *
*     ********************************
*
*   This routine Prints out output in xyz-format
*
*   Entry -
*
*   Exit -
*
      subroutine ion_Print_XYZ(unit)
      implicit none
      integer unit

*     **** local variables ****
      integer ii,taskid,MASTER
      parameter (MASTER=0)


*     **** external functions ****
      character*2 ion_aname
      integer  ion_nion,ion_katm
      real*8   ion_rion
      external ion_aname
      external ion_nion,ion_katm
      external ion_rion


      call Parallel_taskid(taskid)

c     **** ouput xyz format ****
      if (taskid.eq.MASTER) then
        write(unit,1200) 
        write(unit,*) 
        write(unit,*) ion_nion()
        write(unit,*) 
        do ii=1,ion_nion()
           WRITE(unit,1205) ion_aname(ii),
     >                   ion_rion(1,ii)*0.529177d0,
     >                   ion_rion(2,ii)*0.529177d0,
     >                   ion_rion(3,ii)*0.529177d0
        end do
      end if
      return

 1200   FORMAT(//'== XYZ OUTPUT =='/)
 1205   FORMAT(A2,6X,3F14.6)
      end


*     ********************************
*     *                              *
*     *         ion_Print_CIF        *
*     *                              *
*     ********************************
*
*   This routine Prints out output in 
*    Crystallographic Information Files (*.cif)
*
*   Entry -
*
*   Exit -
*
      subroutine ion_Print_CIF(unit)
      implicit none
      integer unit

*     **** local variables ****
      integer i,j,ii,taskid,MASTER
      parameter (MASTER=0)

      real*8 frac(3),a(3,3),b(3,3),volume
      real*8 aa,bb,cc,alpha,beta,gmma,d2,pi
      character*26 dd


*     **** external functions ****
      character*2 ion_aname
      integer  ion_nion,ion_katm
      real*8   ion_rion,lattice_unita
      external ion_aname
      external ion_nion,ion_katm
      external ion_rion,lattice_unita

      call Parallel_taskid(taskid)


*     ***** Determine the unit lattice vectors and distances ******
      do j=1,3
      do i=1,3
        a(i,j) = lattice_unita(i,j)
      end do
      end do

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      volume = a(1,1)*b(1,1)
     >       + a(2,1)*b(2,1)
     >       + a(3,1)*b(3,1)

      volume = 1.0d0/volume
      call dscal(9,volume,b,1)

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      aa = dsqrt(a(1,1)**2 + a(2,1)**2 +a(3,1)**2)
      bb = dsqrt(a(1,2)**2 + a(2,2)**2 +a(3,2)**2)
      cc = dsqrt(a(1,3)**2 + a(2,3)**2 +a(3,3)**2)

      d2 = (a(1,2)-a(1,3))**2 + (a(2,2)-a(2,3))**2 + (a(3,2)-a(3,3))**2
      alpha = (bb*bb + cc*cc - d2)/(2.0d0*bb*cc)
      alpha = dacos(alpha)*180.0d0/pi

      d2 = (a(1,3)-a(1,1))**2 + (a(2,3)-a(2,1))**2 + (a(3,3)-a(3,1))**2
      beta = (cc*cc + aa*aa - d2)/(2.0d0*cc*aa)
      beta = dacos(beta)*180.0d0/pi

      d2 = (a(1,1)-a(1,2))**2 + (a(2,1)-a(2,2))**2 + (a(3,1)-a(3,2))**2
      gmma = (aa*aa + bb*bb - d2)/(2.0d0*aa*bb)
      gmma = dacos(gmma)*180.0d0/pi


      if (taskid.eq.MASTER) then
 
        call util_date(dd)
        write(unit,1200) 
        write(unit,1210) dd(1:24)
        write(unit,1211) 

        write(unit,1220) aa * 0.529177d0
        write(unit,1221) bb * 0.529177d0
        write(unit,1222) cc * 0.529177d0
        write(unit,1223) alpha
        write(unit,1224) beta
        write(unit,1225) gmma

        write(unit,1230)

        write(unit,1240) 
        write(unit,1241) 
        write(unit,1242) 
        write(unit,1243) 
        write(unit,1244) 
        write(unit,1245) 

        do ii=1,ion_nion()
           frac(1) = b(1,1)*ion_rion(1,ii)
     >             + b(2,1)*ion_rion(2,ii)
     >             + b(3,1)*ion_rion(3,ii)
           frac(2) = b(1,2)*ion_rion(1,ii)
     >             + b(2,2)*ion_rion(2,ii)
     >             + b(3,2)*ion_rion(3,ii)
           frac(3) = b(1,3)*ion_rion(1,ii)
     >             + b(2,3)*ion_rion(2,ii)
     >             + b(3,3)*ion_rion(3,ii)
           WRITE(unit,1250) ion_aname(ii),frac

        end do
      end if
      return

 1200 FORMAT('data_nwchem_pspw')
 1210 FORMAT(/'_audit_creation_date   ',A)
 1211 FORMAT(
     > '_audit_creation_method    generated by PSPW module of NWChem')

 1220 FORMAT(//'_cell_length_a   ', F16.4)
 1221 FORMAT(  '_cell_length_b   ', F16.4)
 1222 FORMAT(  '_cell_length_c   ', F16.4)
 1223 FORMAT(  '_cell_angle_alpha', F16.4)
 1224 FORMAT(  '_cell_angle_beta ', F16.4)
 1225 FORMAT(  '_cell_angle_gamma', F16.4)

 1230 FORMAT(/'_symmetry_space_group_name_H-M     P 1  ')

 1240 FORMAT(/'loop_')
 1241 FORMAT('_atom_site_type_symbol')
 1242 FORMAT('_atom_site_label')
 1243 FORMAT('_atom_site_fract_x')
 1244 FORMAT('_atom_site_fract_y')
 1245 FORMAT('_atom_site_fract_z')
 
 1250 FORMAT(A2,6X,3F14.6)
      end
