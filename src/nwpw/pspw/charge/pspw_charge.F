c
c $Id: pspw_charge.F,v 1.5 2002-01-17 01:24:51 bylaska Exp $
c

***********************************************************
*                                                         *
*                   PSPW-Charge module                    *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*    -- developed by Eric J. Bylaska on October 18,2001   *
*                                                         *
***********************************************************
*
*
*  

*     **********************************
*     *	                               *
*     *          pspw_charge_init      *
*     *                                *
*     **********************************
 
      subroutine pspw_charge_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "pspw_charge.fh"


*     **** local variables ****
      logical value
      integer i
      character*16 t

*     **** external functions ****
      logical     ion_chargeexist
      external    ion_chargeexist


*     *************************
*     **** read in charge data ****
*     *************************
      if (ion_chargeexist()) then
    
         value = geom_create(geom,'chargegeometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'chargegeometry')
         value = value.and.geom_ncent(geom,nion)
         if (.not. value) 
     >     call errquit('pspw_charge_init: opening chargegeometry',0)
      else
         nion = 0
      end if


*     **** don't do anything if there are no charges atoms ****
      if (nion.gt.0) then

*        ***** Allocate pspw_charge geometries, and charges *****
         value = MA_alloc_get(mt_dbl,(3*nion),
     >                        'charge_r1',r1(2),r1(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,(nion),
     >                        'charge_charge',charge(2),charge(1))
         if (.not. value) 
     >      call errquit('pspw_charge_init: pushing heap memory',0)


*        **** read in charge geometries, and charges ****
         do i=1,nion
             value = value.and.
     >       geom_cent_get(geom,i,t,
     >                     dbl_mb(r1(1)+(i-1)*3),
     >                     dbl_mb(charge(1)+i-1))
         end do
         if (.not. value) 
     >   call errquit('pspw_charge_init:error reading chargegeometry',0)

        charge_type = 1
      end if

      return
      end



*     **********************************
*     *	                               *
*     *        pspw_charge_destroy       *
*     *                                *
*     **********************************
      subroutine pspw_charge_destroy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      logical value

      if (nion.gt.0) then
        value = geom_destroy(geom)
        if (.not. value) 
     >     call errquit('error destroying chargegeometry', 0)
      end if

      return
      end

*     **********************************
*     *	                               *
*     *        pspw_charge_write       *
*     *                                *
*     **********************************

      subroutine pspw_charge_write(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      logical value
      integer i
      double precision rxyz(3),q
      character*16 t

      integer  control_code
      external control_code

*     *******************************
*     **** write out charge data ****
*     *******************************
      if (nion.gt.0) then
        value = .true.
        do i=1,nion
           value = value.and.geom_cent_get(geom,i,t,rxyz,q)
           value = value.and.
     >             geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
        end do

        value = value.and.geom_rtdb_delete(rtdb,'chargegeometry')
        if ((control_code().eq.1).or.
     >      (control_code().eq.2)) then
        value = value.and.geom_rtdb_store(rtdb,geom,'chargegeometry')
        end if
        value = value.and.geom_destroy(geom)
        if (.not. value) call errquit('error writing chargegeometry', 0)
      end if

      return
      end




*     **********************************
*     *	                               *
*     *         pspw_charge_aname        *
*     *                                *
*     **********************************
      character*2 function pspw_charge_aname(i)
      implicit none
      integer i

#include "stdio.fh"
#include "geom.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' pspw_charge_aname: failed ',i)

      if (.not. geom_tag_to_element(t,symbol,name,q)) then
        symbol = 'bq'
      end if

      pspw_charge_aname = symbol
      return
      end


*     ***************************
*     *                         *
*     *   pspw_charge_charge    *
*     *                         *
*     ***************************
      real*8 function pspw_charge_charge(i)
      implicit none
      integer i

#include "mafdecls.fh"
#include "pspw_charge.fh"


      pspw_charge_charge = dbl_mb(charge(1)+i-1)
      return
      end


*     **********************************
*     *	                               *
*     *         pspw_charge_found        *
*     *                                *
*     **********************************
      logical function pspw_charge_found()
      implicit none

#include "pspw_charge.fh"
 
      logical value

      value = .false.
      if (nion.gt.0) value = .true.

      pspw_charge_found = value
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_charge_end	       *
*     *                                *
*     **********************************

      subroutine pspw_charge_end()
      implicit none

#include "mafdecls.fh"
#include "pspw_charge.fh"

      logical value

      if (nion.gt.0) then
        value = MA_free_heap(r1(2))
        value = value.and.MA_free_heap(charge(2))
        if (.not.value) 
     >    call errquit('pspw_charge_end: freeing heap',0)
      end if
      return
      end

*     **********************************
*     *	                               *
*     *        pspw_charge_nion        *
*     *                                *
*     **********************************
      integer function pspw_charge_nion()
      implicit none

#include "pspw_charge.fh"

      pspw_charge_nion = nion
      return
      end


*     **********************************
*     *	                               *
*     *       pspw_charge_Print        *
*     *                                *
*     **********************************

      subroutine pspw_charge_Print(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_charge_aname
      real*8      pspw_charge_charge
      external    pspw_charge_aname
      external    pspw_charge_charge

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        write(unit,1180)   
        do ii=1,nion
          write(unit,1190) ii,pspw_charge_aname(ii),
     >                     (dbl_mb(r1(1)+3*(ii-1)+k-1),k=1,3),
     >                       pspw_charge_charge(ii)

        end do
      end if

      end if
      return

 1180 FORMAT(/' position of point charges (au): ')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - charge = ',F10.3)
      end




*     **********************************
*     *	                               *
*     *       pspw_charge_PrintXYZ       *
*     *                                *
*     **********************************

      subroutine pspw_charge_PrintXYZ(unit)
      implicit none
      integer unit

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_charge_aname
      external    pspw_charge_aname

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then  
       do ii=1,nion
        write(unit,*) pspw_charge_aname(ii),'      ',
     >               (dbl_mb(r1(1)+3*(ii-1)+k-1)*0.529177d0,k=1,3)
       end do
      end if

      end if

      return
      end


*     **********************************
*     *	                               *
*     *        pspw_charge_rion        *
*     *                                *
*     **********************************
      real*8 function pspw_charge_rion(i,ii)
      implicit none
      integer i,ii

#include "mafdecls.fh"
#include "pspw_charge.fh"

      pspw_charge_rion = dbl_mb(r1(1)+3*(ii-1)+i-1)
      return
      end



*     **********************************
*     *	                               *
*     *     pspw_charge_Generate_V     *
*     *                                *
*     **********************************

      subroutine pspw_charge_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer ii,k
      real*8 x1,y1,z1,q1,r

      if (nion.gt.0) then

         do ii=1,nion
            x1 = dbl_mb(r1(1)+3*(ii-1))
            y1 = dbl_mb(r1(1)+3*(ii-1)+1)
            z1 = dbl_mb(r1(1)+3*(ii-1)+2)
            q1 = dbl_mb(charge(1)+ii-1)
            do k=1,n2ft3d
              r = (rgrid(1,k)-x1)**2
     >          + (rgrid(2,k)-y1)**2
     >          + (rgrid(3,k)-z1)**2
              r = dsqrt(r)
              Vqm(k) = Vqm(k) - q1/r
            end do
         end do   

      end if

      return
      end

*     **********************************
*     *	                               *
*     *     pspw_charge_Energy_ion     *
*     *                                *
*     **********************************

      real*8 function pspw_charge_Energy_ion()
      implicit none

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,ii,ia
      real*8  qi, xi, yi, zi
      real*8  qii,xii,yii,zii
      real*8  r,energy

*     **** external functions ****
      integer   ion_nion,ion_katm
      real*8    ion_rion,psp_zv
      external  ion_nion,ion_katm
      external  ion_rion,psp_zv

      energy = 0.0d0
      do ii=1,ion_nion()
         ia = ion_katm(ii)
         qii = psp_zv(ia)
         xii = ion_rion(1,ii)
         yii = ion_rion(2,ii)
         zii = ion_rion(3,ii)
         do i=1,nion
            xi = dbl_mb(r1(1)+3*(i-1))
            yi = dbl_mb(r1(1)+3*(i-1)+1)
            zi = dbl_mb(r1(1)+3*(i-1)+2)
            qi = dbl_mb(charge(1)+i-1)
            r = (xi-xii)**2 + (yi-yii)**2 + (zi-zii)**2
            r = dsqrt(r)
            energy = energy + qi*qii/r
         end do
      end do

      pspw_charge_Energy_ion = energy
      return
      end
