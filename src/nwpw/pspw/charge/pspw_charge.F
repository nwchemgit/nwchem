c
c $Id: pspw_charge.F,v 1.9 2002-02-27 22:00:46 bylaska Exp $
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

*       **** set charge_type ****
        if (.not.rtdb_get(rtdb,
     >         'nwpw:charge_type',mt_int,1,charge_type))
     >    charge_type = 1
 
        if (charge_type.eq.2) then
           if (.not.rtdb_get(rtdb,
     >         'nwpw:charge_dielectric_d',mt_dbl,1,s_d))
     >       s_d = 0.5d0/0.529177d0

           if (.not.rtdb_get(rtdb,
     >         'nwpw:charge_dielectric_rho',mt_dbl,1,s_rho))
     >       s_rho = 2.5d0/0.529177d0

           if (.not.rtdb_get(rtdb,
     >         'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >       s_sigma = 0.8d0/0.529177d0

        end if
      

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
        if (charge_type.eq.1) then
          write(unit,1160)   
        end if
        if (charge_type.eq.2) then
          write(unit,1170)   
          write(unit,1171) s_d,s_rho,s_sigma   
        end if
        write(unit,1180)   
        do ii=1,nion
          write(unit,1190) ii,pspw_charge_aname(ii),
     >                     (dbl_mb(r1(1)+3*(ii-1)+k-1),k=1,3),
     >                       pspw_charge_charge(ii)

        end do
      end if

      end if
      return

 1160 format(/' Point charges used:')
 1170 format(/' Switching charges used:')
 1171 FORMAT(5X, '(d,rho,sigma) = (',3E11.3,' )  ')
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
      real*8 x1,y1,z1,q1,r,epsilon

*     **** external functions ****
      real*8   util_erf
      external util_erf


      if (nion.gt.0) then

*        **** point charges ****
         if (charge_type.eq.1) then
           do ii=1,nion
              x1 = dbl_mb(r1(1)+3*(ii-1))
              y1 = dbl_mb(r1(1)+3*(ii-1)+1)
              z1 = dbl_mb(r1(1)+3*(ii-1)+2)
              q1 = dbl_mb(charge(1)+ii-1)
              do k=1,n2ft3d
                r = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
                r = dsqrt(r)
                Vqm(k) = Vqm(k) - q1/r
              end do
           end do   

*        **** switching charge potential **** 
         else if (charge_type.eq.2) then
           do ii=1,nion
              x1 = dbl_mb(r1(1)+3*(ii-1))
              y1 = dbl_mb(r1(1)+3*(ii-1)+1)
              z1 = dbl_mb(r1(1)+3*(ii-1)+2)
              q1 = dbl_mb(charge(1)+ii-1)
              do k=1,n2ft3d
                r = (rgrid(1,k)-x1)**2
     >            + (rgrid(2,k)-y1)**2
     >            + (rgrid(3,k)-z1)**2
                r = dsqrt(r)

*               **** define dielectric switching function ****
                if (r.le.s_d) then
                   epsilon = 0.0d0
                else if (r.lt.(s_d+s_rho)) then
                   epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2        
                else
                   epsilon = 1.0d0
                end if

                Vqm(k) = Vqm(k) - epsilon*q1*util_erf(r/s_sigma)/r
              end do
           end do   
         end if

      end if

      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_rho_Fcharge        *
*     *                                *
*     **********************************

      subroutine pspw_charge_rho_Fcharge(n2ft3d,rgrid,rho,dv,fcharge)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fcharge(3,*)

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer i,j
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,verf,yerf,v
      real*8 x,y,z,q,r,epsilon,depsilon,sqrt_pi

*     **** external functions ****
      real*8   util_erf
      external util_erf

      if (nion.gt.0) then
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      do j=1,nion
         x = dbl_mb(r1(1)+3*(j-1))
         y = dbl_mb(r1(1)+3*(j-1)+1)
         z = dbl_mb(r1(1)+3*(j-1)+2)
         q = -dbl_mb(charge(1)+j-1)
         c = 1.0d0/s_sigma
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt( rx**2 + ry**2 + rz**2)

            if (r .gt. 1.0d-8) then

*             **** define dielectric switching function ****
              if (r.le.s_d) then
                 epsilon  = 0.0d0
                 depsilon = 0.0d0
              else if (r.lt.(s_d+s_rho)) then
                 epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2
                 depsilon = 4.0d0*((r-s_d)/s_rho**2)
     >                           *(1.0d0-(r-s_d)**2/s_rho**2)        
              else
                 epsilon  = 1.0d0
                 depsilon = 0.0d0
              end if
              yerf=r*c
              verf = util_erf(yerf)
              v    = epsilon*q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                          - verf)/r**3
     >             + depsilon*q*verf/r**2
            else
              v = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv
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




*     **********************************
*     *	                               *
*     *     pspw_charge_Fion_Fcharge   *
*     *                                *
*     **********************************

*    This routine calculates the forces between the QM ions and point charges.
*
*     Exit - fion: force on the QM ions
*            fcharge: force on the point charges
*
*     Uses - ion_nion,Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_Fion_Fcharge(fion,fcharge)
      implicit none
      real*8 fion(3,*)
      real*8 fcharge(3,*)

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,j,nion_ion
      real*8  qi,ri(3),qj,rj(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      if (nion.gt.0) then      
     
        nion_ion = ion_nion()

        do i=1,nion_ion
           qi = psp_zv(ion_katm(i))
           ri(1) = ion_rion(1,i)
           ri(2) = ion_rion(2,i)
           ri(3) = ion_rion(3,i)
           do j=1,nion
              qj = dbl_mb(charge(1)+j-1)            
              rj(1) = dbl_mb(r1(1)+3*(j-1))
              rj(2) = dbl_mb(r1(1)+3*(j-1)+1)
              rj(3) = dbl_mb(r1(1)+3*(j-1)+2)
              call Q_Electrostatic_Force(ri,qi,fion(1,i),
     >                                 rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_Energy_charge  *
*     *                                *
*     **********************************

      real*8 function pspw_charge_Energy_charge()
      implicit none

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,ii
      real*8  qi, xi, yi, zi
      real*8  qii,xii,yii,zii
      real*8  r,energy

*     **** external functions ****
      logical  nwpw_bqbq
      external nwpw_bqbq

      energy = 0.0d0
      if ((nion.gt.0).and.nwpw_bqbq()) then
        do ii=1,nion
           qii = dbl_mb(charge(1)+ii-1)            
           xii = dbl_mb(r1(1)+3*(ii-1))
           yii = dbl_mb(r1(1)+3*(ii-1)+1)
           zii = dbl_mb(r1(1)+3*(ii-1)+2)
           qii = dbl_mb(charge(1)+ii-1)
           do i=1,ii-1
              xi = dbl_mb(r1(1)+3*(i-1))
              yi = dbl_mb(r1(1)+3*(i-1)+1)
              zi = dbl_mb(r1(1)+3*(i-1)+2)
              qi = dbl_mb(charge(1)+i-1)
              r = (xi-xii)**2 + (yi-yii)**2 + (zi-zii)**2
              r = dsqrt(r)
              energy = energy + qi*qii/r
           end do
        end do
      end if

      pspw_charge_Energy_charge = energy
      return
      end


*     **********************************
*     *	                               *
*     *     pspw_charge_Fcharge        *
*     *                                *
*     **********************************

*    This routine calculates the force between the point charges.
*
*     Exit - fcharge: force on the point charges
*
*     Uses - Q_Electrostatic_Force
*
*     Author - Eric Bylaska

      subroutine pspw_charge_Fcharge(fcharge)
      implicit none
      real*8 fcharge(3,*)

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     **** local variables ****
      integer i,j
      real*8  qi,ri(3),qj,rj(3)

*     **** external functions ****
      logical  nwpw_bqbq
      external nwpw_bqbq

      if ((nion.gt.0).and.(nwpw_bqbq())) then      
     
        do j=1,nion
           qj = dbl_mb(charge(1)+j-1)            
           rj(1) = dbl_mb(r1(1)+3*(j-1))
           rj(2) = dbl_mb(r1(1)+3*(j-1)+1)
           rj(3) = dbl_mb(r1(1)+3*(j-1)+2)
           do i=1,j-1
              qi = dbl_mb(charge(1)+i-1)            
              ri(1) = dbl_mb(r1(1)+3*(i-1))
              ri(2) = dbl_mb(r1(1)+3*(i-1)+1)
              ri(3) = dbl_mb(r1(1)+3*(i-1)+2)
              call Q_Electrostatic_Force(ri,qi,fcharge(1,i),
     >                                   rj,qj,fcharge(1,j))
           end do
        end do

      end if
      return
      end








