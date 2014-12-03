********************************************************
*                                                        *
*              Water pseudopotential module              *
*                                                        *
*          BLCJ Electron-Water pseudopotential           *
*                                                        *
*          Interfaced to nwchem-PSPW code                *
*                                                        *
*    -- developed by Eric J. Bylaska on  July 19,2001    *
*                                                        *
**********************************************************


*     **********************************
*     *	                               *
*     *          Waterpsp_init         *
*     *                                *
*     **********************************

      subroutine Waterpsp_init(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"


*     ***** water common block ****
#include "waterpsp.fh"

      logical mmexist
      common / pspwqmm/ mmexist

*     **** local variables ****
      logical value
      integer i,geom,nion
      double precision q
      character*16 t

      integer  control_version
      external control_version


*     ****************************
*     **** read in water data ****
*     ****************************
      if (mmexist) then
         value = geom_create(geom,'qmmmgeometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'qmmmgeometry')
         value = value.and.geom_ncent(geom,nion)
         if (.not. value) call errquit('opening qmmgeometry',0,
     &       GEOM_ERR)
         nwater = nion/3
      else
         nwater = 0
      end if


      if (nwater.gt.0) then

*        ***** Boundary condtion check ****
         if (control_version().ne.4) then
          value = geom_destroy(geom)
          call errquit(
     >    'Aperiodic boundary conditions must be used with Waterpsp',0,
     &       GEOM_ERR)
         end if

*     ***** allocate waterpsp data structure *****
      value = BA_alloc_get(mt_dbl,(3*nion),
     .                     'rwater2',rwater2(2),rwater2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nion),
     >                     'rwater1',rwater1(2),rwater1(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nion),
     >                     'rwater0',rwater0(2),rwater0(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      do i=1,nion
          value = value.and.
     >            geom_cent_get(geom,i,t,dbl_mb(rwater1(1)+(i-1)*3),q)
      end do
      call dcopy(3*nion,dbl_mb(rwater1(1)),1,dbl_mb(rwater2(1)),1)

      call dcopy(3*nion,0.0d0,0,dbl_mb(rwater0(1)),1)
      value = value.and.geom_vel_get(geom,dbl_mb(rwater0(1)))
      value = value.and.geom_destroy(geom)
      if (.not. value) call errquit('error reading qmmmgeometry', 0,
     &       GEOM_ERR)
      call dcopy(3*nion,0.0d0,0,dbl_mb(rwater0(1)),1)

*     *** Init BLCJ SR Potential and LJparam  ****
      call LJparam_init(rtdb)
      call BLCJ_SR_init(rtdb)
      call BLCJ_LR_init(rtdb)
      end if

*     **** set ekw_total,ekw_count ****
      ekw_total = 0.0d0
      ekw_count = 0

      return
      end


*     **********************************
*     *	                               *
*     *        Waterpsp_write          *
*     *                                *
*     **********************************

      subroutine Waterpsp_write(rtdb)
      implicit none
      integer rtdb

***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer i,geom,nion
      double precision rxyz(3),q
      character*16 t

      if (nwater.gt.0) then
      value = geom_create(geom,'qmmmgeometry')
      value = geom_rtdb_load(rtdb,geom,'qmmmgeometry')
      value = geom_ncent(geom,nion)
      do i=1,nion
          value = geom_cent_get(geom,i,t,rxyz,q)
          value = geom_cent_set(geom,i,t,dbl_mb(rwater1(1)+(i-1)*3),q)
      end do
      value = geom_vel_set(geom,dbl_mb(rwater0(1)))

      value = geom_rtdb_delete(rtdb,'qmmmgeometry')
      value = geom_rtdb_store(rtdb,geom,'qmmmgeometry')
      value = geom_destroy(geom)
      end if

      return
      end

*     **********************************
*     *	                               *
*     *           Waterpsp_found       *
*     *                                *
*     **********************************
      logical function Waterpsp_found()
      implicit none

***** water common block ****
#include "waterpsp.fh"

      logical value
 
      value = .false.
      if (nwater.gt.0) value = .true.

      Waterpsp_found = value
      return
      end

*     **********************************
*     *	                               *
*     *           Waterpsp_end 	       *
*     *                                *
*     **********************************

      subroutine Waterpsp_end()
      implicit none

***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

      logical value

      if (nwater.gt.0) then
        value = BA_free_heap(rwater2(2))
        value = BA_free_heap(rwater1(2))
        value = BA_free_heap(rwater0(2))
        call LJparam_end()
        call BLCJ_SR_end()
      end if

      return
      end

*     **********************************
*     *	                               *
*     *        Waterpsp_nwater         *
*     *                                *
*     **********************************
      integer function Waterpsp_nwater()
      implicit none

*     ***** water common block ****
#include "waterpsp.fh"

      Waterpsp_nwater = nwater
      return
      end

*     **********************************
*     *	                               *
*     *        Waterpsp_ke	           *
*     *                                *
*     **********************************
      real*8 function Waterpsp_ke()
      implicit none

*     ***** water common block ****
#include "waterpsp.fh"

      Waterpsp_ke = ekw
      return
      end

*     **********************************
*     *	                               *
*     *        Waterpsp_Temperature    *
*     *                                *
*     **********************************
      real*8 function Waterpsp_Temperature()
      implicit none

*     ***** water common block ****
#include "waterpsp.fh"

      real*8 kb
      parameter (kb=3.16679d-6)

      Waterpsp_Temperature = 2.0d0*ekw_total/dble(ekw_count)
     >                                      /(9.0d0*nwater-6.0d0)
     >                                      /kb
      return
      end


*     **********************************
*     *	                               *
*     *        Waterpsp_Print          *
*     *                                *
*     **********************************

      subroutine Waterpsp_Print(unit)
      implicit none
      integer unit

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k,i1,i2,i3

      if (nwater.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        write(unit,1180)   
        do ii=1,nwater
          i1 = 9*(ii-1)
          i2 = i1+3
          i3 = i2+3
          write(unit,1190) ii,'O^',(dbl_mb(rwater1(1)+k-1+i1),k=1,3)
          write(unit,1190) ii,'H^',(dbl_mb(rwater1(1)+k-1+i2),k=1,3)
          write(unit,1190) ii,'H^',(dbl_mb(rwater1(1)+k-1+i3),k=1,3)    
        end do
      end if

      end if
      return

 1180 FORMAT(/' position of pseudopotential water:')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
      end

*     **********************************
*     *	                               *
*     *        Waterpsp_Print2         *
*     *                                *
*     **********************************

      subroutine Waterpsp_Print2(unit)
      implicit none
      integer unit

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k,i1,i2,i3

      if (nwater.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        write(unit,1180)   
        do ii=1,nwater
          i1 = 9*(ii-1)
          i2 = i1+3
          i3 = i2+3
          write(unit,1190) ii,'O^',(dbl_mb(rwater0(1)+k-1+i1),k=1,3)
          write(unit,1190) ii,'H^',(dbl_mb(rwater0(1)+k-1+i2),k=1,3)
          write(unit,1190) ii,'H^',(dbl_mb(rwater0(1)+k-1+i3),k=1,3)    
        end do
      end if

      end if
      return

 1180 FORMAT(/' velocity of pseudopotential water:')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
      end


*     **********************************
*     *	                               *
*     *         Waterpsp_PrintXYZ      *
*     *                                *
*     **********************************

      subroutine Waterpsp_PrintXYZ(unit)
      implicit none
      integer unit

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k,i1,i2,i3

      if (nwater.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then  
       do ii=1,nwater
          i1 = 9*(ii-1)
          i2 = i1+3
          i3 = i2+3
        write(unit,*) 'O ',(dbl_mb(rwater1(1)+k-1+i1)*0.529177d0,k=1,3)
        write(unit,*) 'H ',(dbl_mb(rwater1(1)+k-1+i2)*0.529177d0,k=1,3)
        write(unit,*) 'H ',(dbl_mb(rwater1(1)+k-1+i3)*0.529177d0,k=1,3)    
       end do
      end if

      end if

      return
      end

*     **********************************
*     *	                               *
*     *         Waterpsp_rwater        *
*     *                                *
*     **********************************
      real*8 function Waterpsp_rwater(i,j,ii)
      implicit none
      integer i,j,ii

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

      integer index
      index = 9*(ii-1) + 3*(j-1) + (i-1)
      Waterpsp_rwater = dbl_mb(rwater1(1)+index)
      return
      end

*     **********************************
*     *	                               *
*     *         Waterpsp_vwater        *
*     *                                *
*     **********************************
      real*8 function Waterpsp_vwater(i,j,ii)
      implicit none
      integer i,j,ii

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

      integer index
      index = 9*(ii-1) + 3*(j-1) + (i-1)
      Waterpsp_vwater = dbl_mb(rwater0(1)+index)
      return
      end


*     **********************************
*     *	                               *
*     *         Waterpsp_shift         *
*     *                                *
*     **********************************
      subroutine Waterpsp_shift()
      implicit none

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

      call dcopy((9*nwater),dbl_mb(rwater1(1)),1,dbl_mb(rwater0(1)),1)
      call dcopy((9*nwater),dbl_mb(rwater2(1)),1,dbl_mb(rwater1(1)),1)

      return
      end

*     **********************************
*     *	                               *
*     *     Waterpsp_Generate_V        *
*     *                                *
*     **********************************

      subroutine Waterpsp_Generate_V(n2ft3d,rgrid,Vwpsp)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vwpsp(*)

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer ii,i1,i2,i3


      if (nwater.gt.0) then
      do ii=1,nwater
         i1 = 9*(ii-1)
         i2 = i1+3
         i3 = i2+3
         call BLCJ_LR(dbl_mb(rwater1(1)+i1),
     >                dbl_mb(rwater1(1)+i2),
     >                dbl_mb(rwater1(1)+i3),
     >                n2ft3d,rgrid,
     >                Vwpsp)
         call BLCJ_SR(dbl_mb(rwater1(1)+i1),
     >                dbl_mb(rwater1(1)+i2),
     >                dbl_mb(rwater1(1)+i3),
     >                  n2ft3d,rgrid, 
     >                  Vwpsp)
      end do
      end if

      return
      end



*     **********************************
*     *	                               *
*     *      Waterpsp_Energy_Water     *
*     *                                *
*     **********************************

      subroutine Waterpsp_Energy_Water(ewater)
      implicit none
      real*8 ewater

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer ii,jj
      integer i1,i2,i3
      integer j1,j2,j3

*     **** external functions ****
      real*8   BLCJ_Intramolecular,BLCJ_Intermolecular
      external BLCJ_Intramolecular,BLCJ_Intermolecular


      ewater = 0.0d0
      if (nwater.gt.0) then
*     **** calculate Intramolecular contributions ****
      do ii=1,nwater
         i1 = 9*(ii-1)
         i2 = i1+3
         i3 = i2+3
         ewater = ewater 
     >          + BLCJ_Intramolecular(dbl_mb(rwater1(1)+i1),
     >                                dbl_mb(rwater1(1)+i2),
     >                                dbl_mb(rwater1(1)+i3))
   
      end do

*     **** calculate Intermolecular contributions ****
      do ii=1,nwater-1
       i1 = 9*(ii-1)
       i2 = i1+3
       i3 = i2+3
        do jj=ii+1,nwater
         j1 = 9*(jj-1)
         j2 = j1+3
         j3 = j2+3
         ewater = ewater 
     >          + BLCJ_Intermolecular(dbl_mb(rwater1(1)+i1),
     >                                dbl_mb(rwater1(1)+i2),
     >                                dbl_mb(rwater1(1)+i3),
     >                                dbl_mb(rwater1(1)+j1),
     >                                dbl_mb(rwater1(1)+j2),
     >                                dbl_mb(rwater1(1)+j3))
        end do
      end do

      end if

      return
      end

*     **********************************
*     *	                               *
*     *        Waterpsp_Energy_ion     *
*     *                                *
*     **********************************

      subroutine Waterpsp_Energy_ion(ewater)
      implicit none

      real*8 ewater

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer ii,j,i1,i2,i3
      real*8 zv,ec,er,ep
      real*8 ei,si,ew,sw
      real*8 rion(3)

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   psp_zv,ion_rion
      external ion_nion,ion_katm
      external psp_zv,ion_rion

      ewater = 0.0d0    

      if (nwater.gt.0) then
      call LJparam_water(ew,sw)
      do j=1,ion_nion()
        zv = psp_zv(ion_katm(j))
        rion(1) = ion_rion(1,j)
        rion(2) = ion_rion(2,j)
        rion(3) = ion_rion(3,j)
        call LJparam_ion(j,ei,si)        
        do ii=1,nwater
           i1 = 9*(ii-1)
           i2 = i1+3
           i3 = i2+3
           ec = 0.0d0
           call BLCJ_ion_Coulomb(dbl_mb(rwater1(1)+i1),
     >                           dbl_mb(rwater1(1)+i2),
     >                           dbl_mb(rwater1(1)+i3),
     >                           rion,zv,
     >                           ec)
           ep = 0.0d0
           call BLCJ_ion_Polarization(dbl_mb(rwater1(1)+i1),
     >                                dbl_mb(rwater1(1)+i2),
     >                                dbl_mb(rwater1(1)+i3),
     >                                rion,zv,
     >                                ep)
           er = 0.0d0
           call LJ_ion_Repulsion(dbl_mb(rwater1(1)+i1),ew,sw,
     >                           rion,ei,si,er)
            
           ewater = ewater + ec + ep + er
        end do
      end do

      end if

      return
      end

*     **********************************
*     *	                               *
*     *         Waterpsp_Fion          *
*     *                                *
*     **********************************
      subroutine Waterpsp_Fion(fion)
      implicit none
      real*8 fion(3,*)
     
*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer ii,j,i1,i2,i3
      real*8 fc(3),fp(3),fl(3)
      real*8 zv,ei,si,ew,sw
      real*8 rion(3)


*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion


*     **** calculate the force on QM ions ****
      call LJparam_water(ew,sw)
      do j=1,ion_nion()
        zv = psp_zv(ion_katm(j))
        rion(1) = ion_rion(1,j)
        rion(2) = ion_rion(2,j)
        rion(3) = ion_rion(3,j)
        call LJparam_ion(j,ei,si)

        do ii=1,nwater 
          i1 = 9*(ii-1)
          i2 = i1+3
          i3 = i2+3

         call dcopy(3,0.0d0,0,fc,1)
         call dcopy(3,0.0d0,0,fp,1)
         call dcopy(3,0.0d0,0,fl,1)
         call BLCJ_ion_Coulomb_Fion(dbl_mb(rwater1(1) + i1),
     >                              dbl_mb(rwater1(1) + i2),
     >                              dbl_mb(rwater1(1) + i3),
     >                              rion,zv,fc)
         call BLCJ_ion_Polarization_Fion(dbl_mb(rwater1(1)+i1),
     >                                   dbl_mb(rwater1(1)+i2),
     >                                   dbl_mb(rwater1(1)+i3),
     >                                   rion,zv,fp)
         call LJ_ion_Repulsion_Force(dbl_mb(rwater1(1)+i1),ew,sw,
     >                              rion,ei,si,fl)
 
         fion(1,j) = fion(1,j) + fc(1) + fp(1) + fl(1)
         fion(2,j) = fion(2,j) + fc(2) + fp(2) + fl(2)
         fion(3,j) = fion(3,j) + fc(3) + fp(3) + fl(3)
        end do

      end do
      return
      end


*     **********************************
*     *	                               *
*     *         Waterpsp_Fwater        *
*     *                                *
*     **********************************
      subroutine Waterpsp_Fwater(n2ft3d,rgrid,rho,dv,fwater)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fwater(3,3,*)
     

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer j,ii,i1,i2,i3
      integer jj,j1,j2,j3
      real*8 ew,sw,ei,si,zv,rion(3)
      real*8 foww(3),f1ww(3),f2ww(3)
      real*8 foiw(3),f1iw(3),f2iw(3)
      real*8 foew(3),f1ew(3),f2ew(3)
      real*8 foc(3),f1c(3),f2c(3),fop(3)
      real*8 fol(3),f1l(3),f2l(3)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   psp_zv,ion_rion
      external ion_katm,ion_nion
      external psp_zv,ion_rion

      call LJparam_water(ew,sw)

*     **** calculation the force on waters ****
      do ii=1,nwater
         i1 = 9*(ii-1)
         i2 = i1+3
         i3 = i2+3

*        **** Intra water-water interaction ****
         call dcopy(3,0.0d0,0,foww,1)
         call dcopy(3,0.0d0,0,f1ww,1)
         call dcopy(3,0.0d0,0,f2ww,1)
         call BLCJ_Intramolecular_Fwater(dbl_mb(rwater1(1)+i1),
     >                                   dbl_mb(rwater1(1)+i2),
     >                                   dbl_mb(rwater1(1)+i3),
     >                                   foww,f1ww,f2ww)
    
*        **** ion-water interaction ****
         call dcopy(3,0.0d0,0,foiw,1)
         call dcopy(3,0.0d0,0,f1iw,1)
         call dcopy(3,0.0d0,0,f2iw,1)
         do j=1,ion_nion()
            zv = psp_zv(ion_katm(j))
            rion(1) = ion_rion(1,j)
            rion(2) = ion_rion(2,j)
            rion(3) = ion_rion(3,j)
            call LJparam_ion(j,ei,si)

            call dcopy(3,0.0d0,0,foc,1)
            call dcopy(3,0.0d0,0,f1c,1)
            call dcopy(3,0.0d0,0,f2c,1)
            call dcopy(3,0.0d0,0,fop,1)
            call dcopy(3,0.0d0,0,fol,1)
            call BLCJ_ion_Coulomb_Fwater(dbl_mb(rwater1(1)+i1),
     >                                   dbl_mb(rwater1(1)+i2), 
     >                                   dbl_mb(rwater1(1)+i3),
     >                                   rion,zv,foc,f1c,f2c)
            call BLCJ_ion_Polarization_Fwater(dbl_mb(rwater1(1)+i1),
     >                                        dbl_mb(rwater1(1)+i2), 
     >                                        dbl_mb(rwater1(1)+i3),
     >                                        rion,zv,fop)
            call LJ_ion_Repulsion_Force(rion,ei,si,
     >                                  dbl_mb(rwater1(1)+i1),ew,sw,
     >                                  fol)

            foiw(1) = foiw(1) + foc(1) + fop(1) + fol(1)
            foiw(2) = foiw(2) + foc(2) + fop(2) + fol(2)
            foiw(3) = foiw(3) + foc(3) + fop(3) + fol(3)

            f1iw(1) = f1iw(1) + f1c(1)
            f1iw(2) = f1iw(2) + f1c(2)
            f1iw(3) = f1iw(3) + f1c(3)

            f2iw(1) = f2iw(1) + f2c(1)
            f2iw(2) = f2iw(2) + f2c(2)
            f2iw(3) = f2iw(3) + f2c(3)

         end do

*        **** elc-water interaction ****
         call dcopy(3,0.0d0,0,foc,1)
         call dcopy(3,0.0d0,0,f1c,1)
         call dcopy(3,0.0d0,0,f2c,1)
         call dcopy(3,0.0d0,0,fol,1)
         call dcopy(3,0.0d0,0,f1l,1)
         call dcopy(3,0.0d0,0,f2l,1)
         call BLCJ_LR_Fwater(dbl_mb(rwater1(1)+i1),
     >                       dbl_mb(rwater1(1)+i2),
     >                       dbl_mb(rwater1(1)+i3),
     >                       n2ft3d,rgrid,rho,dv,
     >                       foc,f1c,f2c)
         call BLCJ_SR_Fwater(dbl_mb(rwater1(1)+i1),
     >                       dbl_mb(rwater1(1)+i2),
     >                       dbl_mb(rwater1(1)+i3),
     >                       n2ft3d,rgrid,rho,dv,
     >                       fol,f1l,f2l)

          foew(1) = foc(1) + fol(1) 
          foew(2) = foc(2) + fol(2) 
          foew(3) = foc(3) + fol(3) 

          f1ew(1) = f1c(1) + f1l(1) 
          f1ew(2) = f1c(2) + f1l(2) 
          f1ew(3) = f1c(3) + f1l(3) 

          f2ew(1) = f2c(1) + f2l(1) 
          f2ew(2) = f2c(2) + f2l(2) 
          f2ew(3) = f2c(3) + f2l(3) 

*        **** sum up contributions ****
         fwater(1,1,ii) = foww(1) + foiw(1) + foew(1)
         fwater(2,1,ii) = foww(2) + foiw(2) + foew(2)
         fwater(3,1,ii) = foww(3) + foiw(3) + foew(3)
    
         fwater(1,2,ii) = f1ww(1) + f1iw(1) + f1ew(1)
         fwater(2,2,ii) = f1ww(2) + f1iw(2) + f1ew(2)
         fwater(3,2,ii) = f1ww(3) + f1iw(3) + f1ew(3)
    
         fwater(1,3,ii) = f2ww(1) + f2iw(1) + f2ew(1)
         fwater(2,3,ii) = f2ww(2) + f2iw(2) + f2ew(2)
         fwater(3,3,ii) = f2ww(3) + f2iw(3) + f2ew(3)

      end do

*     **** Inter water-water interaction ****
      do ii=1,nwater
         i1 = 9*(ii-1)
         i2 = i1+3
         i3 = i2+3
         do jj=ii+1,nwater
            j1 = 9*(jj-1)
            j2 = j1+3
            j3 = j2+3
            call BLCJ_Intermolecular_Fwater(dbl_mb(rwater1(1)+i1),
     >                                      dbl_mb(rwater1(1)+i2),
     >                                      dbl_mb(rwater1(1)+i3),
     >               fwater(1,1,ii),fwater(1,2,ii),fwater(1,3,ii),
     >                                      dbl_mb(rwater1(1)+j1),
     >                                      dbl_mb(rwater1(1)+j2),
     >                                      dbl_mb(rwater1(1)+j3),
     >               fwater(1,1,jj),fwater(1,2,jj),fwater(1,3,jj))
         end do
      end do

      return
      end



*     **********************************
*     *	                               *
*     *         Waterpsp_Update  	   *
*     *                                *
*     **********************************
      subroutine Waterpsp_Update(algorithm,
     >                           n2ft3d,rgrid,rho,dv,
     >                           dt,fion,deltawater)
      implicit none
#include "errquit.fh"
      integer algorithm
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 dt
      real*8 fion(3,*)
      real*8 deltawater

*     ***** water common block ****
#include "waterpsp.fh"
#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer fwater(2)


*     **** allocate fwater force ***
      value = BA_push_get(mt_dbl,(9*nwater),
     >                     'fwater',fwater(2),fwater(1))
      if (.not.value)
     > call errquit('Waterpsp_SD_Update:error push stack',0, MA_ERR)


*     **** Get Forces ****
      call Waterpsp_Fion(fion)
      call Waterpsp_Fwater(n2ft3d,rgrid,rho,dv,dbl_mb(fwater(1)))
      deltawater = 0.0d0


*     **** steepest descent update of water ****
      if (algorithm.eq.0) then
      call Waterpsp_SD_subUpdate(nwater,dbl_mb(rwater2(1)),
     >                                  dbl_mb(rwater1(1)),
     >                                  dbl_mb(fwater(1)),
     >                                  dt,
     >                                  deltawater)

*     **** Newton step update of water ****
      else if (algorithm.eq.1) then
        call Waterpsp_Newton_subUpdate(nwater,dbl_mb(rwater2(1)),
     >                                  dbl_mb(rwater1(1)),
     >                                  dbl_mb(rwater0(1)),
     >                                  dbl_mb(fwater(1)),
     >                                  dt,ekw)
        ekw_total = ekw_total + ekw
        ekw_count = ekw_count + 1

*     **** Verlet step update of water ****
      else if (algorithm.eq.2) then
        call Waterpsp_Verlet_subUpdate(nwater,dbl_mb(rwater2(1)),
     >                                  dbl_mb(rwater1(1)),
     >                                  dbl_mb(rwater0(1)),
     >                                  dbl_mb(fwater(1)),
     >                                  dt,ekw)
        ekw_total = ekw_total + ekw
        ekw_count = ekw_count + 1
      end if


*     **** deallocate fwater force ***
      value = BA_pop_stack(fwater(2))
      if (.not.value)
     > call errquit('Waterpsp_SD_Update:error pop stack',0, MA_ERR)
      return
	  end



*     **********************************
*     *	                               *
*     *      Waterpsp_SD_subUpdate     *
*     *                                *
*     **********************************
      subroutine Waterpsp_SD_subUpdate(nwater,rwater2,rwater1,fwater,
     >                                 dt,deltawater)
      implicit none
      integer nwater
      real*8 rwater2(3,3,*)
      real*8 rwater1(3,3,*)
      real*8 fwater(3,3,*)
      real*8 dt
      real*8 deltawater

*     **** local variables ***
      integer j,ii
      real*8 dtO,dtH,sum

*     **** steepest descent update of water ****
      dtO = dt/dsqrt(16.0d0*1822.89d0)
      dtH = dt/dsqrt( 1.0d0*1822.89d0)
      do ii=1,nwater
        rwater2(1,1,ii) = rwater1(1,1,ii) + dtO*fwater(1,1,ii)
        rwater2(2,1,ii) = rwater1(2,1,ii) + dtO*fwater(2,1,ii)
        rwater2(3,1,ii) = rwater1(3,1,ii) + dtO*fwater(3,1,ii)

        rwater2(1,2,ii) = rwater1(1,2,ii) + dtH*fwater(1,2,ii)
        rwater2(2,2,ii) = rwater1(2,2,ii) + dtH*fwater(2,2,ii)
        rwater2(3,2,ii) = rwater1(3,2,ii) + dtH*fwater(3,2,ii)
      
        rwater2(1,3,ii) = rwater1(1,3,ii) + dtH*fwater(1,3,ii)
        rwater2(2,3,ii) = rwater1(2,3,ii) + dtH*fwater(2,3,ii)
        rwater2(3,3,ii) = rwater1(3,3,ii) + dtH*fwater(3,3,ii)
      end do

*     **** find maximum force ****
      deltawater = 0.0d0
      do ii=1,nwater
      do j=1,3
        sum = dsqrt(fwater(1,j,ii)**2
     >             +fwater(2,j,ii)**2
     >             +fwater(3,j,ii)**2)
        if (sum.gt.deltawater) deltawater = sum
      end do
      end do

      return
      end

*     **********************************
*     *	                               *
*     *   Waterpsp_Newton_subUpdate    *
*     *                                *
*     **********************************
      subroutine Waterpsp_Newton_subUpdate(nwater,rwater2,rwater1,
     >                                            vwater,fwater,dt,ekw)
      implicit none
      integer nwater
      real*8 rwater2(3,3,*)
      real*8 rwater1(3,3,*)
      real*8 vwater(3,3,*)
      real*8 fwater(3,3,*)
      real*8 dt,ekw

*     **** local variables ****
      integer ii
      real*8  dtO,dtH,MO,MH

*     **** steepest descent update of water ****
      dtO = 0.5d0*dt*dt/(16.0d0*1822.89d0)
      dtH = 0.5d0*dt*dt/( 1.0d0*1822.89d0)
      do ii=1,nwater
        rwater2(1,1,ii) = rwater1(1,1,ii) + dt *vwater(1,1,ii) 
     >                                    + dtO*fwater(1,1,ii)
        rwater2(2,1,ii) = rwater1(2,1,ii) + dt *vwater(2,1,ii) 
     >                                    + dtO*fwater(2,1,ii)
        rwater2(3,1,ii) = rwater1(3,1,ii) + dt *vwater(3,1,ii) 
     >                                    + dtO*fwater(3,1,ii)

        rwater2(1,2,ii) = rwater1(1,2,ii) + dt *vwater(1,2,ii) 
     >                                    + dtH*fwater(1,2,ii)
        rwater2(2,2,ii) = rwater1(2,2,ii) + dt *vwater(2,2,ii) 
     >                                    + dtH*fwater(2,2,ii)
        rwater2(3,2,ii) = rwater1(3,2,ii) + dt *vwater(3,2,ii) 
     >                                    + dtH*fwater(3,2,ii)

        rwater2(1,3,ii) = rwater1(1,3,ii) + dt *vwater(1,3,ii) 
     >                                    + dtH*fwater(1,3,ii)
        rwater2(2,3,ii) = rwater1(2,3,ii) + dt *vwater(2,3,ii) 
     >                                    + dtH*fwater(2,3,ii)
        rwater2(3,3,ii) = rwater1(3,3,ii) + dt *vwater(3,3,ii) 
     >                                    + dtH*fwater(3,3,ii)
      end do

*     ***** find kinetic energy ****
      ekw = 0.0d0
      MO = (16.0d0*1822.89d0)
      MH = ( 1.0d0*1822.89d0)
      do ii=1,nwater
         ekw = ekw + MO*(vwater(1,1,ii)**2
     >                  +vwater(2,1,ii)**2
     >                  +vwater(3,1,ii)**2)
     >             + MH*(vwater(1,2,ii)**2
     >                  +vwater(2,2,ii)**2
     >                  +vwater(3,2,ii)**2)
     >             + MH*(vwater(1,3,ii)**2
     >                  +vwater(2,3,ii)**2
     >                  +vwater(3,3,ii)**2)
      end do
      ekw = 0.5d0*ekw

      return
      end


*     **********************************
*     *	                               *
*     *   Waterpsp_Verlet_subUpdate    *
*     *                                *
*     **********************************
      subroutine Waterpsp_Verlet_subUpdate(nwater,
     >                                 rwater2,rwater1,rwater0,fwater,
     >                                 dt,ekw)
      implicit none
      integer nwater
      real*8 rwater2(3,3,*)
      real*8 rwater1(3,3,*)
      real*8 rwater0(3,3,*)
      real*8 fwater(3,3,*)
      real*8 dt,ekw

*     **** local variables ****
      integer i,j,ii
      real*8 h,dtO,dtH,MO,MH

*     **** steepest descent update of water ****
      dtO = dt*dt/(16.0d0*1822.89d0)
      dtH = dt*dt/( 1.0d0*1822.89d0)
      do ii=1,nwater
        rwater2(1,1,ii) = 2*rwater1(1,1,ii) 
     >                  -   rwater0(1,1,ii) + dtO*fwater(1,1,ii)
        rwater2(2,1,ii) = 2*rwater1(2,1,ii) 
     >                  -   rwater0(2,1,ii) + dtO*fwater(2,1,ii)
        rwater2(3,1,ii) = 2*rwater1(3,1,ii) 
     >                  -   rwater0(3,1,ii) + dtO*fwater(3,1,ii)

        rwater2(1,2,ii) = 2*rwater1(1,2,ii) 
     >                  -   rwater0(1,2,ii) + dtH*fwater(1,2,ii)
        rwater2(2,2,ii) = 2*rwater1(2,2,ii) 
     >                  -   rwater0(2,2,ii) + dtH*fwater(2,2,ii)
        rwater2(3,2,ii) = 2*rwater1(3,2,ii) 
     >                  -   rwater0(3,2,ii) + dtH*fwater(3,2,ii)

        rwater2(1,3,ii) = 2*rwater1(1,3,ii) 
     >                  -   rwater0(1,3,ii) + dtH*fwater(1,3,ii)
        rwater2(2,3,ii) = 2*rwater1(2,3,ii) 
     >                  -   rwater0(2,3,ii) + dtH*fwater(2,3,ii)
        rwater2(3,3,ii) = 2*rwater1(3,3,ii) 
     >                  -   rwater0(3,3,ii) + dtH*fwater(3,3,ii)
      end do

*     **** make rwater0 the velocity - note that velocity is deleted after ****
*     **** Waterpsp_shift call                                             ****
      h = 1.0d0/(2.0d0*dt)
      do ii=1,nwater
      do j=1,3
      do i=1,3
         rwater0(i,j,ii) = h*(rwater2(i,j,ii) - rwater0(i,j,ii))
      end do
      end do
      end do

*     ***** find kinetic energy ****
      ekw = 0.0d0
      MO = (16.0d0*1822.89d0)
      MH = ( 1.0d0*1822.89d0)
      do ii=1,nwater
         ekw = ekw + MO*(rwater0(1,1,ii)**2
     >                  +rwater0(2,1,ii)**2
     >                  +rwater0(3,1,ii)**2)
     >             + MH*(rwater0(1,2,ii)**2
     >                  +rwater0(2,2,ii)**2
     >                  +rwater0(3,2,ii)**2)
     >             + MH*(rwater0(1,3,ii)**2
     >                  +rwater0(2,3,ii)**2
     >                  +rwater0(3,3,ii)**2)
      end do
      ekw = 0.5d0*ekw

      return
      end
c $Id$
