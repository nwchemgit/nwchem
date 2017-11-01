*
* $Id$
*

      subroutine center_geom(cx,cy,cz)
      implicit none
      real*8 cx,cy,cz

*     **** local variables ***
      integer i,nion

*     **** external functions ****
      integer  ion_nion
      real*8   dsum,ion_rion
      external ion_nion
      external dsum,ion_rion

*:::::::::::::::::  geometrical center of the cluster  ::::::::::::::::
      nion = ion_nion()
      cx =0.0d0
      cy =0.0d0
      cz =0.0d0
      do i=1,nion
         cx = cx + ion_rion(1,i)
         cy = cy + ion_rion(2,i)
         cz = cz + ion_rion(3,i)
      end do
      cx = cx/dble(nion)
      cy = cy/dble(nion)
      cz = cz/dble(nion)

      return
      end


      subroutine center_mass(gx,gy,gz)
      implicit none
      real*8   gx,gy,gz

*     **** local variables ****
      integer nion
      integer i
      real*8 am

*     **** external functions ****
      integer  ion_nion
      real*8   ion_amass,ion_rion
      external ion_nion
      external ion_amass,ion_rion

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+ion_amass(i)*ion_rion(1,i)
        gy=gy+ion_amass(i)*ion_rion(2,i)
        gz=gz+ion_amass(i)*ion_rion(3,i)
        am=am+ion_amass(i)
      end do
      gx=gx/am
      gy=gy/am
      gz=gz/am

      return 
      end



      subroutine center_v_geom(cx,cy,cz)
      implicit none
      real*8 cx,cy,cz

*     **** local variables ***
      integer i,nion

*     **** external functions ****
      integer  ion_nion
      real*8   dsum,ion_vion
      external ion_nion
      external dsum,ion_vion

*:::::::::::::::::  geometrical center of the cluster  ::::::::::::::::
      nion = ion_nion()
      cx =0.0d0
      cy =0.0d0
      cz =0.0d0
      do i=1,nion
         cx = cx + ion_vion(1,i)
         cy = cy + ion_vion(2,i)
         cz = cz + ion_vion(3,i)
      end do
      cx = cx/dble(nion)
      cy = cy/dble(nion)
      cz = cz/dble(nion)

      return
      end


      subroutine center_v_mass(gx,gy,gz)
      implicit none
      real*8   gx,gy,gz

*     **** local variables ****
      integer nion
      integer i
      real*8 am

*     **** external functions ****
      integer  ion_nion
      real*8   ion_amass,ion_vion
      external ion_nion
      external ion_amass,ion_vion

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+ion_amass(i)*ion_vion(1,i)
        gy=gy+ion_amass(i)*ion_vion(2,i)
        gz=gz+ion_amass(i)*ion_vion(3,i)
        am=am+ion_amass(i)
      end do
      gx=gx/am
      gy=gy/am
      gz=gz/am

      return 
      end


      subroutine center_F_mass(F,gx,gy,gz)
      implicit none
      real*8   F(3,*)
      real*8   gx,gy,gz

*     **** local variables ****
      integer nion
      integer i,i1
      real*8 am

*     **** external functions ****
      integer  ion_nion
      real*8   ion_amass
      external ion_nion
      external ion_amass

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+ion_amass(i)*F(1,i)
        gy=gy+ion_amass(i)*F(2,i)
        gz=gz+ion_amass(i)*F(3,i)
        am=am+ion_amass(i)
      end do
      gx=gx/am
      gy=gy/am
      gz=gz/am

      return 
      end

      subroutine remove_center_F_mass(F)
      implicit none
      real*8   F(3,*)

*     **** local variables ****
      integer nion
      integer i,i1
      real*8 am
      real*8   gx,gy,gz

*     **** external functions ****
      integer  ion_nion
      real*8   ion_amass
      external ion_nion
      external ion_amass

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+ion_amass(i)*F(1,i)
        gy=gy+ion_amass(i)*F(2,i)
        gz=gz+ion_amass(i)*F(3,i)
        am=am+ion_amass(i)
      end do
      gx=gx/am
      gy=gy/am
      gz=gz/am

      !**** remove center of mass motion ***
!$OMP DO
      do i=1,nion
         F(1,i) = F(1,i) - gx
         F(2,i) = F(2,i) - gy
         F(3,i) = F(3,i) - gz
      end do
!$OMP END DO

      return 
      end



      subroutine remove_center_mass(R2,R1)
      implicit none
      real*8   R2(3,*)
      real*8   R1(3,*)

*     **** local variables ****
      integer nion
      integer i,i1
      real*8   gx,gy,gz
      real*8   hx,hy,hz
      real*8 am0,gx0,gy0,gz0,hx0,hy0,hz0
      common /removecenter_shrd/ am0,gx0,gy0,gz0,hx0,hy0,hz0

*     **** external functions ****
      integer  ion_nion
      real*8   ion_amass
      external ion_nion
      external ion_amass

      nion = ion_nion()
!$OMP MASTER
      hx0 = 0.0d0
      hy0 = 0.0d0
      hz0 = 0.0d0
      gx0 = 0.0d0
      gy0 = 0.0d0
      gz0 = 0.0d0
      am0 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO private(i) reduction(+:am0,gx0,gy0,gz0,hx0,hy0,hz0)
      do i=1,nion
        am0=am0+ion_amass(i)
        gx0=gx0+ion_amass(i)*R1(1,i)
        gy0=gy0+ion_amass(i)*R1(2,i)
        gz0=gz0+ion_amass(i)*R1(3,i)
        hx0=hx0+ion_amass(i)*R2(1,i)
        hy0=hy0+ion_amass(i)*R2(2,i)
        hz0=hz0+ion_amass(i)*R2(3,i)
      end do
!$OMP END DO
!$OMP BARRIER
      hx=hx0/am0
      hy=hy0/am0
      hz=hz0/am0
      gx=gx0/am0
      gy=gy0/am0
      gz=gz0/am0

      !**** remove center of mass motion ***
!$OMP DO private(i)
      do i=1,nion
         R2(1,i) = R2(1,i) - hx + gx
         R2(2,i) = R2(2,i) - hy + gy
         R2(3,i) = R2(3,i) - hz + gz
      end do
!$OMP END DO

      return
      end

