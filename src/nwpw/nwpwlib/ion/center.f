*
* $Id: center.f,v 1.1 2001-08-30 17:56:10 bylaska Exp $
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
      integer  ion_nion,ion_katm
      real*8   psp_amass,ion_rion
      external ion_nion,ion_katm
      external psp_amass,ion_rion

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+psp_amass(ion_katm(i))*ion_rion(1,i)
        gy=gy+psp_amass(ion_katm(i))*ion_rion(2,i)
        gz=gz+psp_amass(ion_katm(i))*ion_rion(3,i)
        am=am+psp_amass(ion_katm(i))
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
      integer  ion_nion,ion_katm
      real*8   psp_amass,ion_vion
      external ion_nion,ion_katm
      external psp_amass,ion_vion

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+psp_amass(ion_katm(i))*ion_vion(1,i)
        gy=gy+psp_amass(ion_katm(i))*ion_vion(2,i)
        gz=gz+psp_amass(ion_katm(i))*ion_vion(3,i)
        am=am+psp_amass(ion_katm(i))
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
      integer i,i1,i2,i3
      real*8 am

*     **** external functions ****
      logical  Waterpsp_found
      integer  Waterpsp_nwater
      integer  ion_nion,ion_katm
      real*8   psp_amass
      external Waterpsp_found
      external Waterpsp_nwater
      external ion_nion,ion_katm
      external psp_amass

      nion = ion_nion()
      gx=0.0d0
      gy=0.0d0
      gz=0.0d0
      am=0.0d0
      do i=1,nion
        gx=gx+psp_amass(ion_katm(i))*F(1,i)
        gy=gy+psp_amass(ion_katm(i))*F(2,i)
        gz=gz+psp_amass(ion_katm(i))*F(3,i)
        am=am+psp_amass(ion_katm(i))
      end do
      if (Waterpsp_found()) then
         do i=1,Waterpsp_nwater()
            i1 = nion + 3*(i-1)+1
            i2 = i1+1
            i3 = i1+2
            gx = gx + (16*1822.89d0)*F(1,i1)
            gy = gy + (16*1822.89d0)*F(2,i1)
            gz = gz + (16*1822.89d0)*F(3,i1)
            am=am+(16*1822.89d0)
            gx = gx + (1 *1822.89d0)*F(1,i2)
            gy = gy + (1 *1822.89d0)*F(2,i2)
            gz = gz + (1 *1822.89d0)*F(3,i2)
            am=am+(1*1822.89d0)
            gx = gx + (1 *1822.89d0)*F(1,i3)
            gy = gy + (1 *1822.89d0)*F(2,i3)
            gz = gz + (1 *1822.89d0)*F(3,i3)
            am=am+(1*1822.89d0)
         end do
      end if
      gx=gx/am
      gy=gy/am
      gz=gz/am

      return 
      end
