*
* $Id: center.f,v 1.3 2001-10-26 02:28:22 bylaska Exp $
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
      logical  pspw_qmmm_found
      integer  pspw_qmmm_nion
      integer  ion_nion
      real*8   ion_amass,pspw_qmmm_amass
      external pspw_qmmm_found
      external pspw_qmmm_nion
      external ion_nion
      external ion_amass,pspw_qmmm_amass

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
      if (pspw_qmmm_found()) then
         do i=1,pspw_qmmm_nion()
            i1 = nion + i
            gx = gx + pspw_qmmm_amass(i)*F(1,i1)
            gy = gy + pspw_qmmm_amass(i)*F(2,i1)
            gz = gz + pspw_qmmm_amass(i)*F(3,i1)
            am=am+pspw_qmmm_amass(i)

         end do
      end if
      gx=gx/am
      gy=gy/am
      gz=gz/am

      return 
      end
