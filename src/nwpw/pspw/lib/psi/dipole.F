*
* $Id: dipole.F,v 1.1 2002-01-14 21:06:28 bylaska Exp $
*

*     ******************************
*     *                            *
*     *    Calculate_Dipole        *
*     *                            *
*     ******************************

      subroutine Calculate_Dipole(ispin,ne,n2ft3d,dn,dipole)
      implicit none
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)
      real*8 dipole(3)

#include "mafdecls.fh"


*     **** local variables ****
      logical value
      integer ii
      integer nx,ny,nz
      real*8 GX,GY,GZ
      real*8 qGX,qGY,qGZ
      real*8 cdx1,cdy1,cdz1
      real*8 cdx2,cdy2,cdz2
      real*8 cdx3,cdy3,cdz3
      real*8 tmass,tcharge,ncharge,pcharge
      real*8 dv

      integer rgrid(2)
      integer rgx(2),rgy(2),rgz(2)
      
      integer taskid,MASTER
      parameter (MASTER=0)

      real*8 autoDebye
      parameter (autoDebye=2.5416d0)
 
*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   ion_amass,psp_zv,ion_rion,lattice_omega
      external ion_katm,ion_nion
      external ion_amass,psp_zv,ion_rion,lattice_omega

      call Parallel_taskid(taskid)

*     ***** calculate the center of mass and the center of ionic charge ****
      GX=0.0d0
      GY=0.0d0
      GZ=0.0d0
      tmass=0.0d0
      DO ii=1,ion_nion()
        tmass=tmass+ion_amass(ii)
        GX=GX+ion_amass(ii)*ion_rion(1,ii)
        GY=GY+ion_amass(ii)*ion_rion(2,ii)
        GZ=GZ+ion_amass(ii)*ion_rion(3,ii)
      END DO
      GX=GX/tmass
      GY=GY/tmass
      GZ=GZ/tmass
      qGX=0.0d0
      qGY=0.0d0
      qGZ=0.0d0
      tcharge=0.0d0
      DO ii=1,ion_nion()
        tcharge=tcharge+psp_zv(ion_katm(ii))
        qGX=qGX+psp_zv(ion_katm(ii))*ion_rion(1,ii)
        qGY=qGY+psp_zv(ion_katm(ii))*ion_rion(2,ii)
        qGZ=qGZ+psp_zv(ion_katm(ii))*ion_rion(3,ii)
      END DO
      qGX=qGX/tcharge
      qGY=qGY/tcharge
      qGZ=qGZ/tcharge

*     **** calculate the center of density ****
      value = MA_push_get(mt_dbl,3*n2ft3d,'rgrid', rgrid(2), rgrid(1))
      value = value.and.
     >        MA_push_get(mt_dbl, n2ft3d,'rgx',rgx(2),rgx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, n2ft3d,'rgy',rgy(2),rgy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, n2ft3d,'rgz',rgz(2),rgz(1))     
      if (.not. value) 
     >   call errquit('Calculate_Dipole: out of stack memory',0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv=lattice_omega()/dble(nx*ny*nz)
      call lattice_r_grid(dbl_mb(rgrid(1)))
      call D3dB_r_Zero(1,dbl_mb(rgx(1)))
      call D3dB_r_Zero(1,dbl_mb(rgy(1)))
      call D3dB_r_Zero(1,dbl_mb(rgz(1)))
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+0),3,dbl_mb(rgx(1)),1)
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+1),3,dbl_mb(rgy(1)),1)
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+2),3,dbl_mb(rgz(1)),1)

      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,1),cdx1)
      call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,1),cdy1)
      call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,1),cdz1)
      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,ispin),cdx2)
      call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,ispin),cdy2)
      call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,ispin),cdz2)
      cdx1 = cdx1*dv
      cdy1 = cdy1*dv
      cdz1 = cdz1*dv
      cdx2 = cdx2*dv
      cdy2 = cdy2*dv
      cdz2 = cdz2*dv

      cdx3=cdx1+cdx2
      cdy3=cdy1+cdy2
      cdz3=cdz1+cdz2
      cdx1=cdx1/ne(1)
      cdy1=cdy1/ne(1)
      cdz1=cdz1/ne(1)
      cdx2=cdx2/ne(ispin)
      cdy2=cdy2/ne(ispin)
      cdz2=cdz2/ne(ispin)
      cdx3=cdx3/dble(ne(1)+ne(ispin))
      cdy3=cdy3/dble(ne(1)+ne(ispin))
      cdz3=cdz3/dble(ne(1)+ne(ispin))

      if (taskid.eq.MASTER) then
        WRITE(*,1200)
        WRITE(*,1220) 'spin up   ',CDX1,CDY1,CDZ1
        WRITE(*,1220) 'spin down ',CDX2,CDY2,CDZ2
        WRITE(*,1220) '     total',CDX3,CDY3,CDZ3
        WRITE(*,1220) 'ionic     ',qGX,qGY,qGZ
      end if

 1200 FORMAT(//'== Center of Charge =='/)
 1220 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')

*     ***** calculate dipole with respect to center of mass ****
      pcharge   = tcharge
      ncharge   = dble(ne(1)+ne(ispin))
      dipole(1) = -ncharge*cdx3 + pcharge*qGX - GX*(pcharge-ncharge)      
      dipole(2) = -ncharge*cdy3 + pcharge*qGY - GY*(pcharge-ncharge)
      dipole(3) = -ncharge*cdz3 + pcharge*qGZ - GZ*(pcharge-ncharge)
      cdx1 = dsqrt(dipole(1)**2 + dipole(2)**2 + dipole(3)**2)
      if (taskid.eq.MASTER) then
         WRITE(*,1230)
         WRITE(*,1231) dipole
         WRITE(*,1232) cdx1,cdx1*autoDebye
      end if
 1230 FORMAT(//'== Molecular Dipole =='/)
 1231 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1232 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')

*     **** pop stack memory ****
      value = value.and.MA_pop_stack(rgz(2))
      value = value.and.MA_pop_stack(rgy(2))
      value = value.and.MA_pop_stack(rgx(2))
      value = value.and.MA_pop_stack(rgrid(2))
      if (.not. value) 
     >   call errquit('Calculate_Dipole: cannot pop stack memory',0)

      return
      end



