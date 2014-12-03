c
c $Id$
c

*     ************************************
*     *                    	         *
*     *    pspw_Lubin_water_analysis     *
*     *                                  *
*     ************************************
      subroutine pspw_Lubin_water_analysis(rtdb,ispin,ne,n2ft3d,dn)
      implicit none
      integer rtdb
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer ma_type,nactive_atoms,h_actlist,l_actlist
      integer ii,nwater


*     **** pspw:Lubin set ****
      if (rtdb_ma_get(rtdb, 'pspw:Lubin', ma_type,
     >        nactive_atoms, h_actlist)) then

*       **** error check ***
        if (mod(nactive_atoms,3).ne.0) 
     >   call errquit(
     >     'pspw_Lubin_Water_analysis: atoms not multiple of 3',911,
     &       INPUT_ERR)

         if (.not.BA_get_index(h_actlist,l_actlist))
     >   call errquit(
     >     'pspw_Lubin_Water_analysis: ma_get_index failed',911, MA_ERR)

         nwater = nactive_atoms/3
         do ii=1,nwater
           call pspw_sub_Lubin(ispin,ne,n2ft3d,dn,
     >                         int_mb(l_actlist + 3*(ii-1)))
         end do

*       **** free heap ****
        if (.not. BA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

      end if

      return
      end

*     ************************************
*     *                    		         *
*     *    pspw_sub_Lubin                *
*     *                            		 *
*     ************************************
      subroutine pspw_sub_Lubin(ispin,ne,n2ft3d,dn,
     >                                     indx)
      implicit none
#include "errquit.fh"
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)
      integer indx(3)

#include "bafdecls.fh"


*     **** Lubin Parameters ****
      real*8 rcut,alpha
      parameter (rcut=4.0d0,alpha=0.1d0)

*     **** Parameters ****
      real*8 autoDebye
      parameter (autoDebye=2.5416d0)

      integer taskid,MASTER
      parameter (MASTER=0)


*     **** local variables ****
      logical value
      integer k,ii
      integer nx,ny,nz
      real*8 GX,GY,GZ
      real*8 qGX,qGY,qGZ
      real*8 cdx1,cdy1,cdz1
      real*8 cdx2,cdy2,cdz2
      real*8 cdx3,cdy3,cdz3
      real*8 tmass,tcharge,ncharge
      real*8 x,dist
      real*8 dv
      real*8 amass(3),zv(3),dipole(3)

*     **** MA stack memory ****
      integer rgrid(2)
      integer rgx(2),rgy(2),rgz(2),fcut(2)
      

 
*     **** external functions ****
      real*8   lattice_omega,util_erfc,ion_rion
      external lattice_omega,util_erfc,ion_rion

      call Parallel_taskid(taskid)


*     ***** calculate the center of mass and the center of ionic charge ****
      amass(1) = 16.0d0
      amass(2) =  1.0d0
      amass(3) =  1.0d0
      zv(1) = 6.0d0
      zv(2) = 1.0d0
      zv(3) = 1.0d0
      GX=0.0d0
      GY=0.0d0
      GZ=0.0d0
      tmass=0.0d0
      DO ii=1,3
        tmass=tmass+amass(ii)
        GX=GX+amass(ii)*ion_rion(1,indx(ii))
        GY=GY+amass(ii)*ion_rion(2,indx(ii))
        GZ=GZ+amass(ii)*ion_rion(3,indx(ii))
      END DO
      GX=GX/tmass
      GY=GY/tmass
      GZ=GZ/tmass
      qGX=0.0d0
      qGY=0.0d0
      qGZ=0.0d0
      tcharge=0.0d0
      DO ii=1,3
        tcharge=tcharge+zv(ii)
        qGX=qGX+zv(ii)*ion_rion(1,indx(ii))
        qGY=qGY+zv(ii)*ion_rion(2,indx(ii))
        qGZ=qGZ+zv(ii)*ion_rion(3,indx(ii))
      END DO
      qGX=qGX/tcharge
      qGY=qGY/tcharge
      qGZ=qGZ/tcharge




*     **** calculate the center of density ****
      value = BA_push_get(mt_dbl,3*n2ft3d,'rgrid', rgrid(2), rgrid(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgx',rgx(2),rgx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgy',rgy(2),rgy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgz',rgz(2),rgz(1))     
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fcut',fcut(2),fcut(1))     
      if (.not. value) 
     >   call errquit('Lubin_analysis: out of stack memory',0, MA_ERR)

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

*     **** define fcut = 0.5*erfc((|r-RO|-rcut)/alpha) ****
      do k=1,n2ft3d
         x = (dbl_mb(rgx(1)+k-1) - ion_rion(1,indx(1)))**2
     >     + (dbl_mb(rgy(1)+k-1) - ion_rion(2,indx(1)))**2
     >     + (dbl_mb(rgz(1)+k-1) - ion_rion(3,indx(1)))**2
         x=dsqrt(x)
         x = (x-rcut)/alpha
         dbl_mb(fcut(1)+k-1) = 0.5d0*util_erfc(x)
      end do

c      call D3dB_rr_Mul(1,dbl_mb(rgx(1)),dbl_mb(fcut(1)),dbl_mb(rgx(1)))
c      call D3dB_rr_Mul(1,dbl_mb(rgy(1)),dbl_mb(fcut(1)),dbl_mb(rgy(1)))
c      call D3dB_rr_Mul(1,dbl_mb(rgz(1)),dbl_mb(fcut(1)),dbl_mb(rgz(1)))
      call D3dB_rr_Mul2(1,dbl_mb(fcut(1)),dbl_mb(rgx(1)))
      call D3dB_rr_Mul2(1,dbl_mb(fcut(1)),dbl_mb(rgy(1)))
      call D3dB_rr_Mul2(1,dbl_mb(fcut(1)),dbl_mb(rgz(1)))

      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,1),cdx1)
      call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,1),cdy1)
      call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,1),cdz1)
      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,ispin),cdx2)
      call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,ispin),cdy2)
      call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,ispin),cdz2)
      cdx3=(cdx1+cdx2)*dv
      cdy3=(cdy1+cdy2)*dv
      cdz3=(cdz1+cdz2)*dv

*     **** calculate ncharge ****
      call D3dB_r_Zero(1,dbl_mb(rgx(1)))
      call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rgx(1)))
      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dbl_mb(fcut(1)),ncharge)
      ncharge = ncharge*dv

*     **** calculate center of charge ****
      cdx3=cdx3/ncharge
      cdy3=cdy3/ncharge
      cdz3=cdz3/ncharge

      if (taskid.eq.MASTER) then
         WRITE(*,100)
         WRITE(*,110)
         WRITE(*,111) rcut
         WRITE(*,112) alpha
         WRITE(*,115)
         WRITE(*,120) ion_rion(1,indx(1)),
     >                ion_rion(2,indx(1)),
     >                ion_rion(3,indx(1))
         WRITE(*,125) ion_rion(1,indx(2)),
     >                ion_rion(2,indx(2)),
     >                ion_rion(3,indx(2))
         WRITE(*,130) ion_rion(1,indx(3)),
     >                ion_rion(2,indx(3)),
     >                ion_rion(3,indx(3))
         WRITE(*,131) ncharge

         WRITE(*,135)
         WRITE(*,140) cdx3,cdy3,cdz3
         WRITE(*,145) qGx, qGy, qGz
         dist = (cdx3-qGx)**2 + (cdy3-qGy)**2 + (cdz3-qGz)**2
         dist = dsqrt(dist)
         write(*,146) dist


*        ***** calculate dipole relative to CM ****
c        ncharge   = 8
         dipole(1) = -ncharge*cdx3 + 8.0d0*qGX -GX*(8.0d0-ncharge)
         dipole(2) = -ncharge*cdy3 + 8.0d0*qGY -GY*(8.0d0-ncharge)
         dipole(3) = -ncharge*cdz3 + 8.0d0*qGZ -GZ*(8.0d0-ncharge)
         cdx1 = dsqrt(dipole(1)**2 + dipole(2)**2 + dipole(3)**2)

         WRITE(*,1230)
         WRITE(*,1231) dipole
         WRITE(*,1232) cdx1,cdx1*autoDebye

         ncharge   = 8
         dipole(1) = -ncharge*cdx3 + 8.0d0*qGX -GX*(8.0d0-ncharge)
         dipole(2) = -ncharge*cdy3 + 8.0d0*qGY -GY*(8.0d0-ncharge)
         dipole(3) = -ncharge*cdz3 + 8.0d0*qGZ -GZ*(8.0d0-ncharge)
         cdx1 = dsqrt(dipole(1)**2 + dipole(2)**2 + dipole(3)**2)

         WRITE(*,1240)
         WRITE(*,1231) dipole
         WRITE(*,1232) cdx1,cdx1*autoDebye

      end if

*     **** pop stack memory ****
      value = value.and.BA_pop_stack(fcut(2))
      value = value.and.BA_pop_stack(rgz(2))
      value = value.and.BA_pop_stack(rgy(2))
      value = value.and.BA_pop_stack(rgx(2))
      value = value.and.BA_pop_stack(rgrid(2))
      if (.not. value) 
     >   call errquit('Lubin analysis: cannot pop stack memory',0,
     &       MA_ERR)

      return

  100 format(//'Lubin Water Analysis')
  110 format( '--------------------')
  111 format('Sphere radius     :',F10.4)
  112 format('Boundary sharpness:',F10.4)

  115 format(/'Water position:') 
  120 format('RO           ','  (',F10.4,',',F10.4,',',F10.4,' )')
  125 format('RH1          ','  (',F10.4,',',F10.4,',',F10.4,' )')
  130 format('RH2          ','  (',F10.4,',',F10.4,',',F10.4,' )')
  131 format('Qe = ',F10.4,' au')  

  135 format(/'Center of charge:')
  140 format('water density','  (',F10.4,',',F10.4,',',F10.4,' )')
  145 format('water ions   ','  (',F10.4,',',F10.4,',',F10.4,' )')
  146 format('charge distance:',F10.4)

 1230 FORMAT(/'Water Dipole:')
 1231 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1232 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1240 FORMAT(/'Water Dipole,ncharge=8.0:')
      end



