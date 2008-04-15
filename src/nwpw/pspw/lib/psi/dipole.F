*
* $Id: dipole.F,v 1.10 2008-04-15 22:59:18 bylaska Exp $
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
#include "errquit.fh"
#include "util.fh"


*     **** local variables ****
      logical value,oprint
      integer ii
      integer nx,ny,nz
      integer n1,n2,n3,ncut
      real*8 GX,GY,GZ
      real*8 qGX,qGY,qGZ
      real*8 cqGX,cqGY,cqGZ,x,y,z,r,rmax
      real*8 cdx1,cdy1,cdz1
      real*8 cdx2,cdy2,cdz2
      real*8 cdx3,cdy3,cdz3
      real*8 tmass,tcharge,ncharge,pcharge
      real*8 dv
      real*8 dipole_crystal(3),dipole_molecule(3)
      real*8 dtmp(3)

      integer rgrid(2)
      integer rgx(2),rgy(2),rgz(2)
      
      integer taskid,MASTER
      parameter (MASTER=0)

      real*8 autoDebye
      parameter (autoDebye=2.5416d0)
 
*     **** external functions ****
      logical  control_print
      integer  ion_katm,ion_nion,control_ncut,control_version
      real*8   ion_amass,psp_zv,ion_rion,lattice_omega
      real*8   lattice_unita
      external control_print
      external ion_katm,ion_nion,control_ncut,control_version
      external ion_amass,psp_zv,ion_rion,lattice_omega
      external lattice_unita

      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))


*     ***** center of mass  ****
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

      !*** crystal center of ionic charge ***
      ncut = 20
      n1 = ncut-2
      n2 = ncut-2
      n3 = ncut-2

      x = n1*lattice_unita(1,1)
      y = n1*lattice_unita(2,1)
      z = n1*lattice_unita(3,1)
      rmax = dsqrt(x*x + y*y + z*z) 

      x = n2*lattice_unita(1,2)
      y = n2*lattice_unita(2,2)
      z = n2*lattice_unita(3,2)
      r = dsqrt(x*x + y*y + z*z) 
      if (r.lt.rmax) rmax = r

      x = n3*lattice_unita(1,3)
      y = n3*lattice_unita(2,3)
      z = n3*lattice_unita(3,3)
      r = dsqrt(x*x + y*y + z*z) 
      if (r.lt.rmax) rmax = r

      cqGX=0.0d0
      cqGY=0.0d0
      cqGZ=0.0d0
      tcharge=0.0d0
      do ii=1,ion_nion()

        do n3= -ncut, ncut
        do n2= -ncut, ncut
        do n1= -ncut, ncut
          x = ion_rion(1,ii) 
     >     + n1*lattice_unita(1,1)
     >     + n2*lattice_unita(1,2)
     >     + n3*lattice_unita(1,3)
          y = ion_rion(2,ii) 
     >     + n1*lattice_unita(2,1)
     >     + n2*lattice_unita(2,2)
     >     + n3*lattice_unita(2,3)
          z = ion_rion(3,ii) 
     >     + n1*lattice_unita(3,1)
     >     + n2*lattice_unita(3,2)
     >     + n3*lattice_unita(3,3)

          r = dsqrt(x*x+y*y+z*z)

          if (r.le.rmax) then
            cqGX=cqGX+psp_zv(ion_katm(ii))*x
            cqGY=cqGY+psp_zv(ion_katm(ii))*y
            cqGZ=cqGZ+psp_zv(ion_katm(ii))*z
            tcharge=tcharge+psp_zv(ion_katm(ii))
          end if
        end do
        end do
        end do
      END DO
      cqGX=cqGX/tcharge
      cqGY=cqGY/tcharge
      cqGZ=cqGZ/tcharge
  
 



      !*** molecular center of ionic charge ***
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
     >   call errquit('Calculate_Dipole: out of stack memory',0, MA_ERR)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv=lattice_omega()/dble(nx*ny*nz)
      call lattice_r_grid_sym(dbl_mb(rgrid(1)))
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+0),3,dbl_mb(rgx(1)),1)
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+1),3,dbl_mb(rgy(1)),1)
      call dcopy(n2ft3d,dbl_mb(rgrid(1)+2),3,dbl_mb(rgz(1)),1)
      call D3dB_r_Zero_Ends(1,dbl_mb(rgx(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(rgy(1)))
      call D3dB_r_Zero_Ends(1,dbl_mb(rgz(1)))

      call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,1),cdx1)
      call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,1),cdy1)
      call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,1),cdz1)
      cdx1 = cdx1*dv
      cdy1 = cdy1*dv
      cdz1 = cdz1*dv


*     *** check for ferromagnetic case ***
      if (ne(ispin).ne.0) then
        call D3dB_rr_dot(1,dbl_mb(rgx(1)),dn(1,ispin),cdx2)
        call D3dB_rr_dot(1,dbl_mb(rgy(1)),dn(1,ispin),cdy2)
        call D3dB_rr_dot(1,dbl_mb(rgz(1)),dn(1,ispin),cdz2)
        cdx2 = cdx2*dv
        cdy2 = cdy2*dv
        cdz2 = cdz2*dv
      else
       cdx2 = 0.0d0
       cdy2 = 0.0d0
       cdz2 = 0.0d0
      end if

      cdx3=cdx1+cdx2
      cdy3=cdy1+cdy2
      cdz3=cdz1+cdz2

      call lattice_mask_sym(dbl_mb(rgrid(1)))
      !cdx1=cdx1/ne(1)
      !cdy1=cdy1/ne(1)
      !cdz1=cdz1/ne(1)
      call D3dB_rr_dot(1,dbl_mb(rgrid(1)),dn(1,1),rmax)
      rmax = rmax*dv
      cdx1=cdx1/rmax
      cdy1=cdy1/rmax
      cdz1=cdz1/rmax
      if (ne(ispin).ne.0) then
        !cdx2=cdx2/ne(ispin)
        !cdy2=cdy2/ne(ispin)
        !cdz2=cdz2/ne(ispin)
        call D3dB_rr_dot(1,dbl_mb(rgrid(1)),dn(1,ispin),rmax)
        rmax = rmax*dv
        cdx2=cdx2/rmax
        cdy2=cdy2/rmax
        cdz2=cdz2/rmax

      end if
      !cdx3=cdx3/dble(ne(1)+ne(ispin))
      !cdy3=cdy3/dble(ne(1)+ne(ispin))
      !cdz3=cdz3/dble(ne(1)+ne(ispin))
      call D3dB_rr_Sum(1,dn(1,1),
     >                   dn(1,ispin),
     >                   dbl_mb(rgrid(1)+n2ft3d))
      call D3dB_rr_dot(1,dbl_mb(rgrid(1)),
     >                   dbl_mb(rgrid(1)+n2ft3d),
     >                    rmax)
      rmax = rmax*dv
      cdx3=cdx3/rmax
      cdy3=cdy3/rmax
      cdz3=cdz3/rmax

      if (oprint) then
        WRITE(*,1200)
        WRITE(*,1220) 'spin up   ',CDX1,CDY1,CDZ1
        if (ne(ispin).ne.0) 
     >    WRITE(*,1220) 'spin down ',CDX2,CDY2,CDZ2
        WRITE(*,1220) '     total',CDX3,CDY3,CDZ3
        WRITE(*,1220) 'ionic     ',qGX,qGY,qGZ
        WRITE(*,1220) 'crystal   ',cqGX,cqGY,cqGZ
      end if
      dtmp(1) = CDX3
      dtmp(2) = CDY3
      dtmp(3) = CDZ3
      call ecce_print1('total dipole',mt_dbl,dtmp,3)
      dtmp(1) = CDX1
      dtmp(2) = CDY1
      dtmp(3) = CDZ1
      call ecce_print1('alpha dipole',mt_dbl,dtmp,3)
      if (ne(ispin).ne.0) then
         dtmp(1) = CDX2
         dtmp(2) = CDY2
         dtmp(3) = CDZ2
         call ecce_print1('beta dipole',mt_dbl,dtmp,3)
      endif
      dtmp(1) = qGX
      dtmp(2) = qGy
      dtmp(3) = qGz
      call ecce_print1('nuclear dipole',mt_dbl,dtmp,3)


 1200 FORMAT(//'== Center of Charge =='/)
 1220 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')

*     ***** calculate crystal dipole with respect to center of cell ****
      pcharge   = tcharge
      ncharge   = dble(ne(1)+ne(ispin))
      dipole_crystal(1) = -ncharge*cdx3 + pcharge*cqGX 
      dipole_crystal(2) = -ncharge*cdy3 + pcharge*cqGY 
      dipole_crystal(3) = -ncharge*cdz3 + pcharge*cqGZ 
      cdx1 = dsqrt( dipole_crystal(1)**2 
     >            + dipole_crystal(2)**2 
     >            + dipole_crystal(3)**2)
      if (oprint) then
         WRITE(*,1240)
         WRITE(*,1231) dipole_crystal
         WRITE(*,1232) cdx1,cdx1*autoDebye
      end if

*     ***** calculate dipole with respect to center of mass ****
      pcharge   = tcharge
      ncharge   = dble(ne(1)+ne(ispin))
      dipole_molecule(1) = -ncharge*cdx3 + pcharge*qGX 
     >                     - GX*(pcharge-ncharge)      
      dipole_molecule(2) = -ncharge*cdy3 + pcharge*qGY 
     >                     - GY*(pcharge-ncharge)
      dipole_molecule(3) = -ncharge*cdz3 + pcharge*qGZ 
     >                     - GZ*(pcharge-ncharge)
      cdx1 = dsqrt( dipole_molecule(1)**2 
     >            + dipole_molecule(2)**2 
     >            + dipole_molecule(3)**2)
      if (oprint) then
         WRITE(*,1230)
         WRITE(*,1231) dipole_molecule
         WRITE(*,1232) cdx1,cdx1*autoDebye
      end if
 1230 FORMAT(//'== Molecular Dipole wrt Center of Mass =='/)
 1231 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1232 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1240 FORMAT(//'== Crystal Dipole =='/)

*     **** pop stack memory ****
      value = value.and.MA_pop_stack(rgz(2))
      value = value.and.MA_pop_stack(rgy(2))
      value = value.and.MA_pop_stack(rgx(2))
      value = value.and.MA_pop_stack(rgrid(2))
      if (.not. value) 
     >   call errquit('Calculate_Dipole: cannot pop stack memory',0,
     &       MA_ERR)

      if (control_version().eq.3) then
         dipole(1) = dipole_crystal(1)
         dipole(2) = dipole_crystal(2)
         dipole(3) = dipole_crystal(3)
      else
         dipole(1) = dipole_molecule(1)
         dipole(2) = dipole_molecule(2)
         dipole(3) = dipole_molecule(3)
      end if

      return
      end



