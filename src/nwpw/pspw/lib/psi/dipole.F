*
* $Id$
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"


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
      value = BA_push_get(mt_dbl,3*n2ft3d,'rgrid', rgrid(2), rgrid(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgx',rgx(2),rgx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgy',rgy(2),rgy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgz',rgz(2),rgz(1))     
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
        write(LuOut,1200)
        write(LuOut,1220) 'spin up   ',CDX1,CDY1,CDZ1
        if (ne(ispin).ne.0) 
     >    write(LuOut,1220) 'spin down ',CDX2,CDY2,CDZ2
        write(LuOut,1220) '     total',CDX3,CDY3,CDZ3
        write(LuOut,1220) 'ionic     ',qGX,qGY,qGZ
        !write(LuOut,1220) 'crystal   ',cqGX,cqGY,cqGZ
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

c*     ***** calculate crystal dipole with respect to center of cell ****
c      pcharge   = tcharge
c      ncharge   = dble(ne(1)+ne(ispin))
c      dipole_crystal(1) = -ncharge*cdx3 + pcharge*cqGX 
c      dipole_crystal(2) = -ncharge*cdy3 + pcharge*cqGY 
c      dipole_crystal(3) = -ncharge*cdz3 + pcharge*cqGZ 
c      cdx1 = dsqrt( dipole_crystal(1)**2 
c     >            + dipole_crystal(2)**2 
c     >            + dipole_crystal(3)**2)
c      if (oprint) then
c         write(LuOut,1240)
c         write(LuOut,1231) dipole_crystal
c         write(LuOut,1232) cdx1,cdx1*autoDebye
c      end if

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
         write(LuOut,1230)
         write(LuOut,1231) dipole_molecule
         write(LuOut,1232) cdx1,cdx1*autoDebye
      end if
 1230 FORMAT(//'== Molecular Dipole wrt Center of Mass =='/)
 1231 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1232 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1240 FORMAT(//'== Crystal Dipole =='/)

*     **** pop stack memory ****
      value = value.and.BA_pop_stack(rgz(2))
      value = value.and.BA_pop_stack(rgy(2))
      value = value.and.BA_pop_stack(rgx(2))
      value = value.and.BA_pop_stack(rgrid(2))
      if (.not. value) 
     >   call errquit('Calculate_Dipole: cannot pop stack memory',0,
     &       MA_ERR)

c      if (control_version().eq.3) then
c         dipole(1) = dipole_crystal(1)
c         dipole(2) = dipole_crystal(2)
c         dipole(3) = dipole_crystal(3)
c      else
c         dipole(1) = dipole_molecule(1)
c         dipole(2) = dipole_molecule(2)
c         dipole(3) = dipole_molecule(3)
c      end if
      dipole(1) = dipole_molecule(1)
      dipole(2) = dipole_molecule(2)
      dipole(3) = dipole_molecule(3)

      return
      end



*     ******************************
*     *                            *
*     * Calculate_Molecular_Dipole *
*     *                            *
*     ******************************

      subroutine Calculate_Molecular_Dipole(ispin,ne,n2ft3d,dn,dipole)
      implicit none
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,ispin)
      real*8 dipole(3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"


*     **** local variables ****
      logical value,oprint
      integer ii
      integer nx,ny,nz
      integer n1,n2,n3,ncut
      real*8 GX,GY,GZ
      real*8 qGX,qGY,qGZ
      real*8 x,y,z,r,rmax
      real*8 cdx1,cdy1,cdz1
      real*8 cdx2,cdy2,cdz2
      real*8 cdx3,cdy3,cdz3
      real*8 tmass,tcharge,ncharge,pcharge
      real*8 dv
      
      integer rgrid(2)
      integer rgx(2),rgy(2),rgz(2)
      
      real*8 autoDebye
      parameter (autoDebye=2.5416d0)
 
*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   ion_amass,psp_zv,ion_rion,lattice_omega
      external ion_katm,ion_nion
      external ion_amass,psp_zv,ion_rion,lattice_omega


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
      value = BA_push_get(mt_dbl,3*n2ft3d,'rgrid',rgrid(2),rgrid(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgx',rgx(2),rgx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgy',rgy(2),rgy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rgz',rgz(2),rgz(1))     
      if (.not. value) 
     > call errquit('Calculate_Molecular Dipole:out of stack',0,MA_ERR)

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


*     ***** calculate dipole with respect to center of mass ****
      pcharge   = tcharge
      ncharge   = dble(ne(1)+ne(ispin))
      dipole(1) = -ncharge*cdx3 + pcharge*qGX 
     >          - GX*(pcharge-ncharge)      
      dipole(2) = -ncharge*cdy3 + pcharge*qGY 
     >          - GY*(pcharge-ncharge)
      dipole(3) = -ncharge*cdz3 + pcharge*qGZ 
     >          - GZ*(pcharge-ncharge)


*     **** pop stack memory ****
      value =           BA_pop_stack(rgz(2))
      value = value.and.BA_pop_stack(rgy(2))
      value = value.and.BA_pop_stack(rgx(2))
      value = value.and.BA_pop_stack(rgrid(2))
      if (.not.value) 
     >call errquit('Calculate_Molecular_Dipole:cannot pop stack memory',
     >             0,MA_ERR)

      return
      end


*     ******************************
*     *                            *
*     *  Calculate_Resta_Dipole    *
*     *                            *
*     ******************************

      subroutine Calculate_Resta_Dipole(doprint,
     >                                  ispin,ne,neq,npack1,nfft3d,psi1,
     >                                  dipole)
      implicit none
      logical doprint
      integer ispin,ne(2),neq(2)
      integer npack1,nfft3d
      complex*16 psi1(npack1,*)
      real*8     dipole(3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid,tmp_len
      parameter (MASTER=0,tmp_len=140)

      real*8 autoDebye
      parameter (autoDebye=2.5416d0)

      logical value,oprint
      integer i,j,ms,n,n2ft3d
      integer shift,rank
      integer info
      integer ii

      integer X(2,6),Xeig(2,6),psi_r(2),psi_r2(2)
      real*8 bv(3,6),wrk(6,6),wts(6),bmat(3,3)
      real*8 b(3),ixmat(3,6)
      real*8 xx,yy,zz,tmp1(tmp_len),maxtime,pcharge
      real*8 dx(2),dy(2),dz(2),nx,ny,nz,tx,ty,tz,alpha,scal
      complex*16 arg,wx
      character*2 labels(6)

*     **** external functions ****
      integer  ion_katm,ion_nion
      external ion_katm,ion_nion
      real*8   ion_rion,psp_zv
      external ion_rion,psp_zv
      real*8   lattice_unitg,lattice_omega
      external lattice_unitg,lattice_omega
      logical  Dneall_w_push_get,Dneall_w_pop_stack,control_print
      external Dneall_w_push_get,Dneall_w_pop_stack,control_print


      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium)
     >        .and.doprint)

*     ***** allocate X,Y,Z  ****
      n2ft3d = 2*nfft3d
      n=ne(1)
      if (n.lt.ne(2)) n=ne(2)
      value = BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,
     >                     'psi_r2',psi_r2(2),psi_r2(1))
      do j=1,6
         value = value.and.Dneall_w_push_get(1,X(1,j))
         value = value.and.
     >           BA_push_get(mt_dcpl,n,'Xeig',Xeig(2,j),Xeig(1,j))
      end do
      if (.not. value) 
     >   call errquit('dipole_resta:out of stack',1,0)

*     **** transform psi1 to realspace ****
      do n=1,neq(1)+neq(2)
         call Pack_c_Copy(1,psi1(1,n),dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do
      call Grsm_gh_fftb(nfft3d,neq(1)+neq(2),dbl_mb(psi_r(1)))



c     *** Silvestrelli G1 ***
      ixmat(1,1)=1.0d0
      ixmat(2,1)=0.0d0
      ixmat(3,1)=0.0d0

c     *** Silvestrelli G4 ***
      ixmat(1,2)=1.0d0
      ixmat(2,2)=1.0d0
      ixmat(3,2)=0.0d0

c     *** Silvestrelli G5 ***
      ixmat(1,3)=1.0d0
      ixmat(2,3)=0.0d0
      ixmat(3,3)=1.0d0

c     *** Silvestrelli G2 ***
      ixmat(1,4)=0.0d0
      ixmat(2,4)=1.0d0
      ixmat(3,4)=0.0d0

c     *** Silvestrelli G6 ***
      ixmat(1,5)=0.0d0
      ixmat(2,5)=1.0d0
      ixmat(3,5)=1.0d0

c     *** Silvestrelli G3 ***
      ixmat(1,6)=0.0d0
      ixmat(2,6)=0.0d0
      ixmat(3,6)=1.0d0

      do i=1,3
         bmat(i,1)=lattice_unitg(1,i)
         bmat(i,2)=lattice_unitg(2,i)
         bmat(i,3)=lattice_unitg(3,i)
      end do

      do i=1,6
         xx=0.0d0
         yy=0.0d0
         zz=0.0d0
         do j=1,3
           xx=xx+bmat(j,1)*ixmat(j,i)
           yy=yy+bmat(j,2)*ixmat(j,i)
           zz=zz+bmat(j,3)*ixmat(j,i)
         end do
         bv(1,i)=xx
         bv(2,i)=yy
         bv(3,i)=zz
      end do

      do i=1,6
         wrk(1,i)=bv(1,i)*bv(1,i)
         wrk(2,i)=bv(1,i)*bv(2,i)
         wrk(3,i)=bv(1,i)*bv(3,i)
         wrk(4,i)=bv(2,i)*bv(2,i)
         wrk(5,i)=bv(2,i)*bv(3,i)
         wrk(6,i)=bv(3,i)*bv(3,i)
         wts(i)=0.0d0
      end do

*     *** scal=(2*pi/L)**2 ***
      scal = lattice_omega()**(1.0d0/3.0d0)
      scal = 8.0*datan(1.0d0)/scal
      scal = scal*scal
      wts(1)=scal
      wts(4)=scal
      wts(6)=scal
      call dgels('N',6,6,1,wrk,6,wts,6,tmp1,tmp_len,info)
      if (info.ne.0) then
        write(*,*)"Illegal argument in call to dgels"
        call flush(6)
      end if
      rank=0
      do i=1,6
         if (dabs(wts(i)).gt.1.e-6) then
           rank=rank+1
           wrk(1,rank)=bv(1,i)
           wrk(2,rank)=bv(2,i)
           wrk(3,rank)=bv(3,i)
           wrk(4,rank)=wts(i)
         end if
      end do
      do i=1,rank
         bv(1,i)=wrk(1,i)
         bv(2,i)=wrk(2,i)
         bv(3,i)=wrk(3,i)
         wts(i)=wrk(4,i)
      end do

      nx=0.0d0
      ny=0.0d0
      nz=0.0d0
      pcharge = 0.0d0
      do i=1,ion_nion()
        j=ion_katm(i)
        nx=nx+psp_zv(j)*ion_rion(1,i)
        ny=ny+psp_zv(j)*ion_rion(2,i)
        nz=nz+psp_zv(j)*ion_rion(3,i)
        pcharge = pcharge + psp_zv(j)
      end do   

      dx(1)=0.0d0
      dx(2)=0.0d0
      dy(1)=0.0d0
      dy(2)=0.0d0
      dz(1)=0.0d0
      dz(2)=0.0d0

      !if (oprint) then
      !   write(*,1260) 
      !   write(*,1261) rank
      !   do i=1,rank
      !      write(*,1262) i,bv(1,i),bv(2,i),bv(3,i),wts(i)
      !   end do
      !end if
      
      do ms=1,ispin
         do i=1,rank
           b(1) = bv(1,i)
           b(2) = bv(2,i)
           b(3) = bv(3,i)
           call silvestrelli_overlap(
     >                   b,ms,ne,neq,
     >                   dbl_mb(psi_r(1)),
     >                   dbl_mb(psi_r2(1)),
     >                   dcpl_mb(X(1,i)))
           call Dneall_w_eigenvalues(ms,dcpl_mb(X(1,i)),
     >                              dcpl_mb(Xeig(1,i)))
         end do
         do i=1,ne(ms)
            shift=(i-1)*ne(ms)+(i-1)
            xx=0.0d0
            yy=0.0d0
            zz=0.0d0
            do j=1,rank
 
               !*** really just want complex eigenvalues of X here ***
               !arg=dcpl_mb(X(1,j)+shift)
               arg=dcpl_mb(Xeig(1,j)+(i-1))
               arg= -wts(j)*datan2(dimag(arg),dble(arg))
               xx=xx+bv(1,j)*arg/scal
               yy=yy+bv(2,j)*arg/scal
               zz=zz+bv(3,j)*arg/scal
            end do
            dx(ms)=dx(ms)+xx
            dy(ms)=dy(ms)+yy
            dz(ms)=dz(ms)+zz
         end do
      end do !* ms *

cccccccccccccccccccccccccccccccccccccccccccccc
c  Molecular dipoles from Resta's theory!
ccccccccccccccccccccccccccccccccccccccccccccc

      tx=nx-dx(1)-dx(ispin)
      ty=ny-dy(1)-dy(ispin)
      tz=nz-dz(1)-dz(ispin)
      xx = dsqrt(tx*tx + ty*ty + tz*tz)
      dipole(1) = tx
      dipole(2) = ty
      dipole(3) = tz

      if (oprint) then
         write(luout,1771)
         write(luout,1772) 'spin up   ',
     >                     dx(1)/dble(ne(1)),
     >                     dy(1)/dble(ne(1)),
     >                     dz(1)/dble(ne(1))
         if (ne(ispin).ne.0)
     >      write(luout,1772) 'spin down ',
     >                        dx(ispin)/dble(ne(ispin)),
     >                        dy(ispin)/dble(ne(ispin)),
     >                        dz(ispin)/dble(ne(ispin))
         write(luout,1772) 'electronic',
     >                      (dx(1)+dx(ispin))/dble(ne(1)+ne(ispin)),
     >                      (dy(1)+dy(ispin))/dble(ne(1)+ne(ispin)),
     >                      (dz(1)+dz(ispin))/dble(ne(1)+ne(ispin))
         write(luout,1772) 'ionic     ',
     >                      nx/pcharge,
     >                      ny/pcharge,
     >                      nz/pcharge
         write(luout,1778)
         write(luout,1774) tx,ty,tz
         write(luout,1775) xx,xx*autoDebye
      end if

      value = .true.
      do j=6,1,-1
         value = value.and.BA_pop_stack(Xeig(2,j))
         value = value.and.Dneall_w_pop_stack(X(1,j))
      end do
      value = value.and.BA_pop_stack(psi_r2(2))
      value = value.and.BA_pop_stack(psi_r(2))
      if (.not. value) 
     >   call errquit('dipole_resta:pop stack',2,0)


      return 
*:::::::::::::::::::::::::::  format  :::::::::::::::::::::::::::::::::

 1771 FORMAT(//'== Center of Charge =='/)
 1772 FORMAT(A10,'  (',F10.4,',',F10.4,',',F10.4,' )')
 1773 FORMAT(//'== Wannier Crystal Dipole =='/)
 1774 FORMAT('mu   =  (',F10.4,',',F10.4,',',F10.4,' ) au')
 1775 FORMAT('|mu| = ',F10.4,' au,   ',F10.4,' Debye')
 1776 FORMAT(/"ELECTRONIC DIPOLES")
 1777 FORMAT("DX =",F11.5," DY= ",F11.5," DZ= ",F11.5)
 1778 FORMAT(//'== Crystal Dipole (Resta) =='/)
 1780 FORMAT("NUCLEAR DIPOLES")
 1785 FORMAT("TOTAL DIPOLES")
 
      end




