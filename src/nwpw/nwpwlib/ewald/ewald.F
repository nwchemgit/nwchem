*
* $Id$
*
      integer function ewald_ncut()
      implicit none 
#include "ewald.fh"
      ewald_ncut = ncut
      return
      end 
ccccccccccccccccccccccccccccccccccccc
      integer function ewald_nida()
      implicit none 
#include "ewald.fh"
      ewald_nida = nida
      return
      end
ccccccccccccccccccccccccccccccccccccc
      integer function ewald_npack()
      implicit none 
#include "ewald.fh"
      ewald_npack = enpack
      return
      end 
ccccccccccccccccccccccccccccccccccccc
      real*8 function ewald_zv(i)
      implicit none
#include "bafdecls.fh"
#include "ewald.fh"
      integer i,ix
      ix=zv(1)+i-1
      ewald_zv=dbl_mb(ix)      
      return
      end
cccccccccccccccccccccccccccccccccccc
      real*8 function ewald_rcut()
      implicit none 
#include "ewald.fh"
      ewald_rcut = rcut
      return
      end 
ccccccccccccccccccccccccccccccccccccccccc
      integer function ewald_nshl3d()
      implicit none 
#include "ewald.fh"
      ewald_nshl3d = nshl3d
      return
      end 
ccccccccccccccccccccccccccccccccccccc
      integer function ewald_rcell_ptr()
      implicit none
#include "ewald.fh"
      ewald_rcell_ptr = rcell(1)
      return
      end
ccccccccccccccccccccccccccccccccccccc
      real*8 function ewald_mandelung()
      implicit none 
#include "ewald.fh"
      ewald_mandelung = alpha
      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function ewald_grid_nx()
      implicit none
#include "ewald.fh"
      ewald_grid_nx=enx
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function ewald_grid_ny()
      implicit none
#include "ewald.fh"
      ewald_grid_ny=eny
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function ewald_grid_nz()
      implicit none
#include "ewald.fh"
      ewald_grid_nz=enz
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      intialization routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mandelung_set(alpha)
      implicit none
      real*8 alpha
      
      integer N
      parameter (N=40)

*     **** local variables ****
      integer n1,n2,n3,nterm
      real*8  rc,rs,epsilon,pi
      real*8  a1,a2,a3,g1,g2,g3,gg
      real*8  unita(3,3),unitg(3,3)
      real*8  alpha1,alpha2,sum,ea,omega

*     **** external functions ****
      real*8   util_erfc
      external util_erfc
      real*8   lattice_unita,lattice_unitg,lattice_omega
      external lattice_unita,lattice_unitg,lattice_omega
      
      pi = 4.0d0*datan(1.0d0)

*     ***** set lattice parameters *****
      omega = lattice_omega()
      do n1=1,3
      do n2=1,3
         unita(n1,n2) = lattice_unita(n1,n2)
         unitg(n1,n2) = lattice_unitg(n1,n2)
      end do
      end do

*     ***** set cutoff radii ****
      rs      = (3.0d0*omega/(4.0d0*pi))**(1.0d0/3.0d0)
      rc      = rs
      epsilon = 1.0d0/rc

*     **** calculate alpha1 *****
      sum = 0.0d0
      do n1=(-N+1),(N-1) 
      do n2=(-N+1),(N-1) 
      do n3=(-N+1),(N-1) 
         nterm=iabs(n1)+iabs(n2)+iabs(n3)
         if (nterm.ne.0) then
            a1 = n1*unita(1,1)
     >         + n2*unita(1,2)
     >         + n3*unita(1,3)

            a2 = n1*unita(2,1)
     >         + n2*unita(2,2)
     >         + n3*unita(2,3)

            a3 = n1*unita(3,1)
     >         + n2*unita(3,2)
     >         + n3*unita(3,3)

            ea = dsqrt(a1*a1 + a2*a2 + a3*a3)

            sum = sum + util_erfc(epsilon*ea)/ea

         end if
      end do
      end do
      end do
      alpha1 = sum

		
*     **** calculate alpha2 *****
      sum = 0.0d0
      do n1=(-N+1),(N-1) 
      do n2=(-N+1),(N-1) 
CDIR$ NOVECTOR
      do n3=(-N+1),(N-1) 
         nterm=iabs(n1)+iabs(n2)+iabs(n3)
         if (nterm.ne.0) then
            g1 = n1*unitg(1,1)
     >         + n2*unitg(1,2)
     >         + n3*unitg(1,3)

            g2 = n1*unitg(2,1)
     >         + n2*unitg(2,2)
     >         + n3*unitg(2,3)

            g3 = n1*unitg(3,1)
     >         + n2*unitg(3,2)
     >         + n3*unitg(3,3)

            gg  = g1*g1 + g2*g2 + g3*g3
            sum = sum +  (4.0d0*pi/gg)* exp(-gg*rc*rc/4.0d0)
	    
         end if
      end do
      end do
      end do
      alpha2 = sum/omega

      sum = alpha1 + alpha2 
     >    - pi*rc*rc/omega - 2.0d0*epsilon/dsqrt(pi)

      alpha = -sum*rs

      return 
      end
cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
      subroutine ewald_init()
      implicit none
#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"
*     **** arguments **********
*     **** local variables ****
      real*8 eps
      parameter(eps=1.0d-12)
      integer i,j,k,l
      integer k1,k2,k3,itmp,i_indx,j_indx,k_indx
      integer enxh,enyh,enzh,enpack0
      integer dutask,indx,xtra
      integer ntyp,ccode
      real*8  pi,pi4,gg,w,ggcut
      real*8  rs,qi
      real*8  zz,z,term
      logical value
      real*8 g1,g2,g3,gg1,gg2,gg3
      real*8 gf1,gf2,gf3,ggf1,ggf2,ggf3
      real*8 unitg(3,3),unita(3,3)
      real*8 unitf(3,3)

*     **** external functions ****
      integer  ewald_strfact_i_indx
      integer  ewald_strfact_j_indx
      integer  ewald_strfact_k_indx
      integer  control_code,ion_nkatm
      integer  control_ncut
      real*8   control_rcut,control_ecut
      integer  ion_nion,ion_katm,ion_nion_qm
      real*8   lattice_omega,lattice_unita
      real*8   lattice_unitg,lattice_unitg_frozen
      real*8 psp_zv,cpsp_zv,paw_basis_ion_charge,mmq_zv
      real*8 paw_basis_core_charge,dsum
      integer control_ewald_ngrid
  
      external ewald_strfact_i_indx
      external ewald_strfact_j_indx
      external ewald_strfact_k_indx
      external control_code,ion_nkatm
      external control_ncut
      external control_rcut,control_ecut
      external ion_nion,ion_katm,ion_nion_qm
      external lattice_omega,lattice_unita                                
      external lattice_unitg,lattice_unitg_frozen
      external psp_zv,cpsp_zv,paw_basis_ion_charge,mmq_zv
      external paw_basis_core_charge,dsum
      external control_ewald_ngrid

      do j=1,3
      do i=1,3
          unita(i,j) = lattice_unita(i,j)
          unitg(i,j) = lattice_unitg(i,j)
          unitf(i,j) = lattice_unitg_frozen(i,j)
      end do
      end do
      call Parallel_np(tnp)
      call Parallel_taskid(tid)
      enx=control_ewald_ngrid(1)
      eny=control_ewald_ngrid(2)
      enz=control_ewald_ngrid(3)
      enxh=enx/2
      enyh=eny/2
      enzh=enz/2

*     **** determine ggcut ****
      g1 = unitf(1,1)*(enxh)
      g2 = unitf(2,1)*(enxh)
      g3 = unitf(3,1)*(enxh)
      gg1 = g1*g1 + g2*g2 + g3*g3

      g1 = unitf(1,2)*(enyh)
      g2 = unitf(2,2)*(enyh)
      g3 = unitf(3,2)*(enyh)
      gg2 = g1*g1 + g2*g2 + g3*g3

      g1 = unitf(1,3)*(enzh)
      g2 = unitf(2,3)*(enzh)
      g3 = unitf(3,3)*(enzh)
      gg3 = g1*g1 + g2*g2 + g3*g3

      ggcut = gg1
      if (gg2.lt.ggcut) ggcut = gg2
      if (gg3.lt.ggcut) ggcut = gg3
      if ((2.0d0*control_ecut()).lt.ggcut) ggcut=2.0d0*control_ecut()

c     **** determine enpack and nida ****
      dutask=0
      enpack = 0
      nida   = 0
      k1 = 0
      k2 = 0
      k3 = 0
      g1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
      g2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
      g3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
      gg=g1*g1+g2*g2+g3*g3
      if ((gg-ggcut).lt.-eps) then
         if (dutask.eq.tid) then
            enpack = enpack + 1
            nida   = nida + 1
         end if 
         dutask = mod(dutask+1,tnp)
      end if
      do k=1,(enzh-1)
         k1 = 0
         k2 = 0
         k3 = k
         g1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         g2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         g3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=g1*g1+g2*g2+g3*g3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) enpack = enpack + 1
            dutask = mod(dutask+1,tnp)
         end if
      end do
      do k=(-enzh+1),(enzh-1)
      do j=1,(enyh-1)
         k1 = 0
         k2 = j
         k3 = k
         g1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         g2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         g3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=g1*g1+g2*g2+g3*g3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) enpack = enpack + 1
            dutask = mod(dutask+1,tnp)
         end if
      end do 
      end do 
      do k=(-enzh+1),(enzh-1)
      do j=(-enyh+1),(enyh-1)
      do i=1,(enxh-1)
         k1=i
         k2=j
         k3=k
         g1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         g2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         g3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=g1*g1+g2*g2+g3*g3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) enpack = enpack + 1
            dutask = mod(dutask+1,tnp)
         end if
      end do
      end do
      end do

c     **** allocate memory for common block ****
      notzero_enpack = (enpack.gt.0)  !** needed for very large numbers cpu
      value = .true.
      if (notzero_enpack) then 
         value=BA_alloc_get(mt_dbl,3*enpack,'eG',eg_hndl,eG(1))
         eG(2)=eG(1)+enpack
         eG(3)=eG(2)+enpack
         value=value.and.BA_alloc_get(mt_dbl,enpack,'vg',vg(2),vg(1))
         value=value.and.BA_alloc_get(mt_dbl,enpack,'vcx',vcx(2),vcx(1))
      end if
      ntyp=ion_nkatm()      
      value=value.and.BA_alloc_get(mt_dbl,ntyp,'zv',zv(2),zv(1))
      if (.not.value) call errquit('ewald_init:out of heap memory',0,0)

      call ewald_strfact_init()


*     **** determine eG ****
      if (notzero_enpack) then
      i_indx = ewald_strfact_i_indx()
      j_indx = ewald_strfact_j_indx()
      k_indx = ewald_strfact_k_indx()
      call dcopy(3*enpack,0.0d0,0,dbl_mb(eG(1)),1)
      dutask=0
      enpack0 = 0
      k1 = 0
      k2 = 0
      k3 = 0
      g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
      g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
      g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)

      gf1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
      gf2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
      gf3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
      gg=gf1*gf1+gf2*gf2+gf3*gf3
      if ((gg-ggcut).lt.-eps) then
         if (dutask.eq.tid) then
            dbl_mb(eG(1)+enpack0) = g1
            dbl_mb(eG(2)+enpack0) = g2
            dbl_mb(eG(3)+enpack0) = g3
            i=k1
            j=k2
            k=k3
            if (i .lt. 0) i = i + enx
            if (j .lt. 0) j = j + eny
            if (k .lt. 0) k = k + enz
            int_mb(i_indx+enpack0) = i+1
            int_mb(j_indx+enpack0) = j+1
            int_mb(k_indx+enpack0) = k+1
            enpack0 = enpack0 + 1
         endif
         dutask = mod(dutask+1,tnp)
      end if

      k1 = 0
      k2 = 0
      do k3=1,(enzh-1)
         g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
         g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
         g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)

         gf1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         gf2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         gf3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=gf1*gf1+gf2*gf2+gf3*gf3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) then
               dbl_mb(eG(1)+enpack0) = g1
               dbl_mb(eG(2)+enpack0) = g2
               dbl_mb(eG(3)+enpack0) = g3
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + enx
               if (j .lt. 0) j = j + eny
               if (k .lt. 0) k = k + enz
               int_mb(i_indx+enpack0) = i+1
               int_mb(j_indx+enpack0) = j+1
               int_mb(k_indx+enpack0) = k+1
               enpack0 = enpack0 + 1
            end if
            dutask = mod(dutask+1,tnp)
         end if
      end do
      k1 = 0
      do k3=(-enzh+1),(enzh-1)
      do k2=1,(enyh-1)
         g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
         g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
         g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)

         gf1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         gf2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         gf3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=gf1*gf1+gf2*gf2+gf3*gf3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) then
               dbl_mb(eG(1)+enpack0) = g1
               dbl_mb(eG(2)+enpack0) = g2
               dbl_mb(eG(3)+enpack0) = g3
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + enx
               if (j .lt. 0) j = j + eny
               if (k .lt. 0) k = k + enz
               int_mb(i_indx+enpack0) = i+1
               int_mb(j_indx+enpack0) = j+1
               int_mb(k_indx+enpack0) = k+1
               enpack0 = enpack0 + 1
            end if
            dutask = mod(dutask+1,tnp)
         end if
      end do 
      end do 
      do k3=(-enzh+1),(enzh-1)
      do k2=(-enyh+1),(enyh-1)
      do k1=1,(enxh-1)
         g1=k1*unitg(1,1)+k2*unitg(1,2)+k3*unitg(1,3)
         g2=k1*unitg(2,1)+k2*unitg(2,2)+k3*unitg(2,3)
         g3=k1*unitg(3,1)+k2*unitg(3,2)+k3*unitg(3,3)

         gf1=k1*unitf(1,1)+k2*unitf(1,2)+k3*unitf(1,3)
         gf2=k1*unitf(2,1)+k2*unitf(2,2)+k3*unitf(2,3)
         gf3=k1*unitf(3,1)+k2*unitf(3,2)+k3*unitf(3,3)
         gg=gf1*gf1+gf2*gf2+gf3*gf3
         if ((gg-ggcut).lt.-eps) then
            if (dutask.eq.tid) then
               dbl_mb(eG(1)+enpack0) = g1
               dbl_mb(eG(2)+enpack0) = g2
               dbl_mb(eG(3)+enpack0) = g3
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + enx
               if (j .lt. 0) j = j + eny
               if (k .lt. 0) k = k + enz
               int_mb(i_indx+enpack0) = i+1
               int_mb(j_indx+enpack0) = j+1
               int_mb(k_indx+enpack0) = k+1
               enpack0 = enpack0 + 1
            end if
            dutask = mod(dutask+1,tnp)
         end if
      end do
      end do
      end do

      call dcopy(enpack,0.0d0,0,dbl_mb(vg(1)),1)
      call dcopy(enpack,0.0d0,0,dbl_mb(vcx(1)),1)
      end if

*     ***** find w *****
      pi  = 4.0d0*datan(1.0d0)
      pi4 = 4.0d0*pi
      ncut = control_ncut()
      rcut = control_rcut()
      if (ncut.le.0)     ncut=1
      if (rcut.le.0.0d0) then
         rs = unita(1,1)**2 + unita(2,1)**2 + unita(3,1)**2
         rs = dsqrt(rs)
         rcut=rs/pi

         rs = unita(1,2)**2 + unita(2,2)**2 + unita(3,2)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w

         rs = unita(1,3)**2 + unita(2,3)**2 + unita(3,3)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w
      end if
      w = 0.25d0*rcut*rcut

      if (notzero_enpack) then
         dbl_mb(vcx(1)) = 0.0d0
         dbl_mb(vg(1))  = 0.0d0
         do k=1+nida,enpack
           g1 = dbl_mb(eG(1)+k-1)
           g2 = dbl_mb(eG(2)+k-1)
           g3 = dbl_mb(eG(3)+k-1)
           gg=g1*g1+g2*g2+g3*g3
           term=pi4/gg
           dbl_mb(vcx(1)+k-1)=term
           dbl_mb(vg(1)+k-1) =term*dexp(-w*gg)
         end do
      end if

ccccccccccc set Mandelung constant       

      call Mandelung_set(alpha)

ccccccccccc set ion charges
      ccode=control_code()
      if (ccode.eq.12) then
         do i=1,ntyp
            dbl_mb(zv(1)+i-1)=mmq_zv(i)
         end do
      else if ((ccode.eq.5).or.(ccode.eq.13).or.(ccode.eq.14)) then
          do i=1,ntyp
             dbl_mb(zv(1)+i-1)=cpsp_zv(i)
          end do
      else
        if ((ccode.eq.6).or.(ccode.eq.7).or.(ccode.eq.8)) then
          do i=1,ntyp
              dbl_mb(zv(1)+i-1)=paw_basis_ion_charge(i)-
     >                 paw_basis_core_charge(i)
          end do
        else
          do i=1,ntyp
             dbl_mb(zv(1)+i-1)=psp_zv(i)
          end do
        end if
      end if

cccccccccc  ewald summation
      rs = (3.0d0*lattice_omega()/pi4)**(1.0d0/3.0d0)

      zz = 0.0d0
      z  = 0.0d0
      do i=1,ion_nion()
         qi=dbl_mb(zv(1)+ion_katm(i)-1)
         zz = zz + qi**2
         z  = z  + qi
      end do

*     ***** t_dsum of vg *****
      if (notzero_enpack) then
         cewald = 2.0d0*dsum(enpack,dbl_mb(vg(1)),1)   
      else
         cewald=0.0d0
      end if
      !if (nida.eq.1) cewald = cewald - dbl_mb(vg(1))  !*** zero because vg(g=0) = 0
      if (tnp.gt.1) call Parallel_SumAll(cewald)

      cewald = -0.5d0*zz*(alpha/rs + cewald/lattice_omega())
     >         -0.5d0*(z*z-zz)*rcut*rcut*pi/lattice_omega()


      zz = 0.0d0
      z  = 0.0d0
      do i=1,ion_nion_qm()
         qi=dbl_mb(zv(1)+ion_katm(i)-1)
         zz = zz + qi**2
         z  = z  + qi
      end do

*     ***** t_dsum of vg *****
      if (notzero_enpack) then
         cewald_qm = 2.0d0*dsum(enpack,dbl_mb(vg(1)),1)   
      else
         cewald_qm = 0.0d0
      end if
      !if (nida.eq.1) cewald_qm = cewald_qm - dbl_mb(vg(1))  !*** zero because vg(g=0) = 0
      if (tnp.gt.1) call Parallel_SumAll(cewald_qm)

      cewald_qm = -0.5d0*zz*(alpha/rs + cewald_qm/lattice_omega())
     >         -0.5d0*(z*z-zz)*rcut*rcut*pi/lattice_omega()

      zz = 0.0d0
      z  = 0.0d0
      do i=ion_nion_qm(1)+1,ion_nion()
         qi=dbl_mb(zv(1)+ion_katm(i)-1)
         zz = zz + qi**2
         z  = z  + qi
      end do

*     ***** t_dsum of vg *****
      if (notzero_enpack) then
         cewald_mm = 2.0d0*dsum(enpack,dbl_mb(vg(1)),1)   
      else
         cewald_mm = 0.0d0
      end if
      !if (nida.eq.1) cewald_mm = cewald_mm - dbl_mb(vg(1))  !*** zero because vg(g=0) = 0
      if (tnp.gt.1) call Parallel_SumAll(cewald_mm)

      cewald_mm = -0.5d0*zz*(alpha/rs + cewald_mm/lattice_omega())
     >         -0.5d0*(z*z-zz)*rcut*rcut*pi/lattice_omega()


cccccccccc allocate and fill rcell
      nshl3d=(2*ncut+1)**3
      value=BA_alloc_get(mt_dbl,(3*nshl3d),'rcell',rcell(2),rcell(1))
      if (.not. value)
     > call errquit('ewald_init:out of heap memory',0,MA_ERR)


*     **** get lattice vectors in real space - define rcell(1) = origin****
      l=1
      dbl_mb(rcell(1)+(l-1))          = 0.0d0
      dbl_mb(rcell(1)+(l-1)+nshl3d)   = 0.0d0
      dbl_mb(rcell(1)+(l-1)+2*nshl3d) = 0.0d0
      do k=-ncut,ncut
      do j=-ncut,ncut
      do i=-ncut,ncut
      if (.not.((i.eq.0).and.(j.eq.0).and.(k.eq.0))) then
        l = l+1
        dbl_mb(rcell(1)+(l-1) )
     >    =i*unita(1,1)+j*unita(1,2)+k*unita(1,3)
        dbl_mb(rcell(1)+(l-1)+nshl3d)
     >    =i*unita(2,1)+j*unita(2,2)+k*unita(2,3)
        dbl_mb(rcell(1)+(l-1)+2*nshl3d)
     >    =i*unita(3,1)+j*unita(3,2)+k*unita(3,3)
      end if
      end do
      end do
      end do

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ewald_end()
      implicit none
#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"
      logical value
      call ewald_strfact_end()
      value = .true.
      if (notzero_enpack) then
         value=value.and.BA_free_heap(eg_hndl)
         value=value.and.BA_free_heap(vcx(2))
         value=value.and.BA_free_heap(vg(2))
      end if
      value=value.and.BA_free_heap(rcell(2))
      value=value.and.BA_free_heap(zv(2))
      if (.not.value) then
        call errquit("ewald end can free heap",0,MA_ERR)
      end if
      return 
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*     ***********************************
*     *                                 *
*     *          ewald_ct_Sqr           *
*     *                                 *
*     ***********************************
      subroutine ewald_ct_Sqr(ng,A,C)
      implicit none
      integer    ng
      complex*16 A(*)
      real*8     C(*)

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,ng
        C(i) = dble(A(i))**2 + dimag(A(i))**2
      end do
!$OMP END DO
      return
      end

      subroutine ewald_tt_2dot(ng,A,B,esum)
      implicit none
      integer    ng
      real*8  A(*)
      real*8  B(*)
      real*8  esum

c     **** local variables ****
      integer i
!OMP MASTER
      esum = 0.0d0
!OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:esum)
      do i=1,ng
        esum = esum + 2.0d0*A(i)*B(i)
      end do
!$OMP END DO
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*     ***********************************
*     *			  		*
*     *		ewald_e			*
*     *				       	*
*     ***********************************
      real*8 function ewald_e()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"


*     **** Error function parameters ****
      logical fast_erf
      real*8 xerf,yerf,yerfc
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)


*     **** local variables ****
      integer dutask,threadid,nthreads,ttid
      integer i,j,ii,l,nion
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  energy,etmp,energyG,ttcr,ttci,zii,etmp0
      complex*16 ttcz
*     **** temporary workspace variables ****
      integer exi(2),s(2),tmp3(2)
      logical value

      common /eiond_ejtmp/ etmp

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      real*8   ddot
      external ion_nion,ion_katm
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      external ddot
      logical  control_fast_erf
      external control_fast_erf
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(21)

*     **** allocate temp workspace ****
      fast_erf = control_fast_erf()
      nion = ion_nion()
      threadid = Parallel_threadid()
      nthreads = Parallel_nthreads()
      ttid = threadid + nthreads*tid
      
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.
     >           BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if (.not. value) 
     >     call errquit('ewald_e:out of stack memory',0,MA_ERR)
     
*        **** get the structure factor ****
         !call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         call Parallel_shared_vector_zero(.true.,(2*enpack),
     >                                    dcpl_mb(s(1)))
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            call DAXPY_OMP(2*enpack,dbl_mb(zv(1)+ion_katm(ii)-1),
     >             dcpl_mb(exi(1)),1,
     >             dcpl_mb(s(1)),1)
         end do

*        **** calculate the ewald energy ****
         call ewald_ct_Sqr(enpack,dcpl_mb(s(1)),dbl_mb(tmp3(1)))
         !etmp=2.0d0*ddot(enpack,dbl_mb(tmp3(1)),1,dbl_mb(vg(1)),1)
         call ewald_tt_2dot(enpack,dbl_mb(tmp3(1)),dbl_mb(vg(1)),etmp)
      else
!$OMP MASTER
         etmp=0.0d0
!$OMP END MASTER
!$OMP BARRIER
      end if
      if (tnp.gt.1) call Parallel_SumAll(etmp)
      energy  = 0.5d0*etmp/lattice_omega() + cewald
!$OMP BARRIER

*     *** made parallel  ****
      dutask = 0
      etmp0 = 0.0d0
!$OMP MASTER
      etmp = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      do i=1,nion-1
      do j=i+1,nion
      if (dutask.eq.ttid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = dbl_mb(zv(1)+ion_katm(i)-1)*dbl_mb(zv(1)+ion_katm(j)-1)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           if ((r.gt.1.0d-6).and.(fast_erf)) then
              xerf=r/rcut
              yerf = (1.0d0
     >              + xerf*(c1 + xerf*(c2
     >              + xerf*(c3 + xerf*(c4
     >              + xerf*(c5 + xerf*c6))))))**4
              !yerfc = (1.0d0 - 1.0d0/yerf**4)
              yerfc = 1.0d0/yerf**4
              etmp0=etmp0+(zz*yerfc/r)
           else if (r.gt.1.0d-6) then
              w = (r/rcut)
              yerfc = util_erfc(w)  !*** needs to be faster
              etmp0=etmp0+(zz*yerfc/r)
           end if
        end do
      end if
      dutask=mod(dutask+1,tnp*nthreads)
      end do
      end do
!$OMP CRITICAL (ewaldp2)
      etmp = etmp + etmp0
!$OMP END CRITICAL (ewaldp2)
!$OMP BARRIER
      if (tnp.gt.1) call Parallel_SumAll(etmp)
      energy = energy + etmp
      

*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_e:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      ewald_e = energy
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f			*
*     *		       			*
*     ***********************************

      subroutine ewald_f(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)

*     **** Error function parameters ****
      logical fast_erf
      real*8 xerf,yerf,yerfc
      real*8 c1,c2,c3,c4,c5,c6
      parameter (c1=0.07052307840d0,c2=0.04228201230d0)
      parameter (c3=0.00927052720d0)
      parameter (c4=0.00015201430d0,c5=0.00027656720d0)
      parameter (c6=0.00004306380d0)

*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii
      common /ewald_f_sw1sw2/ sw1,sw2,sw3

*     **** temporary workspace variables ****
      integer nion
      integer exi(2),s(2),tmp3(2)
      integer fx(2),fy(2),fz(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,ewald_nshl3d
      external ion_nion,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      logical  control_fast_erf
      external control_fast_erf

      call nwpw_timing_start(21)
      nion = ion_nion()
      fast_erf = control_fast_erf()


*     **** allocate temp workspace ****
      value = .true.
      if (notzero_enpack) then
         value = value.and.
     >           BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
      end if
      value = value.and.BA_push_get(mt_dbl,nion,'fx',fx(2),fx(1)) 
      value = value.and.BA_push_get(mt_dbl,nion,'fy',fy(2),fy(1)) 
      value = value.and.BA_push_get(mt_dbl,nion,'fz',fz(2),fz(1)) 
      if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)


      scal2 = 1.0d0/lattice_omega()
      !call dcopy(nion,0.0d0,0,dbl_mb(fx(1)),1)
      !call dcopy(nion,0.0d0,0,dbl_mb(fy(1)),1)
      !call dcopy(nion,0.0d0,0,dbl_mb(fz(1)),1)
      call Parallel_shared_vector_zero(.false.,nion,dbl_mb(fx(1)))
      call Parallel_shared_vector_zero(.false.,nion,dbl_mb(fy(1)))
      call Parallel_shared_vector_zero(.true., nion,dbl_mb(fz(1)))
     
*     **** get the structure factor ****
      if (notzero_enpack) then
         !call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         call Parallel_shared_vector_zero(.true.,2*enpack,dcpl_mb(s(1)))
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            call daxpy_omp(2*enpack,zii,dcpl_mb(exi(1)),1,
     >                     dcpl_mb(s(1)),1)
         end do
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))

c         do i=1,enpack
c            dbl_mb(tmp3(1)+i-1) 
c     >              = ( dble(dcpl_mb(exi(1)+i-1))
c     >                *dimag(dcpl_mb(s(1)+i-1))
c     >              -  dimag(dcpl_mb(exi(1)+i-1))
c     >                 *dble(dcpl_mb(s(1)+i-1))
c     >                )*dbl_mb(vg(1)+i-1)
c         end do
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

!$OMP MASTER
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fx(1)+ii-1) = dbl_mb(fx(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fy(1)+ii-1) = dbl_mb(fy(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fz(1)+ii-1) = dbl_mb(fz(1)+ii-1) 
     >                         +  sum*zii*scal2
!$OMP END MASTER
!$OMP BARRIER
         end do
      end if
    
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
       if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
!$OMP MASTER
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sw1,sw2,sw3)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           if ((r.gt.1.0d-6).and.(fast_erf)) then
              w = r/rcut
              xerf=w
              yerf = (1.0d0
     >              + xerf*(c1 + xerf*(c2
     >              + xerf*(c3 + xerf*(c4
     >              + xerf*(c5 + xerf*c6))))))**4
              !yerf = (1.0d0 - 1.0d0/yerf**4)
              yerfc = 1.0d0/yerf**4
              f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
              sw1=sw1+(x*f)
              sw2=sw2+(y*f)
              sw3=sw3+(z*f)
           else if (r.gt.1.0d-6) then
              w = r/rcut
              yerfc = util_erfc(w)
              f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
              sw1=sw1+(x*f)
              sw2=sw2+(y*f)
              sw3=sw3+(z*f)
           end if
        end do
!$OMP END DO

!$OMP MASTER
        dbl_mb(fx(1)+i-1) = dbl_mb(fx(1)+i-1) + sw1
        dbl_mb(fy(1)+i-1) = dbl_mb(fy(1)+i-1) + sw2
        dbl_mb(fz(1)+i-1) = dbl_mb(fz(1)+i-1) + sw3

        dbl_mb(fx(1)+j-1) = dbl_mb(fx(1)+j-1) - sw1
        dbl_mb(fy(1)+j-1) = dbl_mb(fy(1)+j-1) - sw2
        dbl_mb(fz(1)+j-1) = dbl_mb(fz(1)+j-1) - sw3
!$OMP END MASTER
!$OMP BARRIER

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do
      if (tnp.gt.1) then
          call Parallel_Vector_SumAll(nion,dbl_mb(fx(1)))
          call Parallel_Vector_SumAll(nion,dbl_mb(fy(1)))
          call Parallel_Vector_SumAll(nion,dbl_mb(fz(1)))
      end if

!$OMP DO
      do i=1,nion
         fion(1,i) = fion(1,i) + dbl_mb(fx(1)+i-1)
         fion(2,i) = fion(2,i) + dbl_mb(fy(1)+i-1)
         fion(3,i) = fion(3,i) + dbl_mb(fz(1)+i-1)
      end do
!$OMP END DO
      
*     **** deallocate temp workspace ****
      value =           BA_pop_stack(fz(2))
      value = value.and.BA_pop_stack(fy(2))
      value = value.and.BA_pop_stack(fx(2))
      if (notzero_enpack) then
         value = value.and.BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
      end if
      if (.not. value) 
     >  call errquit('ewald_f:popping stack memory',0,MA_ERR)

      call nwpw_timing_end(21)
      return
      end

      subroutine ewald_f_tmp3(n,e,s,v,t)
      implicit none
      integer n
      complex*16 e(n),s(n)
      real*8 v(n),t(n)

      integer i

!$OMP DO
      do i=1,n
        t(i) = v(i)*(dble(e(i))*dimag(s(i)) - dimag(e(i))*dble(s(i)))
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *					*
*     *		ewald_stress		*
*     *	        			*
*     ***********************************

      subroutine ewald_stress(stress)
      implicit none
      real*8  stress(3,3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

      integer N
      parameter (N=40)

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)

*     **** local variables ****
      logical value
      integer i,ii,j,l,nion
      integer n1,n2,n3
      integer u,v,s,dutask
      real*8 pi,fourpi,scal
      real*8 zz,z
      real*8 Cus(3,3),hm(3,3),energy,sum,ss,rs
      real*8 ea,ax,ay,az,epsilon
      real*8 dx,dy,dz,w,ar,ai
      real*8 unita(3,3),unitg(3,3)
      complex*16 cz
      integer H(2),F(2),tmp1(2),tmp2(2),exi(2),strf(2)

*     **** external functions ****
      integer  ion_katm,ion_nion
      real*8   ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      real*8   util_erfc,ion_rion,ddot
      external ion_katm,ion_nion
      external ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      external util_erfc,ion_rion,ddot

      call nwpw_timing_start(21)
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)
      nion    = ion_nion()
*     *** define hm,unita,unitg ****
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
         unitg(u,v) = lattice_unitg(u,v)
         unita(u,v) = lattice_unita(u,v)
      end do
      end do
 
      zz = 0.0d0
      z  = 0.0d0
      do i=1,nion
         zz = zz + dbl_mb(zv(1)+ion_katm(i)-1)**2
         z  = z  + dbl_mb(zv(1)+ion_katm(i)-1)
      end do

*     **** Miscellaneous contributions - stress from cewald term ****
      do v=1,3
      do u=1,3
         stress(u,v) = 0.5d0*z*z*pi*rcut*rcut/lattice_omega()
     >               *hm(u,v)
      end do
      end do


*     **** G-space contributions ****

*     **** get the structure factor ****
      if (notzero_enpack) then
         value=            BA_push_get(mt_dbl,enpack,'H',H(2),H(1))
         value=value.and.BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1))
         value=value.and.
     >         BA_push_get(mt_dcpl,enpack,'strf',strf(2),strf(1))
         if (.not. value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)


         call dcopy((2*enpack),0.0d0,0,dcpl_mb(strf(1)),1)
         do ii=1,ion_nion()
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            call daxpy(2*enpack,ewald_zv(ion_katm(ii)),
     >                 dcpl_mb(exi(1)),1,dcpl_mb(strf(1)),1)
         end do
         do i=1,enpack
             cz=dcpl_mb(strf(1)+i-1)
             ar=dble(cz)
             ai=dimag(cz)
             dbl_mb(H(1)+i-1)=ar*ar+ai*ai
         end do
         value =           BA_pop_stack(strf(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >    call errquit('ewald_stress:error popping stack',0,MA_ERR)

*        **** calculate the ewald energy ****
c         call dcopy(enpack,dbl_mb(vg(1)),1,dbl_mb(F(1)),1)
         F(1) = vg(1)
         energy=2.0d0*ddot(enpack,dbl_mb(F(1)),1,dbl_mb(H(1)),1)
      else
         energy=0.0d0
      end if
      if (tnp.gt.1) call Parallel_SumAll(energy)
      energy = -0.5d0*energy/lattice_omega()


      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) + energy*hm(u,v)
      end do
      end do
     
*     **** tmp2(G) = F(G)*H(G)/G**2 + F(G)*H(G)*rcut*rcut/4 ****
      if (notzero_enpack) then
         value=BA_push_get(mt_dbl,enpack,'tmp1',tmp1(2),tmp1(1))
         value=value.and.
     >         BA_push_get(mt_dbl,enpack,'tmp2',tmp2(2),tmp2(1))
         if (.not.value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)

         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(F(1)+i-1)*dbl_mb(H(1)+i-1)
         end do
         ss = 0.25d0*rcut*rcut
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         ss = 1.0d0/fourpi
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*dbl_mb(vcx(1)+i-1)
         end do       
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp2(1)+i-1)+dbl_mb(tmp1(1)+i-1)
         end do       
      end if

*     **** calculate Cus ****

      call dcopy(9,0.0d0,0,Cus,1)
      ss =  1.0d0/lattice_omega()
      if (notzero_enpack) then
         do u=1,3
         do s=u,3
            do i=1,enpack
               dbl_mb(tmp1(1)+i-1)=dbl_mb(eG(u)+i-1)*dbl_mb(eG(s)+i-1)
            end do
            sum=2.0d0*ddot(enpack,dbl_mb(tmp1(1)),1,dbl_mb(tmp2(1)),1)
            !if (tnp.gt.1) call Parallel_SumAll(sum)
            Cus(u,s) = ss*sum
         end do
         end do
c      do u=1,3
c      do s=u+1,3
c         Cus(s,u) = Cus(u,s)
c      end do
c      end do
c      do v=1,3
c      do u=1,3
c        do s=1,3
c           stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c        end do
c      end do
c      end do

         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(H(2))
         if (.not.value) 
     >    call errquit(
     >         'ewald_stress:error popping stack memory',0,MA_ERR)
      end if

*     **** R-space contributions ****

*     **** calculate alpha1 - stress from cewald term*****
c        call dcopy(9,0.0d0,0,Cus,1)
        rs      = (3.0d0*lattice_omega()/(4.0d0*pi))**(1.0d0/3.0d0)
        epsilon = 1.0d0/rcut
        sum = 0.0d0
        dutask = 0
        do n1=(-N+1),(N-1) 
        do n2=(-N+1),(N-1) 
        do n3=(-N+1),(N-1) 
           if (dutask.eq.tid) then
           if (.not.((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0))) then
              ax = n1*unita(1,1)
     >           + n2*unita(1,2)
     >           + n3*unita(1,3)

              ay = n1*unita(2,1)
     >           + n2*unita(2,2)
     >           + n3*unita(2,3)

              az = n1*unita(3,1)
     >           + n2*unita(3,2)
     >           + n3*unita(3,3)

              ea = dsqrt(ax*ax + ay*ay + az*az)
              w = ea*epsilon

              ss = util_erfc(w)/ea 
     >           + 2.0d0*epsilon/dsqrt(pi)*dexp(-w*w)
              ss = -(0.5d0*zz)*ss/(ea*ea)
              Cus(1,1) = Cus(1,1) + ss * ax*ax 
              Cus(1,2) = Cus(1,2) + ss * ax*ay 
              Cus(1,3) = Cus(1,3) + ss * ax*az 

              !Cus(2,1) = Cus(2,1) + ss * ay*ax 
              Cus(2,2) = Cus(2,2) + ss * ay*ay 
              Cus(2,3) = Cus(2,3) + ss * ay*az 

              !Cus(3,1) = Cus(3,1) + ss * az*ax 
              !Cus(3,2) = Cus(3,2) + ss * az*ay 
              Cus(3,3) = Cus(3,3) + ss * az*az 

           end if
           end if
           dutask=mod(dutask+1,tnp)
        end do
        end do
        end do


c        do v=1,3
c        do u=1,3
c          do s=1,3
c             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c          end do
c        end do
c        end do


*     **** calculate erfc contribution *****
c      call dcopy(9,0.0d0,0,Cus,1)
      epsilon = 1.0d0/rcut
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
        if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
        do l=1,nshl3d
           ax = dbl_mb(rcell(1)+(l-1))          + dx
           ay = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           az = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           ea = dsqrt(ax*ax+ay*ay+az*az)
           if (ea.gt.1.0d-6) then
              w = ea*epsilon

              ss = -util_erfc(w)/ea
     >           - 2.0d0*epsilon/dsqrt(pi)*exp(-w*w)
              ss = ss/(ea*ea)
              Cus(1,1) = Cus(1,1) + ss * ax*ax * zz
              Cus(1,2) = Cus(1,2) + ss * ax*ay * zz
              Cus(1,3) = Cus(1,3) + ss * ax*az * zz
              Cus(2,2) = Cus(2,2) + ss * ay*ay * zz
              Cus(2,3) = Cus(2,3) + ss * ay*az * zz
              Cus(3,3) = Cus(3,3) + ss * az*az * zz
           end if
        end do
        end if
        dutask=mod(dutask+1,tnp)
      end do
      end do

      if (tnp.gt.1) then
        call Parallel_Vector_SumAll(9,Cus(1,1))
      end if

      do u=1,3
        do s=u+1,3
           Cus(s,u) = Cus(u,s)
        end do
      end do

      do v=1,3
        do u=1,3
          do s=1,3
             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
          end do
        end do
      end do

      call nwpw_timing_end(21)
      return
      end


*     ***********************************
*     *			  		*
*     *		ewald_e_qm		*
*     *				       	*
*     ***********************************
      real*8 function ewald_e_qm()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     **** local variables ****
      integer dutask
      integer i,j,ii,l,nion
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  energy,etmp,energyG,ttcr,ttci,zii
      complex*16 ttcz
*     **** temporary workspace variables ****
      integer exi(2),s(2),tmp3(2)
      logical value

*     **** external functions ****
      integer  ion_nion_qm,ion_katm
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      real*8   ddot
      external ion_nion_qm,ion_katm
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      external ddot

      call nwpw_timing_start(21)

*     **** allocate temp workspace ****
      nion = ion_nion_qm()

      if (notzero_enpack) then      
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.
     >           BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if (.not. value) 
     >     call errquit('ewald_e:out of stack memory',0,MA_ERR)
     
*        **** get the structure factor ****
         call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            call daxpy(2*enpack,dbl_mb(zv(1)+ion_katm(ii)-1),
     >             dcpl_mb(exi(1)),1,
     >             dcpl_mb(s(1)),1)
         end do

*        **** calculate the ewald energy ****
         call ewald_ct_Sqr(enpack,dcpl_mb(s(1)),dbl_mb(tmp3(1)))
         energy=2.0d0*ddot(enpack,dbl_mb(tmp3(1)),1,dbl_mb(vg(1)),1)
      else
         energy=0.0d0
      end if
      if (tnp.gt.1) call Parallel_SumAll(energy)
      energy  = 0.5d0*energy/lattice_omega() + cewald_qm

      
*     *** made parallel  ****
      dutask = 0
      etmp = 0.0d0
      do i=1,nion-1
      do j=i+1,nion
      if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = dbl_mb(zv(1)+ion_katm(i)-1)*dbl_mb(zv(1)+ion_katm(j)-1)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = (r/rcut)
           yerfc = util_erfc(w)  !*** needs to be faster
           etmp=etmp+(zz*yerfc/r)
        end do
      end if
      dutask=mod(dutask+1,tnp)
      end do
      end do
      if (tnp.gt.1) call Parallel_SumAll(etmp)
      energy = energy + etmp
      

*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_e:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      ewald_e_qm = energy
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f_qm		*
*     *		       			*
*     ***********************************

      subroutine ewald_f_qm(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii

*     **** temporary workspace variables ****
      integer nion
      integer exi(2),s(2),tmp3(2)
      integer fx(2),fy(2),fz(2)
      logical value

*     **** external functions ****
      integer  ion_nion_qm,ion_katm,ewald_nshl3d
      external ion_nion_qm,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot

      call nwpw_timing_start(21)
      nion = ion_nion_qm()

*     **** allocate temp workspace ****
      value = .true.
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
      end if
      value = value.and.BA_push_get(mt_dbl,nion,'fx',fx(2),fx(1)) 
      value = value.and.BA_push_get(mt_dbl,nion,'fy',fy(2),fy(1)) 
      value = value.and.BA_push_get(mt_dbl,nion,'fz',fz(2),fz(1)) 
      if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)


      scal2 = 1.0d0/lattice_omega()
      call dcopy(nion,0.0d0,0,dbl_mb(fx(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fy(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fz(1)),1)
     
*     **** get the structure factor ****
      if (notzero_enpack) then
         call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            call daxpy(2*enpack,zii,dcpl_mb(exi(1)),1,dcpl_mb(s(1)),1)
         end do
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))

c         do i=1,enpack
c            dbl_mb(tmp3(1)+i-1) 
c     >              = ( dble(dcpl_mb(exi(1)+i-1))
c     >                *dimag(dcpl_mb(s(1)+i-1))
c     >              -  dimag(dcpl_mb(exi(1)+i-1))
c     >                 *dble(dcpl_mb(s(1)+i-1))
c     >                )*dbl_mb(vg(1)+i-1)
c         end do
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fx(1)+ii-1) = dbl_mb(fx(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fy(1)+ii-1) = dbl_mb(fy(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fz(1)+ii-1) = dbl_mb(fz(1)+ii-1) 
     >                         +  sum*zii*scal2
         end do
      end if
    
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
       if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut
           yerfc = util_erfc(w)
           f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
           sw1=sw1+(x*f)
           sw2=sw2+(y*f)
           sw3=sw3+(z*f)
        end do

        dbl_mb(fx(1)+i-1) = dbl_mb(fx(1)+i-1) + sw1
        dbl_mb(fy(1)+i-1) = dbl_mb(fy(1)+i-1) + sw2
        dbl_mb(fz(1)+i-1) = dbl_mb(fz(1)+i-1) + sw3

        dbl_mb(fx(1)+j-1) = dbl_mb(fx(1)+j-1) - sw1
        dbl_mb(fy(1)+j-1) = dbl_mb(fy(1)+j-1) - sw2
        dbl_mb(fz(1)+j-1) = dbl_mb(fz(1)+j-1) - sw3

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do
      if (tnp.gt.1) then
          call Parallel_Vector_SumAll(nion,dbl_mb(fx(1)))
          call Parallel_Vector_SumAll(nion,dbl_mb(fy(1)))
          call Parallel_Vector_SumAll(nion,dbl_mb(fz(1)))
      end if

      do i=1,nion
         fion(1,i) = fion(1,i) + dbl_mb(fx(1)+i-1)
         fion(2,i) = fion(2,i) + dbl_mb(fy(1)+i-1)
         fion(3,i) = fion(3,i) + dbl_mb(fz(1)+i-1)
      end do
      
*     **** deallocate temp workspace ****
      value =           BA_pop_stack(fz(2))
      value = value.and.BA_pop_stack(fy(2))
      value = value.and.BA_pop_stack(fx(2))
      if (notzero_enpack) then
         value = value.and.BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
      end if
      if (.not. value) 
     >  call errquit('ewald_f:popping stack memory',0,MA_ERR)

      call nwpw_timing_end(21)
      return
      end

*     ***********************************
*     *					*
*     *		ewald_stress_qm		*
*     *	        			*
*     ***********************************

      subroutine ewald_stress_qm(stress)
      implicit none
      real*8  stress(3,3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

      integer N
      parameter (N=40)

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)

*     **** local variables ****
      logical value
      integer i,ii,j,l,nion
      integer n1,n2,n3
      integer u,v,s,dutask
      real*8 pi,fourpi,scal
      real*8 zz,z
      real*8 Cus(3,3),hm(3,3),energy,sum,ss,rs
      real*8 ea,ax,ay,az,epsilon
      real*8 dx,dy,dz,w,ar,ai
      real*8 unita(3,3),unitg(3,3)
      complex*16 cz
      integer H(2),F(2),tmp1(2),tmp2(2),exi(2),strf(2)

*     **** external functions ****
      integer  ion_katm,ion_nion_qm
      real*8   ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      real*8   util_erfc,ion_rion,ddot
      external ion_katm,ion_nion_qm
      external ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      external util_erfc,ion_rion,ddot

      call nwpw_timing_start(21)
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)
      nion    = ion_nion_qm()
*     *** define hm,unita,unitg ****
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
         unitg(u,v) = lattice_unitg(u,v)
         unita(u,v) = lattice_unita(u,v)
      end do
      end do
 
      zz = 0.0d0
      z  = 0.0d0
      do i=1,nion
         zz = zz + dbl_mb(zv(1)+ion_katm(i)-1)**2
         z  = z  + dbl_mb(zv(1)+ion_katm(i)-1)
      end do

*     **** Miscellaneous contributions - stress from cewald_qm term ****
      do v=1,3
      do u=1,3
         stress(u,v) = 0.5d0*z*z*pi*rcut*rcut/lattice_omega()
     >               *hm(u,v)
      end do
      end do


*     **** G-space contributions ****
      if (notzero_enpack) then

*        **** get the structure factor ****
         value=          BA_push_get(mt_dbl,enpack,'H',H(2),H(1))
         value=value.and.BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1))
         value=value.and.
     >         BA_push_get(mt_dcpl,enpack,'strf',strf(2),strf(1))
         if (.not. value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)


         call dcopy((2*enpack),0.0d0,0,dcpl_mb(strf(1)),1)
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            call daxpy(2*enpack,ewald_zv(ion_katm(ii)),
     >                dcpl_mb(exi(1)),1,dcpl_mb(strf(1)),1)
         end do
         do i=1,enpack
             cz=dcpl_mb(strf(1)+i-1)
             ar=dble(cz)
             ai=dimag(cz)
             dbl_mb(H(1)+i-1)=ar*ar+ai*ai
         end do
         value =           BA_pop_stack(strf(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >    call errquit('ewald_stress:error popping stack',0,MA_ERR)

*        **** calculate the ewald energy ****
c         call dcopy(enpack,dbl_mb(vg(1)),1,dbl_mb(F(1)),1)
         F(1) = vg(1)
         energy=2.0d0*ddot(enpack,dbl_mb(F(1)),1,dbl_mb(H(1)),1)
      else
         energy=0.0d0
      end if
      if (tnp.gt.1) call Parallel_SumAll(energy)
      energy = -0.5d0*energy/lattice_omega()


      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) + energy*hm(u,v)
      end do
      end do
     
*     **** tmp2(G) = F(G)*H(G)/G**2 + F(G)*H(G)*rcut*rcut/4 ****
      if (notzero_enpack) then
         value=BA_push_get(mt_dbl,enpack,'tmp1',tmp1(2),tmp1(1))
         value=value.and.
     >         BA_push_get(mt_dbl,enpack,'tmp2',tmp2(2),tmp2(1))
         if (.not.value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)

         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(F(1)+i-1)*dbl_mb(H(1)+i-1)
         end do
         ss = 0.25d0*rcut*rcut
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         ss = 1.0d0/fourpi
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*dbl_mb(vcx(1)+i-1)
         end do       
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp2(1)+i-1)+dbl_mb(tmp1(1)+i-1)
         end do       
      end if

*     **** calculate Cus ****

      call dcopy(9,0.0d0,0,Cus,1)
      ss =  1.0d0/lattice_omega()
      if (notzero_enpack) then
         do u=1,3
         do s=u,3
            do i=1,enpack
               dbl_mb(tmp1(1)+i-1)=dbl_mb(eG(u)+i-1)*dbl_mb(eG(s)+i-1)
            end do
            sum=2.0d0*ddot(enpack,dbl_mb(tmp1(1)),1,dbl_mb(tmp2(1)),1)
            !if (tnp.gt.1) call Parallel_SumAll(sum)
            Cus(u,s) = ss*sum
         end do
         end do
c      do u=1,3
c      do s=u+1,3
c         Cus(s,u) = Cus(u,s)
c      end do
c      end do
c      do v=1,3
c      do u=1,3
c        do s=1,3
c           stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c        end do
c      end do
c      end do

         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(H(2))
         if (.not.value) 
     >    call errquit(
     >         'ewald_stress:error popping stack memory',0,MA_ERR)
      end if

*     **** R-space contributions ****

*     **** calculate alpha1 - stress from cewald_qm term*****
c        call dcopy(9,0.0d0,0,Cus,1)
        rs      = (3.0d0*lattice_omega()/(4.0d0*pi))**(1.0d0/3.0d0)
        epsilon = 1.0d0/rcut
        sum = 0.0d0
        dutask = 0
        do n1=(-N+1),(N-1) 
        do n2=(-N+1),(N-1) 
        do n3=(-N+1),(N-1) 
           if (dutask.eq.tid) then
           if (.not.((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0))) then
              ax = n1*unita(1,1)
     >           + n2*unita(1,2)
     >           + n3*unita(1,3)

              ay = n1*unita(2,1)
     >           + n2*unita(2,2)
     >           + n3*unita(2,3)

              az = n1*unita(3,1)
     >           + n2*unita(3,2)
     >           + n3*unita(3,3)

              ea = dsqrt(ax*ax + ay*ay + az*az)
              w = ea*epsilon

              ss = util_erfc(w)/ea 
     >           + 2.0d0*epsilon/dsqrt(pi)*dexp(-w*w)
              ss = -(0.5d0*zz)*ss/(ea*ea)
              Cus(1,1) = Cus(1,1) + ss * ax*ax 
              Cus(1,2) = Cus(1,2) + ss * ax*ay 
              Cus(1,3) = Cus(1,3) + ss * ax*az 

              !Cus(2,1) = Cus(2,1) + ss * ay*ax 
              Cus(2,2) = Cus(2,2) + ss * ay*ay 
              Cus(2,3) = Cus(2,3) + ss * ay*az 

              !Cus(3,1) = Cus(3,1) + ss * az*ax 
              !Cus(3,2) = Cus(3,2) + ss * az*ay 
              Cus(3,3) = Cus(3,3) + ss * az*az 

           end if
           end if
           dutask=mod(dutask+1,tnp)
        end do
        end do
        end do


c        do v=1,3
c        do u=1,3
c          do s=1,3
c             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c          end do
c        end do
c        end do


*     **** calculate erfc contribution *****
c      call dcopy(9,0.0d0,0,Cus,1)
      epsilon = 1.0d0/rcut
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
        if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
        do l=1,nshl3d
           ax = dbl_mb(rcell(1)+(l-1))          + dx
           ay = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           az = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           ea = dsqrt(ax*ax+ay*ay+az*az)
           w = ea*epsilon

           ss = -util_erfc(w)/ea
     >        - 2.0d0*epsilon/dsqrt(pi)*exp(-w*w)
           ss = ss/(ea*ea)
           Cus(1,1) = Cus(1,1) + ss * ax*ax * zz
           Cus(1,2) = Cus(1,2) + ss * ax*ay * zz
           Cus(1,3) = Cus(1,3) + ss * ax*az * zz
           Cus(2,2) = Cus(2,2) + ss * ay*ay * zz
           Cus(2,3) = Cus(2,3) + ss * ay*az * zz
           Cus(3,3) = Cus(3,3) + ss * az*az * zz
        end do
        end if
        dutask=mod(dutask+1,tnp)
      end do
      end do

      if (tnp.gt.1) then
        call Parallel_Vector_SumAll(9,Cus(1,1))
      end if

      do u=1,3
        do s=u+1,3
           Cus(s,u) = Cus(u,s)
        end do
      end do

      do v=1,3
        do u=1,3
          do s=1,3
             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
          end do
        end do
      end do

      call nwpw_timing_end(21)
      return
      end

*     ***********************************
*     *			  		*
*     *		ewald_e_mm		*
*     *				       	*
*     ***********************************
      real*8 function ewald_e_mm()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     **** local variables ****
      integer dutask
      integer i,j,ii,l,nion,nion_qm,nion_mm
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  energy,etmp,energyG,ttcr,ttci,zii
      complex*16 ttcz
*     **** temporary workspace variables ****
      integer exi(2),s(2),tmp3(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      real*8   ddot
      external ion_nion,ion_nion_qm,ion_nion_mm,ion_katm
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc
      external ddot

      call nwpw_timing_start(21)

*     **** allocate temp workspace ****
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_mm = ion_nion_mm()
      
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.
     >           BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if (.not. value) 
     >     call errquit('ewald_e:out of stack memory',0,MA_ERR)
     
*        **** get the structure factor ****
         call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
            call daxpy(2*enpack,dbl_mb(zv(1)+ion_katm(nion_qm+ii)-1),
     >             dcpl_mb(exi(1)),1,
     >             dcpl_mb(s(1)),1)
         end do

*        **** calculate the ewald energy ****
         call ewald_ct_Sqr(enpack,dcpl_mb(s(1)),dbl_mb(tmp3(1)))
         energy=2.0d0*ddot(enpack,dbl_mb(tmp3(1)),1,dbl_mb(vg(1)),1)
      else
         energy=0.0d0
      end if
      if (tnp.gt.1) call Parallel_SumAll(energy)
      energy  = 0.5d0*energy/lattice_omega() + cewald_mm

      
*     *** made parallel  ****
      dutask = 0
      etmp = 0.0d0
      do i=1,nion_mm-1
      do j=i+1,nion_mm
      if (dutask.eq.tid) then
        dx = ion_rion(1,nion_qm+i) - ion_rion(1,nion_qm+j)
        dy = ion_rion(2,nion_qm+i) - ion_rion(2,nion_qm+j)
        dz = ion_rion(3,nion_qm+i) - ion_rion(3,nion_qm+j)
        zz = dbl_mb(zv(1)+ion_katm(nion_qm+i)-1)
     >      *dbl_mb(zv(1)+ion_katm(nion_qm+j)-1)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = (r/rcut)
           yerfc = util_erfc(w)  !*** needs to be faster
           etmp=etmp+(zz*yerfc/r)
        end do
      end if
      dutask=mod(dutask+1,tnp)
      end do
      end do
      if (tnp.gt.1) call Parallel_SumAll(etmp)
      energy = energy + etmp
      

*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_e:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      ewald_e_mm = energy
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f_mm		*
*     *		       			*
*     ***********************************

      subroutine ewald_f_mm(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii

*     **** temporary workspace variables ****
      integer nion,nion_qm,nion_mm
      integer exi(2),s(2),tmp3(2)
      integer fx(2),fy(2),fz(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_nion_qm,ion_nion_mm,ion_katm,ewald_nshl3d
      external ion_nion,ion_nion_qm,ion_nion_mm,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot

      call nwpw_timing_start(21)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_mm = ion_nion_mm()

*     **** allocate temp workspace ****
      value = .true.
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
      end if
      value = value.and.BA_push_get(mt_dbl,nion_mm,'fx',fx(2),fx(1)) 
      value = value.and.BA_push_get(mt_dbl,nion_mm,'fy',fy(2),fy(1)) 
      value = value.and.BA_push_get(mt_dbl,nion_mm,'fz',fz(2),fz(1)) 
      if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)


      scal2 = 1.0d0/lattice_omega()
      call dcopy(nion_mm,0.0d0,0,dbl_mb(fx(1)),1)
      call dcopy(nion_mm,0.0d0,0,dbl_mb(fy(1)),1)
      call dcopy(nion_mm,0.0d0,0,dbl_mb(fz(1)),1)
     
*     **** get the structure factor ****
      if (notzero_enpack) then
         call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(nion_qm+ii)-1)
            call daxpy(2*enpack,zii,dcpl_mb(exi(1)),1,dcpl_mb(s(1)),1)
         end do
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
   
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

            zii=dbl_mb(zv(1)+ion_katm(nion_qm+ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fx(1)+ii-1) = dbl_mb(fx(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fy(1)+ii-1) = dbl_mb(fy(1)+ii-1) 
     >                         +  sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            dbl_mb(fz(1)+ii-1) = dbl_mb(fz(1)+ii-1) 
     >                         +  sum*zii*scal2
         end do
      end if
    
      dutask=0
      do i=1,nion_mm-1
      do j=i+1,nion_mm
       if (dutask.eq.tid) then
        dx = ion_rion(1,nion_qm+i) - ion_rion(1,nion_qm+j)
        dy = ion_rion(2,nion_qm+i) - ion_rion(2,nion_qm+j)
        dz = ion_rion(3,nion_qm+i) - ion_rion(3,nion_qm+j)
        zz = ewald_zv(ion_katm(nion_qm+i))*ewald_zv(ion_katm(nion_qm+j))
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut
           yerfc = util_erfc(w)
           f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
           sw1=sw1+(x*f)
           sw2=sw2+(y*f)
           sw3=sw3+(z*f)
        end do

        dbl_mb(fx(1)+i-1) = dbl_mb(fx(1)+i-1) + sw1
        dbl_mb(fy(1)+i-1) = dbl_mb(fy(1)+i-1) + sw2
        dbl_mb(fz(1)+i-1) = dbl_mb(fz(1)+i-1) + sw3

        dbl_mb(fx(1)+j-1) = dbl_mb(fx(1)+j-1) - sw1
        dbl_mb(fy(1)+j-1) = dbl_mb(fy(1)+j-1) - sw2
        dbl_mb(fz(1)+j-1) = dbl_mb(fz(1)+j-1) - sw3

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do
      if (tnp.gt.1) then
          call Parallel_Vector_SumAll(nion_mm,dbl_mb(fx(1)))
          call Parallel_Vector_SumAll(nion_mm,dbl_mb(fy(1)))
          call Parallel_Vector_SumAll(nion_mm,dbl_mb(fz(1)))
      end if

      do i=1,nion_mm
         fion(1,nion_qm+i) = fion(1,nion_qm+i) + dbl_mb(fx(1)+i-1)
         fion(2,nion_qm+i) = fion(2,nion_qm+i) + dbl_mb(fy(1)+i-1)
         fion(3,nion_qm+i) = fion(3,nion_qm+i) + dbl_mb(fz(1)+i-1)
      end do
      
*     **** deallocate temp workspace ****
      value =           BA_pop_stack(fz(2))
      value = value.and.BA_pop_stack(fy(2))
      value = value.and.BA_pop_stack(fx(2))
      if (notzero_enpack) then
         value = value.and.BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
      end if
      if (.not. value) 
     >  call errquit('ewald_f:popping stack memory',0,MA_ERR)

      call nwpw_timing_end(21)
      return
      end

*     ***********************************
*     *					*
*     *		ewald_stress_mm		*
*     *	        			*
*     ***********************************

      subroutine ewald_stress_mm(stress)
      implicit none
      real*8  stress(3,3)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

      integer N
      parameter (N=40)

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)

*     **** local variables ****
      logical value
      integer i,ii,j,l,nion,nion_qm,nion_mm
      integer n1,n2,n3
      integer u,v,s,dutask
      real*8 pi,fourpi,scal
      real*8 zz,z
      real*8 Cus(3,3),hm(3,3),energy,sum,ss,rs
      real*8 ea,ax,ay,az,epsilon
      real*8 dx,dy,dz,w,ar,ai
      real*8 unita(3,3),unitg(3,3)
      complex*16 cz
      integer H(2),F(2),tmp1(2),tmp2(2),exi(2),strf(2)

*     **** external functions ****
      integer  ion_katm,ion_nion,ion_nion_qm,ion_nion_mm
      real*8   ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      real*8   util_erfc,ion_rion,ddot
      external ion_katm,ion_nion,ion_nion_qm,ion_nion_mm
      external ewald_zv,lattice_unitg,lattice_unita,lattice_omega
      external util_erfc,ion_rion,ddot

      call nwpw_timing_start(21)
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_mm = ion_nion_mm()
*     *** define hm,unita,unitg ****
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
         unitg(u,v) = lattice_unitg(u,v)
         unita(u,v) = lattice_unita(u,v)
      end do
      end do
 
      zz = 0.0d0
      z  = 0.0d0
      do i=1,nion_mm
         zz = zz + dbl_mb(zv(1)+ion_katm(nion_qm+i)-1)**2
         z  = z  + dbl_mb(zv(1)+ion_katm(nion_qm+i)-1)
      end do

*     **** Miscellaneous contributions - stress from cewald_qm term ****
      do v=1,3
      do u=1,3
         stress(u,v) = 0.5d0*z*z*pi*rcut*rcut/lattice_omega()
     >               *hm(u,v)
      end do
      end do


*     **** G-space contributions ****

*     **** get the structure factor ****
      if (notzero_enpack) then
         value=BA_push_get(mt_dbl,enpack,'H',H(2),H(1))
         value=value.and.BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1))
         value=value.and.
     >         BA_push_get(mt_dcpl,enpack,'strf',strf(2),strf(1))
         if (.not. value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)


         call dcopy((2*enpack),0.0d0,0,dcpl_mb(strf(1)),1)
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
            call daxpy(2*enpack,ewald_zv(ion_katm(nion_qm+ii)),
     >                dcpl_mb(exi(1)),1,
     >                dcpl_mb(strf(1)),1)
         end do
         do i=1,enpack
             cz=dcpl_mb(strf(1)+i-1)
             ar=dble(cz)
             ai=dimag(cz)
             dbl_mb(H(1)+i-1)=ar*ar+ai*ai
         end do
         value =           BA_pop_stack(strf(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >    call errquit('ewald_stress:error popping stack',0,MA_ERR)

*        **** calculate the ewald energy ****
c         call dcopy(enpack,dbl_mb(vg(1)),1,dbl_mb(F(1)),1)
         F(1) = vg(1)
         energy=2.0d0*ddot(enpack,dbl_mb(F(1)),1,dbl_mb(H(1)),1)
      else
         energy=0.0d0
      end if
      if (tnp.gt.1) call Parallel_SumAll(energy)
      energy = -0.5d0*energy/lattice_omega()


      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) + energy*hm(u,v)
      end do
      end do
     
*     **** tmp2(G) = F(G)*H(G)/G**2 + F(G)*H(G)*rcut*rcut/4 ****
      if (notzero_enpack) then
         value=BA_push_get(mt_dbl,enpack,'tmp1',tmp1(2),tmp1(1))
         value=value.and.
     >         BA_push_get(mt_dbl,enpack,'tmp2',tmp2(2),tmp2(1))
         if (.not.value) 
     >     call errquit('ewald_stress:out of stack memory',0,MA_ERR)

         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(F(1)+i-1)*dbl_mb(H(1)+i-1)
         end do
         ss = 0.25d0*rcut*rcut
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         ss = 1.0d0/fourpi
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*ss
         end do       
         do i=1,enpack
            dbl_mb(tmp1(1)+i-1)=dbl_mb(tmp1(1)+i-1)*dbl_mb(vcx(1)+i-1)
         end do       
         do i=1,enpack
            dbl_mb(tmp2(1)+i-1)=dbl_mb(tmp2(1)+i-1)+dbl_mb(tmp1(1)+i-1)
         end do       
      end if

*     **** calculate Cus ****

      call dcopy(9,0.0d0,0,Cus,1)
      ss =  1.0d0/lattice_omega()
      if (notzero_enpack) then
         do u=1,3
         do s=u,3
            do i=1,enpack
               dbl_mb(tmp1(1)+i-1)=dbl_mb(eG(u)+i-1)*dbl_mb(eG(s)+i-1)
            end do
            sum=2.0d0*ddot(enpack,dbl_mb(tmp1(1)),1,dbl_mb(tmp2(1)),1)
            !if (tnp.gt.1) call Parallel_SumAll(sum)
            Cus(u,s) = ss*sum
         end do
         end do
c      do u=1,3
c      do s=u+1,3
c         Cus(s,u) = Cus(u,s)
c      end do
c      end do
c      do v=1,3
c      do u=1,3
c        do s=1,3
c           stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c        end do
c      end do
c      end do

         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(H(2))
         if (.not.value) 
     >    call errquit(
     >         'ewald_stress:error popping stack memory',0,MA_ERR)
      end if

*     **** R-space contributions ****

*     **** calculate alpha1 - stress from cewald_qm term*****
c        call dcopy(9,0.0d0,0,Cus,1)
        rs      = (3.0d0*lattice_omega()/(4.0d0*pi))**(1.0d0/3.0d0)
        epsilon = 1.0d0/rcut
        sum = 0.0d0
        dutask = 0
        do n1=(-N+1),(N-1) 
        do n2=(-N+1),(N-1) 
        do n3=(-N+1),(N-1) 
           if (dutask.eq.tid) then
           if (.not.((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0))) then
              ax = n1*unita(1,1)
     >           + n2*unita(1,2)
     >           + n3*unita(1,3)

              ay = n1*unita(2,1)
     >           + n2*unita(2,2)
     >           + n3*unita(2,3)

              az = n1*unita(3,1)
     >           + n2*unita(3,2)
     >           + n3*unita(3,3)

              ea = dsqrt(ax*ax + ay*ay + az*az)
              w = ea*epsilon

              ss = util_erfc(w)/ea 
     >           + 2.0d0*epsilon/dsqrt(pi)*dexp(-w*w)
              ss = -(0.5d0*zz)*ss/(ea*ea)
              Cus(1,1) = Cus(1,1) + ss * ax*ax 
              Cus(1,2) = Cus(1,2) + ss * ax*ay 
              Cus(1,3) = Cus(1,3) + ss * ax*az 

              !Cus(2,1) = Cus(2,1) + ss * ay*ax 
              Cus(2,2) = Cus(2,2) + ss * ay*ay 
              Cus(2,3) = Cus(2,3) + ss * ay*az 

              !Cus(3,1) = Cus(3,1) + ss * az*ax 
              !Cus(3,2) = Cus(3,2) + ss * az*ay 
              Cus(3,3) = Cus(3,3) + ss * az*az 

           end if
           end if
           dutask=mod(dutask+1,tnp)
        end do
        end do
        end do


c        do v=1,3
c        do u=1,3
c          do s=1,3
c             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
c          end do
c        end do
c        end do


*     **** calculate erfc contribution *****
c      call dcopy(9,0.0d0,0,Cus,1)
      epsilon = 1.0d0/rcut
      dutask=0
      do i=1,nion_mm-1
      do j=i+1,nion_mm
        if (dutask.eq.tid) then
        dx = ion_rion(1,nion_qm+i) - ion_rion(1,nion_qm+j)
        dy = ion_rion(2,nion_qm+i) - ion_rion(2,nion_qm+j)
        dz = ion_rion(3,nion_qm+i) - ion_rion(3,nion_qm+j)
        zz = ewald_zv(ion_katm(nion_qm+i))*ewald_zv(ion_katm(nion_qm+j))
        do l=1,nshl3d
           ax = dbl_mb(rcell(1)+(l-1))          + dx
           ay = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           az = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           ea = dsqrt(ax*ax+ay*ay+az*az)
           w = ea*epsilon

           ss = -util_erfc(w)/ea
     >        - 2.0d0*epsilon/dsqrt(pi)*exp(-w*w)
           ss = ss/(ea*ea)
           Cus(1,1) = Cus(1,1) + ss * ax*ax * zz
           Cus(1,2) = Cus(1,2) + ss * ax*ay * zz
           Cus(1,3) = Cus(1,3) + ss * ax*az * zz
           Cus(2,2) = Cus(2,2) + ss * ay*ay * zz
           Cus(2,3) = Cus(2,3) + ss * ay*az * zz
           Cus(3,3) = Cus(3,3) + ss * az*az * zz
        end do
        end if
        dutask=mod(dutask+1,tnp)
      end do
      end do

      if (tnp.gt.1) then
        call Parallel_Vector_SumAll(9,Cus(1,1))
      end if

      do u=1,3
        do s=u+1,3
           Cus(s,u) = Cus(u,s)
        end do
      end do

      do v=1,3
        do u=1,3
          do s=1,3
             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
          end do
        end do
      end do

      call nwpw_timing_end(21)
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f_async		*
*     *		       			*
*     ***********************************

      subroutine ewald_f_async(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii
      common /ewald_f_async_sw1s2/ sw1,sw2,sw3

*     **** temporary workspace variables ****
      integer nion
      integer exi(2),s(2),tmp3(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,ewald_nshl3d
      external ion_nion,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot

      call nwpw_timing_start(21)
      nion = ion_nion()
      scal2 = 1.0d0/lattice_omega()

*     **** allocate temp workspace ****
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)

*        **** get the structure factor ****
         !call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         call Parallel_shared_vector_zero(.true.,2*enpack,dcpl_mb(s(1)))
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            call daxpy_omp(2*enpack,zii,dcpl_mb(exi(1)),1,
     >                     dcpl_mb(s(1)),1)
         end do
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))

c         do i=1,enpack
c            dbl_mb(tmp3(1)+i-1) 
c     >              = ( dble(dcpl_mb(exi(1)+i-1))
c     >                *dimag(dcpl_mb(s(1)+i-1))
c     >              -  dimag(dcpl_mb(exi(1)+i-1))
c     >                 *dble(dcpl_mb(s(1)+i-1))
c     >                )*dbl_mb(vg(1)+i-1)
c         end do
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

!$OMP MASTER
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            fion(1,ii) = fion(1,ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            fion(2,ii) = fion(2,ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            fion(3,ii) = fion(3,ii) + sum*zii*scal2
!$OMP END MASTER
!$OMP BARRIER
         end do
      end if
    
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
       if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
!$OMP MASTER
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sw1,sw2,sw3)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           if (r.gt.1.0d-6) then
              w = r/rcut
              yerfc = util_erfc(w)
              f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
              sw1=sw1+(x*f)
              sw2=sw2+(y*f)
              sw3=sw3+(z*f)
           end if
        end do
!$OMP END DO

!$OMP MASTER
        fion(1,i) = fion(1,i) + sw1
        fion(2,i) = fion(2,i) + sw2
        fion(3,i) = fion(3,i) + sw3

        fion(1,j) = fion(1,j) - sw1
        fion(2,j) = fion(2,j) - sw2
        fion(3,j) = fion(3,j) - sw3
!$OMP END MASTER
!$OMP BARRIER

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do

      
*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_f_async:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f_qm_async	*
*     *		       			*
*     ***********************************

      subroutine ewald_f_qm_async(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii
      common /ewald_f_qm_async_sw1sw2/ sw1,sw2,sw3

*     **** temporary workspace variables ****
      integer nion
      integer exi(2),s(2),tmp3(2)
      logical value

*     **** external functions ****
      integer  ion_nion_qm,ion_katm,ewald_nshl3d
      external ion_nion_qm,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot

      call nwpw_timing_start(21)
      nion = ion_nion_qm()
      scal2 = 1.0d0/lattice_omega()

*     **** allocate temp workspace ****
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)

*        **** get the structure factor ****
         !call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         call Parallel_shared_vector_zero(.true.,2*enpack,dcpl_mb(s(1)))
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            call daxpy_omp(2*enpack,zii,dcpl_mb(exi(1)),1,
     >                     dcpl_mb(s(1)),1)
         end do
         do ii=1,nion
            call ewald_strfac(ii,dcpl_mb(exi(1)))

c         do i=1,enpack
c            dbl_mb(tmp3(1)+i-1) 
c     >              = ( dble(dcpl_mb(exi(1)+i-1))
c     >                *dimag(dcpl_mb(s(1)+i-1))
c     >              -  dimag(dcpl_mb(exi(1)+i-1))
c     >                 *dble(dcpl_mb(s(1)+i-1))
c     >                )*dbl_mb(vg(1)+i-1)
c         end do
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

!$OMP MASTER
            zii=dbl_mb(zv(1)+ion_katm(ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            fion(1,ii) = fion(1,ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            fion(2,ii) = fion(2,ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            fion(3,ii) = fion(3,ii) + sum*zii*scal2
!$OMP END MASTER
!$OMP BARRIER
         end do
      end if
    
      dutask=0
      do i=1,nion-1
      do j=i+1,nion
       if (dutask.eq.tid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = ewald_zv(ion_katm(i)) * ewald_zv(ion_katm(j))
!$OMP MASTER
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sw1,sw2,sw3)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut
           yerfc = util_erfc(w)
           f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
           sw1=sw1+(x*f)
           sw2=sw2+(y*f)
           sw3=sw3+(z*f)
        end do
!$OMP END DO

!$OMP MASTER
        fion(1,i) = fion(1,i) + sw1
        fion(2,i) = fion(2,i) + sw2
        fion(3,i) = fion(3,i) + sw3

        fion(1,j) = fion(1,j) - sw1
        fion(2,j) = fion(2,j) - sw2
        fion(3,j) = fion(3,j) - sw3
!$OMP END MASTER

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do

      
*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_f:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      return
      end


*     ***********************************
*     *	        			*
*     *		ewald_f_mm_async	*
*     *		       			*
*     ***********************************

      subroutine ewald_f_mm_async(fion)
      implicit none
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3,zii
      common /ewald_f_mm_async_sw1sw2/ sw1,sw2,sw3

*     **** temporary workspace variables ****
      integer nion,nion_qm,nion_mm
      integer exi(2),s(2),tmp3(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_nion_qm,ion_nion_mm,ion_katm,ewald_nshl3d
      external ion_nion,ion_nion_qm,ion_nion_mm,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot
      external lattice_omega,ewald_zv,dsum,ion_rion,util_erfc,ddot

      call nwpw_timing_start(21)
      nion    = ion_nion()
      nion_qm = ion_nion_qm()
      nion_mm = ion_nion_mm()
      scal2 = 1.0d0/lattice_omega()

*     **** allocate temp workspace ****
      if (notzero_enpack) then
         value = BA_push_get(mt_dcpl,enpack,'exi',exi(2),exi(1)) 
         value = value.and.BA_push_get(mt_dcpl,enpack,'s',s(2),s(1)) 
         value = value.and.
     >           BA_push_get(mt_dbl, enpack,'tmp3',tmp3(2),tmp3(1)) 
         if(.not.value) call errquit('ewald_f:out of stack',0,MA_ERR)
     
*        **** get the structure factor ****
         !call dcopy((2*enpack),0.0d0,0,dcpl_mb(s(1)),1)
         call Parallel_shared_vector_zero(.true.,2*enpack,dcpl_mb(s(1)))
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
            zii=dbl_mb(zv(1)+ion_katm(nion_qm+ii)-1)
            call daxpy_omp(2*enpack,zii,dcpl_mb(exi(1)),1,
     >                     dcpl_mb(s(1)),1)
         end do
         do ii=1,nion_mm
            call ewald_strfac(nion_qm+ii,dcpl_mb(exi(1)))
   
            call ewald_f_tmp3(enpack,
     >                        dcpl_mb(exi(1)),dcpl_mb(s(1)),
     >                        dbl_mb(vg(1)),dbl_mb(tmp3(1)))

!$OMP MASTER
            zii=dbl_mb(zv(1)+ion_katm(nion_qm+ii)-1)
            sum=2.0d0*ddot(enpack,dbl_mb(eG(1)),1,dbl_mb(tmp3(1)),1)
            fion(1,nion_qm+ii) = fion(1,nion_qm+ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(2)),1,dbl_mb(tmp3(1)),1)
            fion(2,nion_qm+ii) = fion(2,nion_qm+ii) + sum*zii*scal2
            sum=2.0d0*ddot(enpack,dbl_mb(eG(3)),1,dbl_mb(tmp3(1)),1)
            fion(3,nion_qm+ii) = fion(3,nion_qm+ii) + sum*zii*scal2
!$OMP END MASTER
!$OMP BARRIER
         end do
      end if
    
      dutask=0
      do i=1,nion_mm-1
      do j=i+1,nion_mm
       if (dutask.eq.tid) then
        dx = ion_rion(1,nion_qm+i) - ion_rion(1,nion_qm+j)
        dy = ion_rion(2,nion_qm+i) - ion_rion(2,nion_qm+j)
        dz = ion_rion(3,nion_qm+i) - ion_rion(3,nion_qm+j)
        zz = ewald_zv(ion_katm(nion_qm+i))*ewald_zv(ion_katm(nion_qm+j))
!$OMP MASTER
        sw1=0.0d0
        sw2=0.0d0
        sw3=0.0d0  
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sw1,sw2,sw3)
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut
           yerfc = util_erfc(w)
           f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
           sw1=sw1+(x*f)
           sw2=sw2+(y*f)
           sw3=sw3+(z*f)
        end do
!$OMP END DO

!$OMP MASTER
        fion(1,nion_qm+i) = fion(1,nion_qm+i) + sw1
        fion(2,nion_qm+i) = fion(2,nion_qm+i) + sw2
        fion(3,nion_qm+i) = fion(3,nion_qm+i) + sw3

        fion(1,nion_qm+j) = fion(1,nion_qm+j) - sw1
        fion(2,nion_qm+j) = fion(2,nion_qm+j) - sw2
        fion(3,nion_qm+j) = fion(3,nion_qm+j) - sw3
!$OMP END MASTER
!$OMP BARRIER

       end if
       dutask = mod((dutask+1),tnp)
      end do
      end do

*     **** deallocate temp workspace ****
      if (notzero_enpack) then
         value =           BA_pop_stack(tmp3(2))
         value = value.and.BA_pop_stack(s(2))
         value = value.and.BA_pop_stack(exi(2))
         if (.not. value) 
     >     call errquit('ewald_f:popping stack memory',0,MA_ERR)
      end if

      call nwpw_timing_end(21)
      return
      end




*     ***********************************
*     *	        			*
*     *		ewald_efg		*
*     *		       			*
*     ***********************************

      subroutine ewald_efg(efg_ion)
      implicit none
      real*8  efg_ion(3,3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     **** local variables ****
      integer dutask
      integer ii,jj,kk,l,k,u,v,nion
      real*8  w,dx,dy,dz,x,y,z,r,qj
      real*8  yerfc,kfac
      real*8  aconst,bconst,scal2,sqrt_pi
      integer mu,nu
      real*8  gvec(3),gg,phase,cosphase
      real*8  rji(3),termgg,pi,expfac

*     **** external functions ****
      integer  ion_nion,ion_katm,ewald_nshl3d
      external ion_nion,ion_katm,ewald_nshl3d
      real*8   lattice_omega,ewald_zv,ion_rion,util_erfc
      external lattice_omega,ewald_zv,ion_rion,util_erfc

*     **** preliminaries ****
      call nwpw_timing_start(21)
      nion = ion_nion()
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      pi = sqrt_pi*sqrt_pi

*     **** initialize efg_ion ****
      call dcopy(9*nion,0.0d0,0,efg_ion,1)

*     **** reciprocal space part of efg_ion ****
*     **** Honma eq. 4.8 ****
      aconst = 4.d0*pi/3.d0/lattice_omega()
      if (notzero_enpack) then  ! eliminate g=0 term
         do k=1+nida,enpack
           gvec(1) = dbl_mb(eG(1)+k-1)
           gvec(2) = dbl_mb(eG(2)+k-1)
           gvec(3) = dbl_mb(eG(3)+k-1)
           gg=gvec(1)*gvec(1)+gvec(2)*gvec(2)+gvec(3)*gvec(3)
           expfac = aconst*dbl_mb(vg(1)+k-1)/dbl_mb(vcx(1)+k-1)  ! see ewald_init for exp factors
           do ii=1,nion
            do jj=1,nion
             qj = ewald_zv(ion_katm(jj))
             rji(1) = ion_rion(1,jj) - ion_rion(1,ii)
             rji(2) = ion_rion(2,jj) - ion_rion(2,ii)
             rji(3) = ion_rion(3,jj) - ion_rion(3,ii)
             phase = gvec(1)*rji(1)+gvec(2)*rji(2)+gvec(3)*rji(3)
             cosphase = dcos(phase)
             do mu= 1,3
              do nu= 1,3
               termgg = -3.d0*gvec(mu)*gvec(nu)/gg
               if (mu == nu) termgg = 1.d0 + termgg

               !**** extra factor of 2 needed because its an integral ****
               efg_ion(mu,nu,ii) = efg_ion(mu,nu,ii) 
     &                           + 2.0d0*qj*cosphase*termgg*expfac
              end do
             end do
            end do ! jj
           end do ! ii
         end do  ! k
      end if   ! notzero_enpack
!
*     **** real-space part of efg_ion ****
*     **** Honma eq. 2.10 ****
      bconst = 4.0d0/(3.0d0*sqrt_pi*rcut**3)
      dutask=0
      do ii=1,nion
      do jj=1,nion
       if (dutask.eq.tid) then
        dx = ion_rion(1,ii) - ion_rion(1,jj)
        dy = ion_rion(2,ii) - ion_rion(2,jj)
        dz = ion_rion(3,ii) - ion_rion(3,jj)
        qj = ewald_zv(ion_katm(jj))
        do l=1,nshl3d
         x = dbl_mb(rcell(1)+(l-1))          + dx
         y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
         z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
         r = dsqrt(x*x+y*y+z*z)
         if (r.gt.1.0d-6) then
          w = r/rcut
          yerfc = util_erfc(w)
          kfac = -bconst*((1.0d0+3.0d0/(2.0d0*w*w))*exp(-w*w) 
     >            + 3.0d0*sqrt_pi*yerfc/(4.0d0*w**3))
          efg_ion(1,1,ii)=efg_ion(1,1,ii)+qj*kfac*(1.d0-3.d0*(x*x)/r**2)
          efg_ion(2,1,ii)=efg_ion(2,1,ii)+qj*kfac*(    -3.d0*(y*x)/r**2)
          efg_ion(3,1,ii)=efg_ion(3,1,ii)+qj*kfac*(    -3.d0*(z*x)/r**2)
          efg_ion(1,2,ii)=efg_ion(1,2,ii)+qj*kfac*(    -3.d0*(x*y)/r**2)
          efg_ion(2,2,ii)=efg_ion(2,2,ii)+qj*kfac*(1.d0-3.d0*(y*y)/r**2)
          efg_ion(3,2,ii)=efg_ion(3,2,ii)+qj*kfac*(    -3.d0*(z*y)/r**2)
          efg_ion(1,3,ii)=efg_ion(1,3,ii)+qj*kfac*(    -3.d0*(x*z)/r**2)
          efg_ion(2,3,ii)=efg_ion(2,3,ii)+qj*kfac*(    -3.d0*(y*z)/r**2)
          efg_ion(3,3,ii)=efg_ion(3,3,ii)+qj*kfac*(1.d0-3.d0*(z*z)/r**2)
         end if
        end do !** lattice vectors **
       end if ! dutask
       dutask = mod((dutask+1),tnp)
      end do !**jj**
      end do !**ii**
!
      if (tnp.gt.1) then
        call Parallel_Vector_SumAll(9*nion,efg_ion)
      end if
!
      call nwpw_timing_end(21)
!
      return
      end
