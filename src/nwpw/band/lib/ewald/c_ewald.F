*
* $Id: c_ewald.F,v 1.1 2001-12-13 03:53:32 bylaska Exp $
*
      integer function c_ewald_ncut()
      implicit none 


*     **** common block for c_ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut
   
      c_ewald_ncut = ncut
      return
      end 

      real*8 function c_ewald_rcut()
      implicit none 

*     **** common block for c_ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut
   
      c_ewald_rcut = rcut
      return
      end 

      integer function c_ewald_nshl3d()
      implicit none 


*     **** common block for c_ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut
   
      c_ewald_nshl3d = (2*ncut+1)**3
      return
      end 


      real*8 function c_ewald_mandelung()
      implicit none 


*     **** common block for ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut
   
      c_ewald_mandelung = alpha
      return
      end 


      subroutine c_ewald_end()
      implicit none 

#include "mafdecls.fh"

*     **** common block for ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut

      logical value

      value = MA_free_heap(vg(2))
      value = MA_free_heap(rcell(2))

      return
      end

 
      subroutine c_ewald_init()
      implicit none 

#include "mafdecls.fh"


*     **** common block for ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut

*     **** local variables ****
      integer nfft3d,G(3)
      integer nshl3d
      integer i,j,k,l
      real*8  pi,fourpi,gg,w
      real*8  rs
      real*8  zz,z
      integer taskid,pzero,qzero,zero
      integer nx,ny,nxh
      logical value
      real*8 kv(3),ecut

*     **** external functions ****
      integer  control_ncut
      real*8   control_rcut,control_ecut
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega,lattice_unita,cpsp_zv

      external control_ncut
      external control_rcut,control_ecut
      external ion_nion,ion_katm,G_indx
      external lattice_omega,lattice_unita,cpsp_zv


*     **** allocate vg memory ****
      call C3dB_nfft3d(1,nfft3d)
      value = MA_alloc_get(mt_dbl,nfft3d,'vg',vg(2),vg(1))
      if (.not. value) call errquit('out of heap memory',0)

      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

*     **** get constants ****
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi

      call Parallel_taskid(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      nxh=nx/2

*     ***** find the G==0 index ******
      i=0
      j=0
      k=0
      call C3dB_ktoqp(1,k+1,qzero,pzero)
      zero = (qzero-1)*(nx)*ny
     >     + j*(nx)
     >     + i+1
      

*     ***** initialize common block and find w *****
      ncut = control_ncut()
      rcut = control_rcut()
      if (ncut.le.0)     ncut=1
      if (rcut.le.0.0d0) then
         rs = lattice_unita(1,1)**2
     >      + lattice_unita(2,1)**2
     >      + lattice_unita(3,1)**2
         rs = dsqrt(rs)
         rcut=rs/pi
         
         rs = lattice_unita(1,2)**2
     >      + lattice_unita(2,2)**2
     >      + lattice_unita(3,2)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w

         rs = lattice_unita(1,3)**2
     >      + lattice_unita(2,3)**2
     >      + lattice_unita(3,3)**2
         rs = dsqrt(rs)
         w=rs/pi
         if (w.lt.rcut) rcut = w
      end if 

      w      = 0.25d0*rcut*rcut


*     ***** initialize Vg  *****
      do i=1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
      
         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            dbl_mb(vg(1)+i-1) = 0.0d0
         else
            dbl_mb(vg(1)+i-1) = (fourpi/gg)*exp(-w*gg)
         end if
      end do

      kv(1) = 0.0d0
      kv(2) = 0.0d0
      kv(3) = 0.0d0
      ecut = control_ecut()
      call cloak_set(kv,ecut)
      call cloak_R(dbl_mb(vg(1)))

    
*     **** set the Mandelung constant ****
      call mandelung_set(alpha)

      
*     **** ewald summation ****
      rs = (3.0d0*lattice_omega()/fourpi)**(1.0d0/3.0d0)
      zz = 0.0d0
      z  = 0.0d0
      do i=1,ion_nion()
         zz = zz + cpsp_zv(ion_katm(i))**2
         z  = z  + cpsp_zv(ion_katm(i))
      end do
      call C3dB_r_dsum(1,dbl_mb(vg(1)),cewald)
      cewald = -0.5d0*zz*(alpha/rs + cewald/lattice_omega())
     >         -0.5d0*(z*z-zz)*rcut*rcut*pi/lattice_omega()
     
*     **** allocate rcell memory ****
      nshl3d=(2*ncut+1)**3
      value = MA_alloc_get(mt_dbl,(3*nshl3d),'rcell',rcell(2),
     >                                           rcell(1))
      if (.not. value) call errquit('out of heap memory',0)

      
*     **** get lattice vectors in real space ****
      l=0
      do k=-ncut,ncut
        do j=-ncut,ncut
          do i=-ncut,ncut
             l = l+1
             dbl_mb(rcell(1)+ (l-1) ) 
     >                = i*lattice_unita(1,1)
     >                + j*lattice_unita(1,2)
     >                + k*lattice_unita(1,3)
             dbl_mb(rcell(1)+(l-1)+nshl3d) 
     >                = i*lattice_unita(2,1)
     >                + j*lattice_unita(2,2)
     >                + k*lattice_unita(2,3)
             dbl_mb(rcell(1)+(l-1)+2*nshl3d) 
     >                = i*lattice_unita(3,1)
     >                + j*lattice_unita(3,2)
     >                + k*lattice_unita(3,3)

          end do
        end do
      end do
      

      return
      end

*     ***********************************
*     *				        			*
*     *				c_ewald_e				*
*     *				        			*
*     ***********************************
      real*8 function c_ewald_e()
      implicit none

#include "mafdecls.fh"



*     **** common block for ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut


*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,l
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  energy,etmp

*     **** temporary workspace variables ****
c     complex*16  exi(nfft3d)
c     complex*16    s(nfft3d)
c     real*8     tmp3(nfft3d*2)
      integer nfft3d,nshl3d
      integer exi(2),s(2),tmp3(2),ft(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,c_ewald_nshl3d
      real*8   lattice_omega,cpsp_zv,dsum,ion_rion,util_erfc
      external ion_nion,ion_katm,c_ewald_nshl3d
      external lattice_omega,cpsp_zv,dsum,ion_rion,util_erfc

      call Parallel_np(np)
      call Parallel_taskid(taskid)

*     **** allocate temp workspace ****
      call C3dB_nfft3d(1,nfft3d)
      nshl3d = c_ewald_nshl3d()
      value = MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'s',s(2),s(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'tmp3',tmp3(2),tmp3(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (3*nshl3d),'ft',ft(2),ft(1)) 
      if (.not. value) call errquit('out of stack memory',0)
     
*     **** get the structure factor ****
      call dcopy((2*nfft3d),0.0d0,0,dcpl_mb(s(1)),1)
      do ii=1,ion_nion()
         call cstrfac(ii,dcpl_mb(exi(1)))
         call C3dB_cc_daxpy(1,cpsp_zv(ion_katm(ii)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(s(1)))

      end do

      
*     **** calculate the ewald energy ****
      call C3dB_cr_Sqr(1,dcpl_mb(s(1)),dbl_mb(tmp3(1)))
      call C3dB_rr_dot(1,dbl_mb(tmp3(1)),dbl_mb(vg(1)),energy)
      energy = 0.5d0*energy/lattice_omega() + cewald

*     *** need to make parallel ****
      dutask = 0
      etmp = 0.0d0
      do i=1,ion_nion()-1
      do j=i+1,ion_nion()
      if (dutask.eq.taskid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = cpsp_zv(ion_katm(i)) * cpsp_zv(ion_katm(j))
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+nshl3d)   + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut

c          erfc=1.0d0/(1.0d0+w*(b1+w*(b2+w*(b3
c    >                   +w*(b4+w*(b5+w*b6))))))**4
c          dbl_mb(ft(1)+(l-1))=zz*erfc**4/r
           yerfc = util_erfc(w)
           dbl_mb(ft(1)+(l-1))=zz*yerfc/r
        end do
        etmp = etmp + dsum(nshl3d,dbl_mb(ft(1)),1)
      end if
      dutask = mod(dutask+1,np)
      end do
      end do
      if (np.gt.1) call C3dB_SumAll(etmp)
      energy = energy + etmp
      

*     **** deallocate temp workspace ****
      value = MA_pop_stack(ft(2))
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(s(2))
      value = MA_pop_stack(exi(2))

      c_ewald_e = energy
      return
      end


*     ***********************************
*     *				        			*
*     *				c_ewald_f				*
*     *				        			*
*     ***********************************

      subroutine c_ewald_f(fion)
      implicit none
      real*8  fion(3,*)

#include "mafdecls.fh"

*     **** common block for ewald.f ****
      integer    ncut
      real*8     rcut,cewald,alpha
      integer    vg(2),rcell(2)
      common / c_ewald_block / vg,rcell,cewald,alpha,rcut,ncut


*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)


*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,l,ii
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  sum,scal2,f
      real*8  sw1,sw2,sw3

*     **** temporary workspace variables ****
c     complex*16  exi(nfft3d)
c     complex*16    s(nfft3d)
c     real*8     tmp3(nfft3d*2)
      integer nfft3d,nshl3d,nion
      integer exi(2),s(2),tmp3(2),ft(2)
      integer fx(2),fy(2),fz(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,G_indx,c_ewald_nshl3d
      external ion_nion,ion_katm,G_indx,c_ewald_nshl3d
      real*8   lattice_omega,cpsp_zv,dsum,ion_rion,util_erfc
      external lattice_omega,cpsp_zv,dsum,ion_rion,util_erfc

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nion = ion_nion()

*     **** allocate temp workspace ****
      call C3dB_nfft3d(1,nfft3d)
      nshl3d = c_ewald_nshl3d()
      value = MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'s',s(2),s(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'tmp3',tmp3(2),tmp3(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (3*nshl3d),'ft',ft(2),ft(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fx',fx(2),fx(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fy',fy(2),fy(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fz',fz(2),fz(1)) 
      if (.not. value) call errquit('out of heap memory',0)


      scal2 = 1.0d0/lattice_omega()
      call dcopy(nion,0.0d0,0,dbl_mb(fx(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fy(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fz(1)),1)
     
*     **** get the structure factor ****
      call dcopy((2*nfft3d),0.0d0,0,dcpl_mb(s(1)),1)
      do ii=1,nion
         call cstrfac(ii,dcpl_mb(exi(1)))
         call C3dB_cc_daxpy(1,cpsp_zv(ion_katm(ii)),dcpl_mb(exi(1)),
     >                                           dcpl_mb(s(1)))
      end do

      do ii=1,nion
         call cstrfac(ii,dcpl_mb(exi(1)))

         do i=1,nfft3d
            dbl_mb(tmp3(1)+i-1) 
     >              = ( dble(dcpl_mb(exi(1)+i-1))
     >                *dimag(dcpl_mb(s(1)+i-1))
     >              -  dimag(dcpl_mb(exi(1)+i-1))
     >                 *dble(dcpl_mb(s(1)+i-1))
     >                )*dbl_mb(vg(1)+i-1)
         end do

         call C3dB_rr_idot(1,dbl_mb(G_indx(1)),dbl_mb(tmp3(1)),sum)
*        fion(1,ii) = fion(1,ii) + sum*cpsp_zv(ion_katm(ii))*scal2
         dbl_mb(fx(1)+ii-1) = dbl_mb(fx(1)+ii-1) 
     >                      +  sum*cpsp_zv(ion_katm(ii))*scal2

         call C3dB_rr_idot(1,dbl_mb(G_indx(2)),dbl_mb(tmp3(1)),sum)
*        fion(2,ii) = fion(2,ii) + sum*cpsp_zv(ion_katm(ii))*scal2
         dbl_mb(fy(1)+ii-1) = dbl_mb(fy(1)+ii-1) 
     >                      +  sum*cpsp_zv(ion_katm(ii))*scal2

         call C3dB_rr_idot(1,dbl_mb(G_indx(3)),dbl_mb(tmp3(1)),sum)
*        fion(3,ii) = fion(3,ii) + sum*cpsp_zv(ion_katm(ii))*scal2
         dbl_mb(fz(1)+ii-1) = dbl_mb(fz(1)+ii-1) 
     >                      +  sum*cpsp_zv(ion_katm(ii))*scal2
      end do
    
      

      dutask=0
      do i=1,nion-1
      do j=i+1,nion
        if (dutask.eq.taskid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
        zz = cpsp_zv(ion_katm(i)) * cpsp_zv(ion_katm(j))
        do l=1,nshl3d
           x = dbl_mb(rcell(1)+(l-1))          + dx
           y = dbl_mb(rcell(1)+(l-1)+  nshl3d) + dy
           z = dbl_mb(rcell(1)+(l-1)+2*nshl3d) + dz
           r = dsqrt(x*x+y*y+z*z)
           w = r/rcut

c          erfc=(1.0d0+w*(b1+w*(b2+w*(b3
c    >                   +w*(b4+w*(b5+w*b6))))))**4
c          erfc = 1.0d0/erfc**4
           yerfc = util_erfc(w)
           f = zz*(yerfc+cerfc*w*dexp(-w*w))/r**3
           dbl_mb(ft(1)+(l-1))         =x*f
           dbl_mb(ft(1)+(l-1)+nshl3d)  =y*f
           dbl_mb(ft(1)+(l-1)+2*nshl3d)=z*f
        end do
        sw1 = dsum(nshl3d,dbl_mb(ft(1)),1)
        sw2 = dsum(nshl3d,dbl_mb(ft(1)+  nshl3d),1)
        sw3 = dsum(nshl3d,dbl_mb(ft(1)+2*nshl3d),1)

*       fion(1,i) = fion(1,i) + sw1
*       fion(2,i) = fion(2,i) + sw2
*       fion(3,i) = fion(3,i) + sw3
*       fion(1,j) = fion(1,j) - sw1
*       fion(2,j) = fion(2,j) - sw2
*       fion(3,j) = fion(3,j) - sw3

        dbl_mb(fx(1)+i-1) = dbl_mb(fx(1)+i-1) + sw1
        dbl_mb(fy(1)+i-1) = dbl_mb(fy(1)+i-1) + sw2
        dbl_mb(fz(1)+i-1) = dbl_mb(fz(1)+i-1) + sw3

        dbl_mb(fx(1)+j-1) = dbl_mb(fx(1)+j-1) - sw1
        dbl_mb(fy(1)+j-1) = dbl_mb(fy(1)+j-1) - sw2
        dbl_mb(fz(1)+j-1) = dbl_mb(fz(1)+j-1) - sw3

      end if
      dutask = mod((dutask+1),np)
      end do
      end do

      if (np.gt.1) call C3dB_Vector_SumAll(nion,dbl_mb(fx(1)))
      if (np.gt.1) call C3dB_Vector_SumAll(nion,dbl_mb(fy(1)))
      if (np.gt.1) call C3dB_Vector_SumAll(nion,dbl_mb(fz(1)))
      do i=1,nion
         fion(1,i) = fion(1,i) + dbl_mb(fx(1)+i-1)
         fion(2,i) = fion(2,i) + dbl_mb(fy(1)+i-1)
         fion(3,i) = fion(3,i) + dbl_mb(fz(1)+i-1)
      end do
      
*     **** deallocate temp workspace ****
      value = MA_pop_stack(fz(2))
      value = MA_pop_stack(fy(2))
      value = MA_pop_stack(fx(2))
      value = MA_pop_stack(ft(2))
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(s(2))
      value = MA_pop_stack(exi(2))

      return
      end

