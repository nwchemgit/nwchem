*
* $Id: ewald.F,v 1.14 2009-02-07 03:50:50 bylaska Exp $
*
      integer function ewald_ncut()
      implicit none 
#include "ewald.fh"
      ewald_ncut = ncut
      return
      end 

      real*8 function ewald_rcut()
      implicit none 
#include "ewald.fh"
      ewald_rcut = rcut
      return
      end 

      integer function ewald_nshl3d()
      implicit none 
#include "ewald.fh"
      ewald_nshl3d = (2*ncut+1)**3
      return
      end 

      integer function ewald_rcell_ptr()
      implicit none
#include "ewald.fh"
      ewald_rcell_ptr = rcell(1)
      return
      end


      real*8 function ewald_mandelung()
      implicit none 
#include "ewald.fh"
      ewald_mandelung = alpha
      return
      end 


      subroutine mandelung_set(alpha)
      implicit none
      real*8 alpha
      
      integer N
      parameter (N=40)

*     **** local variables ****
      integer n1,n2,n3
      real*8  rc,rs,epsilon,pi
      real*8  a1,a2,a3,g1,g2,g3,gg
      real*8  unita(3,3),unitg(3,3)
      real*8  alpha1,alpha2,sum,ea,omega

*     ****  expansion coefficient of the error function ****
c     real*8 cerfc,b1,b2,b3,b4,b5,b6,yerfc,xerfc
c     parameter (cerfc=1.128379167d0)
c     parameter (b1=0.0705230784d0,b2=0.0422820123d0,b3=0.0092705272d0)
c     parameter (b4=0.0001520143d0,b5=0.0002765672d0,b6=0.0000430638d0)
     
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
         if (.not.((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0))) then
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

*           xerfc = epsilon*ea
*           yerfc=(1.0d0+xerfc*(b1+xerfc*(b2+xerfc*(b3
*    >                 +xerfc*(b4+xerfc*(b5+xerfc*b6))))))**4
*           yerfc = 1.0d0/yerfc**4
*           sum = sum + yerfc/ea

c          x = epsilon*ea
c          yerfc=1.0d0/((1.0d0+x*(7.05230784d-02
c    >      +x*(4.22820123d-02+
c    >      + x*(9.2705272d-03+x*(1.520143d-04
c    >      +x*(2.765672d-04+
c    >      + x*4.30638d-05))))))**16)
c           sum = sum + yerfc/ea

          
         end if
      end do
      end do
      end do
      alpha1 = sum


		
*     **** calculate alpha2 *****
      sum = 0.0d0
      do n1=(-N+1),(N-1) 
      do n2=(-N+1),(N-1) 
      do n3=(-N+1),(N-1) 
         if (.not.((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0))) then
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

      subroutine ewald_end()
      implicit none 

#include "mafdecls.fh"
#include "ewald.fh"

      logical value

      value = MA_free_heap(vg(2))
      value = MA_free_heap(rcell(2))

      return
      end

 
      subroutine ewald_init()
      implicit none 

#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

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

*     **** external functions ****
      integer  control_ncut
      real*8   control_rcut
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega,lattice_unita,psp_zv

      external control_ncut
      external control_rcut
      external ion_nion,ion_katm,G_indx
      external lattice_omega,lattice_unita,psp_zv


*     **** allocate vg memory ****
      call D3dB_nfft3d(1,nfft3d)
      value = MA_alloc_get(mt_dbl,nfft3d,'vg',vg(2),vg(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

*     **** get constants ****
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi

      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      nxh=nx/2

*     ***** find the G==0 index ******
      i=0
      j=0
      k=0
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)
      

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
c      call mask_R(0,dbl_mb(vg(1)))
      call Pack_t_pack(0,dbl_mb(vg(1)))

    
*     **** set the Madelung constant ****
      call mandelung_set(alpha)

      
*     **** ewald summation ****
      rs = (3.0d0*lattice_omega()/fourpi)**(1.0d0/3.0d0)
      zz = 0.0d0
      z  = 0.0d0
      do i=1,ion_nion()
         zz = zz + psp_zv(ion_katm(i))**2
         z  = z  + psp_zv(ion_katm(i))
      end do
c      call D3dB_t_dsum(1,dbl_mb(vg(1)),cewald)
      call Pack_t_dsum(0,dbl_mb(vg(1)),cewald)
      cewald = -0.5d0*zz*(alpha/rs + cewald/lattice_omega())
     >         -0.5d0*(z*z-zz)*rcut*rcut*pi/lattice_omega()
     


*     **** allocate rcell memory ****
      nshl3d=(2*ncut+1)**3
      value = MA_alloc_get(mt_dbl,(3*nshl3d),'rcell',rcell(2),
     >                                           rcell(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      
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
             end if

          end do
        end do
      end do
      

      return
      end

*     ***********************************
*     *			  		*
*     *		ewald_e			*
*     *				       	*
*     ***********************************
      real*8 function ewald_e()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
c     real*8 cerfc,b1,b2,b3,b4,b5,b6
c     parameter (cerfc=1.128379167d0)
c     parameter (b1=0.0705230784d0,b2=0.0422820123d0,b3=0.0092705272d0)
c     parameter (b4=0.0001520143d0,b5=0.0002765672d0,b6=0.0000430638d0)


*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,l
      real*8  w,dx,dy,dz,x,y,z,r,zz
      real*8  yerfc
      real*8  energy,etmp,energyG

*     **** temporary workspace variables ****
c     complex*16  exi(nfft3d)
c     complex*16    s(nfft3d)
c     real*8     tmp3(nfft3d*2)
      integer npack0,nshl3d
      integer exi(2),s(2),tmp3(2),ft(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,ewald_nshl3d
      real*8   lattice_omega,psp_zv,dsum,ion_rion,util_erfc
      external ion_nion,ion_katm,ewald_nshl3d
      external lattice_omega,psp_zv,dsum,ion_rion,util_erfc

      call nwpw_timing_start(21)
      call Parallel_np(np)
      call Parallel_taskid(taskid)

*     **** allocate temp workspace ****
c      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      nshl3d = ewald_nshl3d()
      
      value = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,'s',s(2),s(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, npack0,'tmp3',tmp3(2),tmp3(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (3*nshl3d),'ft',ft(2),ft(1)) 
      if (.not. value) 
     >  call errquit('ewald_e:out of stack memory',0,MA_ERR)
     
*     **** get the structure factor ****
      call dcopy((2*npack0),0.0d0,0,dcpl_mb(s(1)),1)
      do ii=1,ion_nion()
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_cc_daxpy(0,psp_zv(ion_katm(ii)),
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(s(1)))

      end do

      
*     **** calculate the ewald energy ****
      call Pack_ct_Sqr(0,dcpl_mb(s(1)),dbl_mb(tmp3(1)))
      call Pack_tt_dot(0,dbl_mb(tmp3(1)),dbl_mb(vg(1)),energy)
      energyG = 0.5d0*energy/lattice_omega()
      energy  = 0.5d0*energy/lattice_omega() + cewald

      
*     *** need to make parallel ****
      dutask = 0
      etmp = 0.0d0
      do i=1,ion_nion()-1
      do j=i+1,ion_nion()
      if (dutask.eq.taskid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
c        call lattice_min_difference(dx,dy,dz)  
        zz = psp_zv(ion_katm(i)) * psp_zv(ion_katm(j))
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
      if (np.gt.1) call Parallel_SumAll(etmp)
      energy = energy + etmp
      

*     **** deallocate temp workspace ****
      value = MA_pop_stack(ft(2))
      value = value.and.MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(s(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('ewald_e:popping stack memory',0,MA_ERR)

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

#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)
c     real*8 cerfc,b1,b2,b3,b4,b5,b6
c     parameter (b1=0.0705230784d0,b2=0.0422820123d0,b3=0.0092705272d0)
c     parameter (b4=0.0001520143d0,b5=0.0002765672d0,b6=0.0000430638d0)


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
      integer npack0,nshl3d,nion
      integer exi(2),s(2),tmp3(2),ft(2)
      integer fx(2),fy(2),fz(2)
      logical value

*     **** external functions ****
      integer  ion_nion,ion_katm,Pack_G_indx,ewald_nshl3d
      external ion_nion,ion_katm,Pack_G_indx,ewald_nshl3d
      real*8   lattice_omega,psp_zv,dsum,ion_rion,util_erfc
      external lattice_omega,psp_zv,dsum,ion_rion,util_erfc

      call nwpw_timing_start(21)
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nion = ion_nion()

*     **** allocate temp workspace ****
c      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      nshl3d = ewald_nshl3d()
      value = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1)) 
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,'s',s(2),s(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, npack0,'tmp3',tmp3(2),tmp3(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (3*nshl3d),'ft',ft(2),ft(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fx',fx(2),fx(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fy',fy(2),fy(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl, (nion),'fz',fz(2),fz(1)) 
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)


      scal2 = 1.0d0/lattice_omega()
      call dcopy(nion,0.0d0,0,dbl_mb(fx(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fy(1)),1)
      call dcopy(nion,0.0d0,0,dbl_mb(fz(1)),1)
     
*     **** get the structure factor ****
      call dcopy((2*npack0),0.0d0,0,dcpl_mb(s(1)),1)
      do ii=1,nion
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_cc_daxpy(0,psp_zv(ion_katm(ii)),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(s(1)))
      end do

      do ii=1,nion
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         do i=1,npack0
            dbl_mb(tmp3(1)+i-1) 
     >              = ( dble(dcpl_mb(exi(1)+i-1))
     >                *dimag(dcpl_mb(s(1)+i-1))
     >              -  dimag(dcpl_mb(exi(1)+i-1))
     >                 *dble(dcpl_mb(s(1)+i-1))
     >                )*dbl_mb(vg(1)+i-1)
         end do

*        fion(1,ii) = fion(1,ii) + sum*psp_zv(ion_katm(ii))*scal2
*        fion(2,ii) = fion(2,ii) + sum*psp_zv(ion_katm(ii))*scal2
*        fion(3,ii) = fion(3,ii) + sum*psp_zv(ion_katm(ii))*scal2

         call Pack_tt_idot(0,dbl_mb(Pack_G_indx(0,1)),
     >                       dbl_mb(tmp3(1)),sum)
         dbl_mb(fx(1)+ii-1) = dbl_mb(fx(1)+ii-1) 
     >                      +  sum*psp_zv(ion_katm(ii))*scal2

         call Pack_tt_idot(0,dbl_mb(Pack_G_indx(0,2)),
     >                       dbl_mb(tmp3(1)),sum)
         dbl_mb(fy(1)+ii-1) = dbl_mb(fy(1)+ii-1) 
     >                      +  sum*psp_zv(ion_katm(ii))*scal2

         call Pack_tt_idot(0,dbl_mb(Pack_G_indx(0,3)),
     >                       dbl_mb(tmp3(1)),sum)
         dbl_mb(fz(1)+ii-1) = dbl_mb(fz(1)+ii-1) 
     >                      +  sum*psp_zv(ion_katm(ii))*scal2
      end do
    
      

      dutask=0
      do i=1,nion-1
      do j=i+1,nion
        if (dutask.eq.taskid) then
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
c        call lattice_min_difference(dx,dy,dz)  
        zz = psp_zv(ion_katm(i)) * psp_zv(ion_katm(j))
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
      if (np.gt.1) then
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
      value = MA_pop_stack(fz(2))
      value = MA_pop_stack(fy(2))
      value = MA_pop_stack(fx(2))
      value = MA_pop_stack(ft(2))
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(s(2))
      value = MA_pop_stack(exi(2))

      call nwpw_timing_end(21)
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "ewald.fh"

      integer N
      parameter (N=40)

*     **** common block used for coulomb.f ****
*     real*8 vc(nfft3d)
*     common / vc_block / vc
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     ****  expansion coefficient of the error function ****
      real*8 cerfc
      parameter (cerfc=1.128379167d0)
c     real*8 cerfc,b1,b2,b3,b4,b5,b6
c     parameter (b1=0.0705230784d0,b2=0.0422820123d0,b3=0.0092705272d0)
c     parameter (b4=0.0001520143d0,b5=0.0002765672d0,b6=0.0000430638d0)

*     **** local variables ****
      logical value
      integer npack0
      integer i,ii,j,l
      integer n1,n2,n3
      integer u,v,s
      real*8 pi,fourpi,scal
      real*8 zz,z
      real*8 Cus(3,3),hm(3,3),energy,sum,ss,rs
      real*8 ea,ax,ay,az,epsilon
      real*8 dx,dy,dz,w
      real*8 unita(3,3),unitg(3,3)

      integer G(2,3),H(2),F(2),tmp1(2),tmp2(2),exi(2),strf(2)
      integer nshl3d

*     **** external functions ****
      integer  ion_katm,ion_nion,Pack_G_indx,ewald_nshl3d
      real*8   psp_zv,lattice_unitg,lattice_unita,lattice_omega
      real*8   util_erfc,ion_rion
      external ion_katm,ion_nion,Pack_G_indx,ewald_nshl3d
      external psp_zv,lattice_unitg,lattice_unita,lattice_omega
      external util_erfc,ion_rion

      call nwpw_timing_start(21)
      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)

*     *** define hm,unita,unitg ****
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
         unitg(u,v) = lattice_unitg(u,v)
         unita(u,v) = lattice_unita(u,v)
      end do
      end do

c      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      zz = 0.0d0
      z  = 0.0d0
      do i=1,ion_nion()
         zz = zz + psp_zv(ion_katm(i))**2
         z  = z  + psp_zv(ion_katm(i))
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
      value =           MA_push_get(mt_dbl,npack0,'H',H(2),H(1))
      value = value.and.MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      value = value.and.
     >      MA_push_get(mt_dcpl,npack0,'strf',strf(2),strf(1))
      if (.not. value) call errquit('out of stack memory',0,MA_ERR)

      call dcopy((2*npack0),0.0d0,0,dcpl_mb(strf(1)),1)
      do ii=1,ion_nion()
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
c         ***!call Pack_c_pack(0,dcpl_mb(exi(1)))
         call Pack_cc_daxpy(0,psp_zv(ion_katm(ii)),dcpl_mb(exi(1)),
     >                                             dcpl_mb(strf(1)))
      end do
      call Pack_ct_Sqr(0,dcpl_mb(strf(1)),dbl_mb(H(1)))
      value =           MA_pop_stack(strf(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0,MA_ERR)

*     **** calculate the ewald energy ****
c      value = MA_push_get(mt_dbl,nfft3d,'F',F(2),F(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c      call dcopy(nfft3d,dbl_mb(vg(1)),    1,dbl_mb(F(1)),  1)
c      call Pack_t_Pack(0,dbl_mb(F(1)))
      F(1) = vg(1)

      call Pack_tt_dot(0,dbl_mb(F(1)),dbl_mb(H(1)),energy)
      energy = -0.5d0*energy/lattice_omega()


      do v=1,3
      do u=1,3
         stress(u,v) = stress(u,v) + energy*hm(u,v)
      end do
      end do
     
*     **** tmp2(G) = F(G)*H(G)/G**2 + F(G)*H(G)*rcut*rcut/4 ****
      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call Pack_tt_Mul(0,dbl_mb(F(1)),
     >                   dbl_mb(H(1)),
     >                   dbl_mb(tmp1(1)))
      ss = 0.25d0*rcut*rcut
      call Pack_t_SMul(0,ss,dbl_mb(tmp1(1)),
     >                      dbl_mb(tmp2(1)))
      ss = 1.0d0/fourpi
c      call Pack_t_SMul(0,ss,dbl_mb(tmp1(1)),
c     >                      dbl_mb(tmp1(1)))
c      call Pack_tt_Mul(0,dbl_mb(tmp1(1)),
c     >                   dbl_mb(vc_indx),
c     >                   dbl_mb(tmp1(1)))
c      call Pack_tt_Sum(0,dbl_mb(tmp1(1)),
c     >                   dbl_mb(tmp2(1)),
c     >                   dbl_mb(tmp2(1)))
      call Pack_t_SMul1(0,ss,dbl_mb(tmp1(1)))
      call Pack_tt_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp1(1)))
      call Pack_tt_Sum2(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)))


*     **** calculate Cus ****
c      value =           MA_push_get(mt_dbl,nfft3d,
c     >                     'G1',G(2,1),G(1,1))
c      value = value.and.MA_push_get(mt_dbl,nfft3d,
c     >                     'G2',G(2,2),G(1,2))
c      value = value.and.MA_push_get(mt_dbl,nfft3d,
c     >                     'G3',G(2,3),G(1,3))
c      if (.not. value) call errquit('out of stack  memory',0, MA_ERR)
c      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(G(1,1)),1)
c      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(G(1,2)),1)
c      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(G(1,3)),1)
c      call Pack_t_pack(0,dbl_mb(G(1,1)))
c      call Pack_t_pack(0,dbl_mb(G(1,2)))
c      call Pack_t_pack(0,dbl_mb(G(1,3)))

      call dcopy(9,0.0d0,0,Cus,1)
c     ss = -1.0d0/lattice_omega()
      ss =  1.0d0/lattice_omega()
      do u=1,3
      do s=u,3
         call Pack_tt_Mul(0,dbl_mb(Pack_G_indx(0,u)),
     >                      dbl_mb(Pack_G_indx(0,s)),
     >                      dbl_mb(tmp1(1)))
         call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)
         Cus(u,s) = ss*sum
      end do
      end do
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

c      value =           MA_pop_stack(G(2,3))
c      value = value.and.MA_pop_stack(G(2,2))
c      value = value.and.MA_pop_stack(G(2,1))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
c      value = value.and.MA_pop_stack(F(2))
      value = value.and.MA_Pop_stack(H(2))
      if (.not. value) call errquit('error popping stack memory',0,
     &       MA_ERR)


*     **** R-space contributions ****


*     **** calculate alpha1 - stress from cewald term*****
        call dcopy(9,0.0d0,0,Cus,1)
        rs      = (3.0d0*lattice_omega()/(4.0d0*pi))**(1.0d0/3.0d0)
        epsilon = 1.0d0/rcut
        sum = 0.0d0
        do n1=(-N+1),(N-1) 
        do n2=(-N+1),(N-1) 
        do n3=(-N+1),(N-1) 
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

              Cus(2,1) = Cus(2,1) + ss * ay*ax 
              Cus(2,2) = Cus(2,2) + ss * ay*ay 
              Cus(2,3) = Cus(2,3) + ss * ay*az 

              Cus(3,1) = Cus(3,1) + ss * az*ax 
              Cus(3,2) = Cus(3,2) + ss * az*ay 
              Cus(3,3) = Cus(3,3) + ss * az*az 

           end if
        end do
        end do
        end do

c       do u=1,3
c       do s=u+1,3
c          Cus(s,u) = Cus(u,s)
c       end do
c       end do

        do v=1,3
        do u=1,3
          do s=1,3
             stress(u,v) = stress(u,v) + Cus(u,s)*hm(s,v)
          end do
        end do
        end do

*     *** need to make parallel ****
*     **** calculate erfc contribution *****
      nshl3d = ewald_nshl3d()
      call dcopy(9,0.0d0,0,Cus,1)
      epsilon = 1.0d0/rcut
      do i=1,ion_nion()-1
      do j=i+1,ion_nion()
        dx = ion_rion(1,i) - ion_rion(1,j)
        dy = ion_rion(2,i) - ion_rion(2,j)
        dz = ion_rion(3,i) - ion_rion(3,j)
c        call lattice_min_difference(dx,dy,dz)  
        zz = psp_zv(ion_katm(i)) * psp_zv(ion_katm(j))
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
      end do
      end do
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
