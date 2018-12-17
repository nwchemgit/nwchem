*
* $Id$
*

*     *********************************************
*     *                                           *
*     *               kbpp_calc_Gmax              *
*     *                                           *
*     *********************************************

      double precision function kbpp_calc_Gmax(nfft1,nfft2,nfft3,unita)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3)

*     *** local variables ****
      integer k1,k2,k3
      double precision qmax,gx,gy,gz,gg,ggmax
      double precision unitg(3,3)

      call get_unitg(unita,unitg)
      k1 = nfft1/2 
      k2 = nfft2/2 
      k3 = nfft3/2 
      ggmax = 0.0d0
      gx = unitg(1,1)*k1
      gy = unitg(2,1)*k1
      gz = unitg(3,1)*k1
      gg = gx*gx + gy*gy + gz*gz
      if (gg.gt.ggmax) ggmax = gg
      gx = unitg(1,2)*k2
      gy = unitg(2,2)*k2
      gz = unitg(3,2)*k2
      gg = gx*gx + gy*gy + gz*gz
      if (gg.gt.ggmax) ggmax = gg
      gx = unitg(1,1)*k3
      gy = unitg(2,1)*k3
      gz = unitg(3,1)*k3
      gg = gx*gx + gy*gy + gz*gz
      if (gg.gt.ggmax) ggmax = gg

      qmax = dsqrt(ggmax)

      qmax = qmax 

      kbpp_calc_Gmax = qmax
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_calc_dGmin             *
*     *                                           *
*     *********************************************

      double precision function kbpp_calc_dGmin(unita)
      implicit none
      double precision unita(3,3)

*     *** local variables ****
      double precision gx,gy,gz,q,dqmin
      double precision unitg(3,3)

      call get_unitg(unita,unitg)

      gx = unitg(1,1) 
      gy = unitg(2,1)
      gz = unitg(3,1)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      dqmin = q

      gx = unitg(1,2)
      gy = unitg(2,2)
      gz = unitg(3,2)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      dqmin = q

      gx = unitg(1,3)
      gy = unitg(2,3)
      gz = unitg(3,3)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      dqmin = q

*     *** make fine dqmin ****
      dqmin = 0.01d0*dqmin

      kbpp_calc_dGmin = dqmin
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_calc_nray              *
*     *                                           *
*     *********************************************
      integer function kbpp_calc_nray(nfft1,nfft2,nfft3,unita)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3)

*     ***** local variables ****
      double precision dG,Gmax
      integer nray

*     ***** external functions ****
      real*8   kbpp_calc_dGmin,kbpp_calc_Gmax
      external kbpp_calc_dGmin,kbpp_calc_Gmax

      dG   = kbpp_calc_dGmin(unita)
      Gmax = kbpp_calc_Gmax(nfft1,nfft2,nfft3,unita) + 2.0d0
      nray = Gmax/dG + 1.0d0
      if (nray.lt.10) nray = 10

      kbpp_calc_nray = nray
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_generate_G_ray         *
*     *                                           *
*     *********************************************
      subroutine kbpp_generate_G_ray(nfft1,nfft2,nfft3,unita,G_ray)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3)
      double precision G_ray(*)

*     **** local variables ***
      integer i,nray
      double precision dG
*     ***** external functions ****
      real*8   kbpp_calc_dGmin
      integer  kbpp_calc_nray
      external kbpp_calc_dGmin
      external kbpp_calc_nray

      dG   = kbpp_calc_dGmin(unita)
      nray = kbpp_calc_nray(nfft1,nfft2,nfft3,unita)
      do i=1,nray
       G_ray(i) = dG*dble(i-1)
      end do
      return
      end


*     *********************************************
*     *                                           *
*     *               kbpp_filter_ray             *
*     *                                           *
*     *********************************************
      subroutine kbpp_filter_ray(nray,G_ray,ecut,V_ray)
      implicit none
      integer nray
      double precision G_ray(*)
      double precision ecut
      double precision V_ray(*)

*     **** local variables ****
      !real*8 ncut,eps
      !parameter (eps=1.0d-12,ncut=15.0d0)
      integer ncut
      parameter (ncut=15)
      integer i
      double precision g,qmax,fac

c      qmax = dsqrt(ecut+ecut)
c     >      /(-log(1.0d0-(1.0d0-eps)**(1.0d0/ncut)))**(1.0d0/ncut)
      qmax = dsqrt(ecut+ecut) 
      do i=1,nray
         g = G_ray(i)
         if (g.gt.(qmax-0.2d0)) then
            fac = 1.0d0 - (1.0d0-exp(-(g/qmax)**ncut))**ncut
            V_ray(i) = V_ray(i)*fac
         end if
      end do
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_filter                 *
*     *                                           *
*     *********************************************
      subroutine kbpp_filter(nfft1,nfft2,nfft3,G,ecut,vl)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision G(nfft1/2+1,nfft2,nfft3,3)
      double precision ecut
      double precision vl(nfft1/2+1,nfft2,nfft3)

*     **** local variables ****
      !real*8 ncut,eps
      !parameter (eps=1.0d-12,ncut=15.0d0)
      integer ncut
      parameter (ncut=15)
      integer i,j,k
      double precision q,qmax,fac

c      qmax = dsqrt(ecut+ecut)
c     >      /(-log(1.0d0-(1.0d0-eps)**(1.0d0/ncut)))**(1.0d0/ncut)
      qmax = dsqrt(ecut+ecut)
      do k=1,nfft3
      do j=1,nfft2
      do i=1,(nfft1/2+1)
         q=DSQRT(G(i,j,k,1)**2
     >          +G(i,j,k,2)**2
     >          +G(i,j,k,3)**2)

         if (q.gt.(qmax-0.2d0)) then
            fac = 1.0d0 - (1.0d0-exp(-(q/qmax)**ncut))**ncut
            vl(i,j,k) = vl(i,j,k)*fac
         end if
      end do
      end do
      end do
      return
      end
