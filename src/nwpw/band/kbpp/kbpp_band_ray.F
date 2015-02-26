*
* $Id$
*

*     *********************************************
*     *                                           *
*     *               kbpp_band_calc_Gmax         *
*     *                                           *
*     *********************************************

      double precision function kbpp_band_calc_Gmax(nfft1,nfft2,nfft3,G)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)

*     *** local variables ****
      integer k1,k2,k3
      double precision q,qmax

      qmax = 0.0d0
      do k3=1,nfft3
      do k2=1,nfft2
      do k1=1,nfft1
        q=DSQRT(G(k1,k2,k3,1)**2
     >         +G(k1,k2,k3,2)**2
     >         +G(k1,k2,k3,3)**2)
        if (q.gt.qmax) qmax = q
      end do
      end do
      end do

*     *** make Gmax = Gmax + 2.0d0 ****
      qmax = qmax + 2.0d0

      kbpp_band_calc_Gmax = qmax
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_band_calc_dGmin        *
*     *                                           *
*     *********************************************

      double precision function kbpp_band_calc_dGmin(nfft1,nfft2,nfft3,
     >                                               G)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)

*     *** local variables ****
      double precision gx,gy,gz,q,dqmin

      gx=G(2,1,1,1) - G(1,1,1,1)
      gy=G(2,1,1,2) - G(1,1,1,2)
      gz=G(2,1,1,3) - G(1,1,1,3)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      dqmin = q

      gx=G(1,2,1,1) - G(1,1,1,1)
      gy=G(1,2,1,2) - G(1,1,1,2)
      gz=G(1,2,1,3) - G(1,1,1,3)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      if (q.lt.dqmin) dqmin = q

      gx=G(1,1,2,1) - G(1,1,1,1)
      gy=G(1,1,2,2) - G(1,1,1,2)
      gz=G(1,1,2,3) - G(1,1,1,3)
      q = dsqrt(gx**2 + gy**2 + gz**2)
      if (q.lt.dqmin) dqmin = q

*     *** make fine dqmin ****
      dqmin = 0.01d0*dqmin

      kbpp_band_calc_dGmin = dqmin
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_band_calc_nray         *
*     *                                           *
*     *********************************************
      integer function kbpp_band_calc_nray(nfft1,nfft2,nfft3,G)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)

*     ***** local variables ****
      double precision qmax
      double precision dG,Gmax
      integer nray

*     ***** external functions ****
      real*8   kbpp_band_calc_dGmin,kbpp_band_calc_Gmax
      external kbpp_band_calc_dGmin,kbpp_band_calc_Gmax

      dG   = kbpp_band_calc_dGmin(nfft1,nfft2,nfft3,G)
      Gmax = kbpp_band_calc_Gmax(nfft1,nfft2,nfft3,G)
      nray = Gmax/dG + 1.0d0
      if (nray.lt.10) nray = 10

      kbpp_band_calc_nray = nray
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_band_generate_G_ray    *
*     *                                           *
*     *********************************************
      subroutine kbpp_band_generate_G_ray(nfft1,nfft2,nfft3,G,G_ray)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision G(nfft1,nfft2,nfft3,3)
      double precision G_ray(*)

*     **** local variables ***
      integer i,nray
      double precision dG
*     ***** external functions ****
      real*8   kbpp_band_calc_dGmin,kbpp_band_calc_Gmax
      integer  kbpp_band_calc_nray
      external kbpp_band_calc_dGmin,kbpp_band_calc_Gmax
      external kbpp_band_calc_nray

      dG   = kbpp_band_calc_dGmin(nfft1,nfft2,nfft3,G)
      nray = kbpp_band_calc_nray(nfft1,nfft2,nfft3,G)

      do i=1,nray
       G_ray(i) = dG*dble(i-1)
      end do
      return
      end


*     *********************************************
*     *                                           *
*     *               kbpp_band_filter_ray        *
*     *                                           *
*     *********************************************
      subroutine kbpp_band_filter_ray(nray,G_ray,ecut,V_ray)
      implicit none
      integer nray
      double precision G_ray(*)
      double precision ecut
      double precision V_ray(*)

*     **** local variables ****
      integer ncut
      parameter (ncut=15)
      integer i
      double precision g,qmax,fac

      qmax = dsqrt(ecut+ecut)
      do i=1,nray
         g = G_ray(i)
         if (g.gt.qmax) then
            fac = 1.0d0 - (1.0d0-exp(-(g/qmax)**ncut))**ncut
            V_ray(i) = V_ray(i)*fac
         end if
      end do
      return
      end

*     *********************************************
*     *                                           *
*     *               kbpp_band_filter            *
*     *                                           *
*     *********************************************
      subroutine kbpp_band_filter(nfft1,nfft2,nfft3,kvec,G,ecut,vl)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision kvec(3)
      double precision G(nfft1,nfft2,nfft3,3)
      double precision ecut
      double precision vl(nfft1,nfft2,nfft3)

*     **** local variables ****
      !real*8 ncut,eps
      !parameter (eps=1.0d-12,ncut=15.0d0)
      integer ncut
      parameter (ncut=15)
      integer i,j,k
      double precision q,qmax,fac
      double precision gx,gy,gz

c      qmax = dsqrt(ecut+ecut)
c     >      /(-log(1.0d0-(1.0d0-eps)**(1.0d0/ncut)))**(1.0d0/ncut)
      qmax = dsqrt(ecut+ecut)
      do k=1,nfft3
      do j=1,nfft2
      do i=1,nfft1
         gx = G(i,j,k,1) + kvec(1)
         gy = G(i,j,k,2) + kvec(2)
         gz = G(i,j,k,3) + kvec(3)
         q = dsqrt(gx*gx + gy*gy + gz*gz)

         if (q.gt.(qmax-0.2d0)) then
            fac = 1.0d0 - (1.0d0-exp(-(q/qmax)**ncut))**ncut
            vl(i,j,k) = vl(i,j,k)*fac
         end if
      end do
      end do
      end do
      return
      end



