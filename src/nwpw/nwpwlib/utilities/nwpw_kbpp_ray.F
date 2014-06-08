
*     *********************************************
*     *                                           *
*     *          nwpw_kbpp_calc_Gmax              *
*     *                                           *
*     *********************************************

      double precision function nwpw_kbpp_calc_Gmax(nfft1,nfft2,nfft3,
     >                                              unita)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3)

*     *** local variables ****
      integer k1,k2,k3
      double precision qmax,gx,gy,gz
      double precision unitg(3,3)

      call nwpw_kbpp_get_unitg(unita,unitg)
      k1 = nfft1/2
      k2 = nfft2/2
      k3 = nfft3/2
      gx = k1*unitg(1,1) + k2*unitg(1,2) +k3*unitg(1,3)
      gy = k1*unitg(2,1) + k2*unitg(2,2) +k3*unitg(2,3)
      gz = k1*unitg(3,1) + k2*unitg(3,2) +k3*unitg(3,3)

      qmax = dsqrt(gx**2 + gy**2 + gz**2)

*     *** make Gmax = Gmax + 2.0d0 ****
      qmax = qmax + 2.0d0

      nwpw_kbpp_calc_Gmax = qmax
      return
      end


*     *********************************************
*     *                                           *
*     *            nwpw_kbpp_calc_dGmin           *
*     *                                           *
*     *********************************************

      double precision function nwpw_kbpp_calc_dGmin(unita)
      implicit none
      double precision unita(3,3)

*     *** local variables ****
      double precision gx,gy,gz,q,dqmin
      double precision unitg(3,3)

      call nwpw_kbpp_get_unitg(unita,unitg)

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

      nwpw_kbpp_calc_dGmin = dqmin
      return
      end


*     *********************************************
*     *                                           *
*     *             nwpw_kbpp_calc_nray           *
*     *                                           *
*     *********************************************
      integer function nwpw_kbpp_calc_nray(nfft1,nfft2,nfft3,unita)
      implicit none
      integer nfft1,nfft2,nfft3
      double precision unita(3,3)

*     ***** local variables ****
      double precision dG,Gmax
      integer nray

*     ***** external functions ****
      real*8   nwpw_kbpp_calc_dGmin,nwpw_kbpp_calc_Gmax
      external nwpw_kbpp_calc_dGmin,nwpw_kbpp_calc_Gmax

      dG   = nwpw_kbpp_calc_dGmin(unita)
      Gmax = nwpw_kbpp_calc_Gmax(nfft1,nfft2,nfft3,unita)
      nray = Gmax/dG + 1.0d0
      if (nray.lt.10) nray = 10

      nwpw_kbpp_calc_nray = nray
      return
      end


*     *********************************************
*     *                                           *
*     *           nwpw_kbpp_generate_G_ray        *
*     *                                           *
*     *********************************************
      subroutine nwpw_kbpp_generate_G_ray(nfft1,nfft2,nfft3,unita,G_ray)
      implicit none
      integer nfft1,nfft2,nfft3 
      double precision unita(3,3)
      double precision G_ray(*)

*     **** local variables ***
      integer i,nray
      double precision dG
*     ***** external functions ****
      real*8   nwpw_kbpp_calc_dGmin,nwpw_kbpp_calc_Gmax
      integer  nwpw_kbpp_calc_nray
      external nwpw_kbpp_calc_dGmin,nwpw_kbpp_calc_Gmax
      external nwpw_kbpp_calc_nray

      dG   = nwpw_kbpp_calc_dGmin(unita)
      nray = nwpw_kbpp_calc_nray(nfft1,nfft2,nfft3,unita)
      do i=1,nray
       G_ray(i) = dG*dble(i-1)
      end do
      return
      end


*     *********************************************
*     *                                           *
*     *             nwpw_kbpp_filter_ray          *
*     *                                           *
*     *********************************************
      subroutine nwpw_kbpp_filter_ray(nray,G_ray,ecut,V_ray)
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


      subroutine nwpw_kbpp_get_unitg(unita,unitg)
      implicit none

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors                                *
*               in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             unita  --- primitive vectors in coordination space             *
*                                                                            *
*     Outputs:                                                               *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  dscal from BLAS                                              *
*                                                                            *
*     Last modification:  3/30/99  by Eric Bylaska                           *
*                                                                            *
******************************************************************************


*     ------------------
*     argument variables
*     ------------------
      double precision unita(3,3), unitg(3,3)

*     ---------------
*     local variables
*     ---------------
      double precision volume
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)

*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume=dabs(volume)

      return
      end





c $Id$
