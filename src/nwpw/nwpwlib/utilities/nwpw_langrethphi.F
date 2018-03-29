*     ************************************************
*     *                                              *
*     *             nwpw_langrethphi_init            *
*     *                                              *
*     ************************************************
*     
      subroutine nwpw_langrethphi_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

c     **** Langreth kernel common block ****
      integer Na
      parameter (Na=256)
      real*8 amin,amax
      integer a(2),a2(2),aweights(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2),Wab(2)
      common /langrethphi_common/ Wab,a,a2,aweights,
     >                            cos_a,sin_a,nu,nu1,amin,amax

*     ***** local variables ****
      logical value

      value = BA_alloc_get(mt_dbl,Na*Na,'lg_Wab',Wab(2),Wab(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_a',a(2),a(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_a2',a2(2),a2(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_aw',aweights(2),aweights(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_cs',cos_a(2),cos_a(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_sn',sin_a(2),sin_a(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_nu',nu(2),nu(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,Na,'lg_nu1',nu1(2),nu1(1))
      if (.not.value)
     >   call errquit('nwpw_langrethphi_init: alloc heap',0,MA_ERR)



      call nwpw_langrethphi_setup(Na,
     >                         dbl_mb(Wab(1)),
     >                         dbl_mb(a(1)),
     >                         dbl_mb(a2(1)),
     >                         dbl_mb(aweights(1)),
     >                         dbl_mb(cos_a(1)),
     >                         dbl_mb(sin_a(1)),amin,amax)
      return
      end

*     ************************************************
*     *                                              *
*     *             nwpw_langrethphi_setup           *
*     *                                              *
*     ************************************************
*     
      subroutine nwpw_langrethphi_setup(Na,Wab,a,a2,aweights,
     >                                  cos_a,sin_a,amin,amax)
      implicit none
      real*8 Wab(Na,Na)
      real*8 a(Na),a2(Na),aweights(Na)
      real*8 cos_a(Na),sin_a(Na)
      real*8 amin,amax

*     **** local variables ****
      integer i,j

      call nwpw_langrethph_GL(datan(amin),datan(amax),Na,a,aweights)
      do i=1,Na
         a(i) = dtan(a(i))
         a2(i) = a(i)**2
         aweights(i) = aweights(i)*(1.0d0+a2(i))
         cos_a(i) = dcos(a(i))
         sin_a(i) = dsin(a(i))
      end do
      do i=1,Na
      do j=1,Na
         Wab(i,j) = 2.0d0 * aweights(i)*aweights(j) *  (
     >       (3.0d0-a2(i))*a(j)*cos_a(j)*sin_a(i)  +
     >       (3.0d0-a2(j))*a(i)*cos_a(i)*sin_a(j)  +
     >       (a2(i)+a2(j)-3.0d0)*sin_a(i)*sin_a(j) -
     >       3.0d0*a(i)*a(j)*cos_a(i)*cos_a(j) )   / (a(i)*a(j))
      end do
      end do


      return
      end 

*     ************************************************
*     *                                              *
*     *             nwpw_langrethphi_end             *
*     *                                              *
*     ************************************************
*     

      subroutine nwpw_langrethphi_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

c     **** Langreth kernel common block ****
      integer Na
      parameter (Na=256)
      real*8 amin,amax
      integer a(2),a2(2),aweights(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2),Wab(2)
      common /langrethphi_common/ Wab,a,a2,aweights,
     >                            cos_a,sin_a,nu,nu1,amin,amax

c     **** local variables ****
      logical value

      value = BA_free_heap(Wab(2))
      value = value.and.BA_free_heap(a(2))
      value = value.and.BA_free_heap(a2(2))
      value = value.and.BA_free_heap(aweights(2))
      value = value.and.BA_free_heap(cos_a(2))
      value = value.and.BA_free_heap(sin_a(2))
      value = value.and.BA_free_heap(nu(2))
      value = value.and.BA_free_heap(nu1(2))
      if (.not.value)
     >   call errquit('nwpw_langrethphi_end: free heap',0,MA_ERR)
      return
      end

*     ************************************************
*     *                                              *
*     *           nwpw_langrethphi                   *
*     *                                              *
*     ************************************************
*     
      real*8 function nwpw_langrethphi(d1,d2)
      implicit none
      real*8 d1,d2

#include "bafdecls.fh"
#include "errquit.fh"

c     **** Langreth kernel common block ****
      integer Na
      parameter (Na=256)
      real*8 amin,amax
      integer a(2),a2(2),aweights(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2),Wab(2)
      common /langrethphi_common/ Wab,a,a2,aweights,
     >                            cos_a,sin_a,nu,nu1,amin,amax

*     **** external functions ****
      real*8   nwpw_langrethphi_tmp
      external nwpw_langrethphi_tmp

      nwpw_langrethphi = nwpw_langrethphi_tmp(Na,dbl_mb(Wab(1)),
     >                                        dbl_mb(a2(1)),
     >                                        dbl_mb(nu(1)),
     >                                        dbl_mb(nu1(1)),
     >                                        d1,d2)
      return
      end

*     ************************************************
*     *                                              *
*     *           nwpw_langrethphi_tmp               *
*     *                                              *
*     ************************************************
*     
      real*8 function nwpw_langrethphi_tmp(Na,Wab,a2,nu,nu1,d1,d2)
      implicit none
      integer Na
      real*8  Wab(Na,Na),a2(Na),nu(Na),nu1(Na)
      real*8  d1,d2

*     **** local variables ****
      real*8 small
      parameter (small=1.0d-12)
      integer i,j
      real*8  w,x,y,z,T,pi,gamma,d1s,d2s,phi_value

      pi    = 4.0d0*datan(1.0d0)
      gamma = 4.0d0*pi/9.0d0
      d1s = d1**2
      d2s = d2**2

      phi_value = 0.0d0
      if ((d1.gt.small).or.(d2.gt.small)) then
         do i=1,Na
            if ((a(i).le.small).and.(d1.gt.small)) then
               nu(i) = 9.0d0/8.0d0*d1s/pi
            else if ((d1.le.small).or.
     >               (((a2(i)*gamma)/d1s).gt.700.0d0)) then
               nu(i) = a2(i)/2.0d0
            else
               nu(i) = a2(i)
     >                /((1.0d0-dexp(-(a2(i)*gamma)/d1s))*2.0d0)
            end if

            if ((a(i).le.small).and.(d2.gt.small)) then
               nu1(i) = 9.0d0/8.0d0*d2s/pi
            else if ((d2.le.small).or.
     >               (((a2(i)*gamma)/d2s).gt.700.0d0)) then
               nu1(i) = a2(i)/2.0d0
            else
               nu1(i) = a2(i)
     >                /((1.0d0-dexp(-(a2(i)*gamma)/d2s))*2.0d0)
            end if
         end do
         do i=1,Na
         do j=1,Na
             w = nu(i)
             x = nu(j)
             y = nu1(i)
             z = nu1(j)
             T = (1.0d0/(w+x) + 1.0d0/(y+z))
     >          *(1.0d0/((w+y)*(x+z)) + 1.0d0/((w+z)*(y+x)))
             phi_value = phi_value + T * Wab(i,j)
         end do
         end do
         phi_value = phi_value/(pi**2)
      end if

      nwpw_langrethphi_tmp = phi_value
      return
      end



*     ************************************************
*     *                                              *
*     *           nwpw_langrethphi_GL                *
*     *                                              *
*     ************************************************
*     
      subroutine nwpw_langrethphi_GL(amin,amax,Npoints,a,weights)
      implicit none
      real*8 amin,amax
      integer Npoints
      real*8 a(*),weights(*)

c     **** local varaiables ****    
      logical done
      real*8  poly1,poly2,poly3   !! Temporary storage for Legendre Polynomials
      integer i,j,N               !! Indexing variables 
      real*8  root,dpdx,last_root !! The value of the root of a given Legendre polynomial,
      real*8  midpoint,length     !! The middle of the x-range and the length to that point
      real*8  pi

      pi = 4.0d0*datan(1.0d0)
      N = (Npoints + 1)/2
      midpoint = 0.5d0*(amin+amax)
      length   = 0.5d0*(amax-amin)

      do i=1,N
         root = dcos(dble(pi*(i-0.25d0)/(Npoints+0.5d0)))

         done = .false.
         do while (.not.done)
            poly1 = 1.0d0
            poly2 = 0.0d0
            do j = 1, Npoints
               poly3 = poly2
               poly2 = poly1
               poly1 = ((2.0d0*j-1.0d0)*root*poly2-(j-1.0d0)
     >                  *poly3)/j
            end do
            dpdx=Npoints*(root*poly1-poly2)/(root**2-1.0d0)

            last_root = root
            root = last_root - poly1/dpdx

            if (dabs(root-last_root).le.1.0d-14) done = .true.
         end do

         a(i)           = midpoint-length*root
         a(Npoints+1-i) = midpoint+length*root

         weights(i) = 2.0d0*length/((1.0d0-root**2)*dpdx**2)
         weights(Npoints+1-i) = weights(i)
      end do

      return
      end





