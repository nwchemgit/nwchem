c $Id$

*     ***********************************
*     *                                 *
*     *          nwpw_gaunt_init        *
*     *                                 *
*     ***********************************
      subroutine nwpw_gaunt_init(iscmplx,lmax0)
      implicit none
      logical iscmplx
      integer lmax0

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_gaunt common block ****
      integer nwpw_gaunt_lmax(2),
     >        nwpw_gaunt_sizel2(2),
     >        nwpw_gaunt_coeff(2,2)
      logical nwpw_gaunt_iscmplx
      common / nwpw_gaunt_common / nwpw_gaunt_coeff,
     >                             nwpw_gaunt_lmax,
     >                             nwpw_gaunt_sizel2,
     >                             nwpw_gaunt_iscmplx

*     **** local variables ****
      integer sizel,i,j,k,indx,ic
      integer l1,l2,l3
      integer m1,m2,m3

*     **** external functions ****
      real*8   nwpw_gaunt_sub
      external nwpw_gaunt_sub

      nwpw_gaunt_iscmplx = iscmplx
      ic = 1
      if (iscmplx) ic = 2

      nwpw_gaunt_lmax(ic)    = lmax0
      sizel                  = (nwpw_gaunt_lmax(ic)+1)
      nwpw_gaunt_sizel2(ic)  = sizel**2

      if (.not.BA_alloc_get(mt_dbl,
     >nwpw_gaunt_sizel2(ic)*nwpw_gaunt_sizel2(ic)*nwpw_gaunt_sizel2(ic),
     >                    'nwpw_gaunt_coeff',
     >                     nwpw_gaunt_coeff(ic,2),
     >                     nwpw_gaunt_coeff(ic,1)))
     >   call errquit('nwpw_gaunt_init:out of heap memory',0,MA_ERR)

c      write(*,*) "nwpw_gaunt_lmax   =",nwpw_gaunt_lmax
c      write(*,*) "nwpw_gaunt_iscmplx=",nwpw_gaunt_iscmplx
c
      do l1=0,nwpw_gaunt_lmax(ic)
      do l2=0,nwpw_gaunt_lmax(ic)
      do l3=0,nwpw_gaunt_lmax(ic)

        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
           i = l1**2 + (l1+m1)
           j = l2**2 + (l2+m2)
           k = l3**2 + (l3+m3)
           indx = i
     >          + j*nwpw_gaunt_sizel2(ic) 
     >          + k*nwpw_gaunt_sizel2(ic)*nwpw_gaunt_sizel2(ic)

           dbl_mb(nwpw_gaunt_coeff(ic,1)+indx) 
     >     = nwpw_gaunt_sub(iscmplx,l1,m1,l2,m2,l3,m3)

        end do
        end do
        end do

      end do
      end do
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          nwpw_gaunt_end         *
*     *                                 *
*     ***********************************
      subroutine nwpw_gaunt_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_gaunt common block ****
      integer nwpw_gaunt_lmax(2),
     >        nwpw_gaunt_sizel2(2),
     >        nwpw_gaunt_coeff(2,2)
      logical nwpw_gaunt_iscmplx
      common / nwpw_gaunt_common / nwpw_gaunt_coeff,
     >                             nwpw_gaunt_lmax,
     >                             nwpw_gaunt_sizel2,
     >                             nwpw_gaunt_iscmplx

      if (.not.BA_free_heap(nwpw_gaunt_coeff(1,2)))
     >   call errquit('nwpw_gaunt_end:freeing heap memory',0,MA_ERR)

      if (nwpw_gaunt_iscmplx) then
         if (.not.BA_free_heap(nwpw_gaunt_coeff(2,2)))
     >      call errquit('nwpw_gaunt_end:freeing heap memory',0,MA_ERR)
      end if

      return
      end

*     ***********************************
*     *                                 *
*     *          nwpw_gaunt             *
*     *                                 *
*     ***********************************
      real*8 function nwpw_gaunt(iscmplx,
     >                           l1,m1,l2,m2,l3,m3)
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_gaunt common block ****
      integer nwpw_gaunt_lmax(2),
     >        nwpw_gaunt_sizel2(2),
     >        nwpw_gaunt_coeff(2,2)
      logical nwpw_gaunt_iscmplx
      common / nwpw_gaunt_common / nwpw_gaunt_coeff,
     >                             nwpw_gaunt_lmax,
     >                             nwpw_gaunt_sizel2,
     >                             nwpw_gaunt_iscmplx

*     **** local variables ****
      integer i,j,k,indx,ic
      real*8  gg

*     **** external functions ****
      real*8   nwpw_gaunt_sub
      external nwpw_gaunt_sub

      !*** initialize complex gaunt coefficients if called ***
      if ((iscmplx).and.(.not.nwpw_gaunt_iscmplx)) then
         call nwpw_gaunt_init(iscmplx,nwpw_gaunt_lmax(1))
      end if

      ic = 1
      if (iscmplx) ic = 2

      call nwpw_timing_start(51)
      if ((l1.gt.nwpw_gaunt_lmax(ic)).or.
     >    (l2.gt.nwpw_gaunt_lmax(ic)).or.
     >    (l3.gt.nwpw_gaunt_lmax(ic))) then
        !write(*,*) "calling gaunt_sub=",l1,m1,l2,m2,l3,m3
        gg = nwpw_gaunt_sub(iscmplx,l1,m1,l2,m2,l3,m3)
      else
         i = l1**2 + (l1+m1)
         j = l2**2 + (l2+m2)
         k = l3**2 + (l3+m3)
         indx = i
     >        + j*nwpw_gaunt_sizel2(ic)
     >        + k*nwpw_gaunt_sizel2(ic)*nwpw_gaunt_sizel2(ic)
         gg = dbl_mb(nwpw_gaunt_coeff(ic,1)+indx)

      end if
      call nwpw_timing_end(51)
      nwpw_gaunt = gg

      return
      end
c

*     ***********************************
*     *                                 *
*     *          nwpw_gaunt_sub         *
*     *                                 *
*     ***********************************
*
* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi 
*     |     |
*   = |     | dconjg(Y_lm(theta,phi)) * Y_l1m1(theta,phi) * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | Y_lm(theta,phi) * dconjg(Y_l1m1(theta,phi)) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*   = sqrt( (2*l+1)*(2*l2+1)/(4*pi*(2*l1+1)) ) * Clebsch(l l2 l1; m m2 m1) * Clebsch(l l2 l1; 0 0 0)
*
*
* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | T_lm(theta,phi)) * T_l1m1(theta,phi) * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
      real*8 function nwpw_gaunt_sub(iscmplx,
     >                           l1,m1,l2,m2,l3,m3)
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     
*     *** local variables ****
      double precision pi
      double precision piover2,twopi,fourpi
c      parameter(pi = 3.14159265358979323846264338327950288419d0)
c      parameter(piover2= pi*0.5d0)
c      parameter(twopi  = pi*2.0d0)
c      parameter(fourpi = pi*4.0d0)

c     work arrays for integration
      double precision x(100),w(100),coeff

      logical tp,tm
      integer mm(3),newn,n,itmp
      real*8   rtheta_lm,ytheta_lm
      external rtheta_lm,ytheta_lm

      pi = 4.0d0*datan(1.0d0)
      twopi   = 2.0d0*pi
      fourpi  = 4.0d0*pi
      piover2 = pi/2.0d0

      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) 
     > write(*,*)'Invalid parameter in gen_gaunt_coeff, negative l'
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2))
     1 write(*,*) 'Invalid parameter in gen_gaunt_coeff, m > l'

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gaunt_sub = 0.0d0
            return 
         else
           coeff = twopi
         endif
      else
         mm(1) = m1
         mm(2) = m2
         mm(3) = m3
         n = 3
         do while (n>1) 
            newn = 0
            do i=1,n-1
               if (mm(i).lt.mm(i+1)) then
                  itmp    = mm(i) 
                  mm(i)   = mm(i+1)
                  mm(i+1) = itmp
                  newn    = i
               end if
            end do
            n = newn
         end do
         tp = (abs(mm(1)).eq.(abs(mm(2))+abs(mm(3))) )
         tm = (abs(mm(1)).eq.abs((abs(mm(2))-abs(mm(3)))) )

         if ((mm(1).gt.0).and.(mm(2).gt.0).and.(mm(3).gt.0).and.tp) then
            coeff = piover2
         else if ((mm(1).gt.0).and.(mm(2).gt.0)
     >                        .and.(mm(3).eq.0).and.tp) then
            coeff = pi
         else if ((mm(1).gt.0).and.(mm(2).lt.0) 
     >                        .and.(mm(3).lt.0).and.tp) then
            coeff = -piover2
         else if ((mm(1).gt.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.tm) then
            coeff = piover2
         else if ((mm(1).eq.0).and.(mm(2).eq.0)
     >                        .and.(mm(3).eq.0)) then
            coeff = twopi
         else if ((mm(1).eq.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.(tm.or.tp)) then
            coeff = pi
         else
            nwpw_gaunt_sub = 0.0d0
            return 
         end if
      endif

      !**** Check the triangle rule ****
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         nwpw_gaunt_sub = 0.0d0
         return 
      endif

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         nwpw_gaunt_sub = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gaunt_sub = 1.0d0/dsqrt(fourpi)
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL nwpw_gauss_weights(x1,x2,x,w,order)
      nwpw_gaunt_sub = 0.0d0
      if (iscmplx) then
         do i = 1, order
            nwpw_gaunt_sub = nwpw_gaunt_sub
     >                 + w(i)*ytheta_lm(l1,m1,x(i)) 
     >                       *ytheta_lm(l2,m2,x(i))
     >                       *ytheta_lm(l3,m3,x(i))
         end do
      else
         do i = 1, order
            nwpw_gaunt_sub = nwpw_gaunt_sub
     >                 + w(i)*rtheta_lm(l1,m1,x(i)) 
     >                       *rtheta_lm(l2,m2,x(i))
     >                       *rtheta_lm(l3,m3,x(i))
         end do
      end if

      nwpw_gaunt_sub = nwpw_gaunt_sub*coeff

      return
      end 


c*     ***********************************
c*     *                                 *
c*     *          nwpw_gaunt2_init       *
c*     *                                 *
c*     ***********************************
c      subroutine nwpw_gaunt2_init(iscmplx,lmax0)
c      implicit none
c      logical iscmplx
c      integer lmax0
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt2_common / nwpw_gaunt_coeff,
c     >                             nwpw_gaunt_lmax,
c     >                             nwpw_gaunt_sizel2,
c     >                             nwpw_gaunt_iscmplx
c
c*     **** local variables ****
c      integer sizel,i,j,k,indx
c      integer l1,l2,l3
c      integer m1,m2,m3
c
c*     **** external functions ****
c      real*8   nwpw_gaunt2_sub
c      external nwpw_gaunt2_sub
c
c      nwpw_gaunt_iscmplx = iscmplx
c      nwpw_gaunt_lmax    = lmax0
c      sizel              = (nwpw_gaunt_lmax+1)
c      nwpw_gaunt_sizel2  = sizel**2
c
c      if (.not.BA_alloc_get(mt_dbl,
c     >      nwpw_gaunt_sizel2*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2,
c     >                    'nwpw_gaunt2_coeff',
c     >                     nwpw_gaunt_coeff(2),
c     >                     nwpw_gaunt_coeff(1)))
c     >   call errquit('nwpw_gaunt2_init:out of heap memory',0,MA_ERR)
c
c      do l1=0,nwpw_gaunt_lmax
c      do l2=0,nwpw_gaunt_lmax
c      do l3=0,nwpw_gaunt_lmax
c
c        do m1=-l1,l1
c        do m2=-l2,l2
c        do m3=-l3,l3
c           i = l1**2 + (l1+m1)
c           j = l2**2 + (l2+m2)
c           k = l3**2 + (l3+m3)
c           indx = i
c     >          + j*nwpw_gaunt_sizel2 
c     >          + k*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2
c
c           dbl_mb(nwpw_gaunt_coeff(1)+indx) 
c     >     = nwpw_gaunt2_sub(iscmplx,l1,m1,l2,m2,l3,m3)
c
c        end do
c        end do
c        end do
c
c      end do
c      end do
c      end do
c
c      return
c      end
c
c*     ***********************************
c*     *                                 *
c*     *          nwpw_gaunt2_end        *
c*     *                                 *
c*     ***********************************
c      subroutine nwpw_gaunt2_end()
c      implicit none
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt2_common / nwpw_gaunt_coeff,
c     >                              nwpw_gaunt_lmax,
c     >                              nwpw_gaunt_sizel2,
c     >                              nwpw_gaunt_iscmplx
c
c      if (.not.BA_free_heap(nwpw_gaunt_coeff(2)))
c     >   call errquit('nwpw_gaunt2_end:freeing heap memory',0,MA_ERR)
c
c      return
c      end
c
c*     ***********************************
c*     *                                 *
c*     *          nwpw_gaunt2            *
c*     *                                 *
c*     ***********************************
c      real*8 function nwpw_gaunt2(iscmplx,
c     >                           l1,m1,l2,m2,l3,m3)
c      implicit none
c      logical iscmplx
c      integer l1,l2,l3
c      integer m1,m2,m3
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt2_common / nwpw_gaunt_coeff,
c     >                             nwpw_gaunt_lmax,
c     >                             nwpw_gaunt_sizel2,
c     >                             nwpw_gaunt_iscmplx
c
c*     **** local variables ****
c      integer i,j,k,indx
c      real*8  gg
c
c*     **** external functions ****
c      real*8   nwpw_gaunt2_sub
c      external nwpw_gaunt2_sub
c
c      if ((l1.gt.nwpw_gaunt_lmax).or.
c     >    (l2.gt.nwpw_gaunt_lmax).or.
c     >    (l3.gt.nwpw_gaunt_lmax).or.
c     >    (iscmplx.neqv.nwpw_gaunt_iscmplx)) then
c        gg = nwpw_gaunt2_sub(iscmplx,l1,m1,l2,m2,l3,m3)
c      else
c         i = l1**2 + (l1+m1)
c         j = l2**2 + (l2+m2)
c         k = l3**2 + (l3+m3)
c         indx = i
c     >        + j*nwpw_gaunt_sizel2
c     >        + k*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2
c         gg = dbl_mb(nwpw_gaunt_coeff(1)+indx)
c
c      end if
c      nwpw_gaunt2 = gg
c
c      return
c      end

*     ***********************************
*     *                                 *
*     *          nwpw_gaunt2            *
*     *                                 *
*     ***********************************
*
* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt2(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi 
*     |     |
*   = |     | dconjg(dY_lm(theta,phi)/dtheta) * dY_l1m1(theta,phi)/dtheta * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | dY_lm(theta,phi)/dtheta * dconjg(dY_l1m1(theta,phi)/dtheta) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*
* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt2(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | dT_lm(theta,phi))/dtheta * dT_l1m1(theta,phi)/dtheta * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
      real*8 function nwpw_gaunt2(iscmplx,
     >                            l1,m1,l2,m2,l3,m3)
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     
*     *** local variables ****
      double precision pi
      double precision piover2,twopi,fourpi
c      parameter(pi = 3.14159265358979323846264338327950288419d0)
c      parameter(piover2= pi*0.5d0)
c      parameter(twopi  = pi*2.0d0)
c      parameter(fourpi = pi*4.0d0)

c     work arrays for integration
      double precision x(100),w(100),coeff

      logical tp,tm
      integer mm(3),newn,n,itmp
      real*8   rtheta_lm,rtheta2_lm,ytheta_lm,ytheta2_lm
      external rtheta_lm,rtheta2_lm,ytheta_lm,ytheta2_lm

      pi = 4.0d0*datan(1.0d0)
      twopi   = 2.0d0*pi
      fourpi  = 4.0d0*pi
      piover2 = pi/2.0d0

      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) 
     > write(*,*)'Invalid parameter in nwpw_gaunt2, negative l'

      If (l1.lt.(abs(m1)+1).or.l2.lt.(abs(m2)+1).or.l3.lt.abs(m3)) then
        nwpw_gaunt2 = 0.0d0
        return
      end if

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gaunt2 = 0.0d0
            return 
         else
           coeff = twopi
         endif
      else
         mm(1) = m1
         mm(2) = m2
         mm(3) = m3
         n = 3
         do while (n>1) 
            newn = 0
            do i=1,n-1
               if (mm(i).lt.mm(i+1)) then
                  itmp    = mm(i) 
                  mm(i)   = mm(i+1)
                  mm(i+1) = itmp
                  newn    = i
               end if
            end do
            n = newn
         end do
         tp = (abs(mm(1)).eq.(abs(mm(2))+abs(mm(3))) )
         tm = (abs(mm(1)).eq.abs((abs(mm(2))-abs(mm(3)))) )

         if ((mm(1).gt.0).and.(mm(2).gt.0).and.(mm(3).gt.0).and.tp) then
            coeff = piover2
         else if ((mm(1).gt.0).and.(mm(2).gt.0)
     >                        .and.(mm(3).eq.0).and.tp) then
            coeff = pi
         else if ((mm(1).gt.0).and.(mm(2).lt.0) 
     >                        .and.(mm(3).lt.0).and.tp) then
            coeff = -piover2
         else if ((mm(1).gt.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.tm) then
            coeff = piover2
         else if ((mm(1).eq.0).and.(mm(2).eq.0)
     >                        .and.(mm(3).eq.0)) then
            coeff = twopi
         else if ((mm(1).eq.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.(tm.or.tp)) then
            coeff = pi
         else
            nwpw_gaunt2 = 0.0d0
            return 
         end if
      endif

      !**** Check the triangle rule ****
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL nwpw_gauss_weights(x1,x2,x,w,order)
      nwpw_gaunt2 = 0.0d0
      if (iscmplx) then
         do i = 1, order
            nwpw_gaunt2 = nwpw_gaunt2
     >                 + w(i)*ytheta2_lm(l1,m1,x(i)) 
     >                       *ytheta2_lm(l2,m2,x(i))
     >                       *ytheta_lm(l3,m3,x(i))
         end do
      else
         do i = 1, order
            nwpw_gaunt2 = nwpw_gaunt2
     >                 + w(i)*rtheta2_lm(l1,m1,x(i)) 
     >                       *rtheta2_lm(l2,m2+1,x(i))
     >                       *rtheta_lm(l3,m3,x(i))
         end do
      end if

      nwpw_gaunt2 = nwpw_gaunt2*coeff

      return
      end 


*     ***********************************
*     *                                 *
*     *          nwpw_gaunt3_init       *
*     *                                 *
*     ***********************************
c      subroutine nwpw_gaunt3_init(iscmplx,lmax0)
c      implicit none
c      logical iscmplx
c      integer lmax0
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt3_common / nwpw_gaunt_coeff,
c     >                             nwpw_gaunt_lmax,
c     >                             nwpw_gaunt_sizel2,
c     >                             nwpw_gaunt_iscmplx
c
c*     **** local variables ****
c      integer sizel,i,j,k,indx
c      integer l1,l2,l3
c      integer m1,m2,m3
c
c*     **** external functions ****
c      real*8   nwpw_gaunt3_sub
c      external nwpw_gaunt3_sub
c
c      nwpw_gaunt_iscmplx = iscmplx
c      nwpw_gaunt_lmax    = lmax0
c      sizel              = (nwpw_gaunt_lmax+1)
c      nwpw_gaunt_sizel2  = sizel**2
c
c      if (.not.BA_alloc_get(mt_dbl,
c     >      nwpw_gaunt_sizel2*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2,
c     >                    'nwpw_gaunt3_coeff',
c     >                     nwpw_gaunt_coeff(2),
c     >                     nwpw_gaunt_coeff(1)))
c     >   call errquit('nwpw_gaunt3_init:out of heap memory',0,MA_ERR)
c
c      do l1=0,nwpw_gaunt_lmax
c      do l2=0,nwpw_gaunt_lmax
c      do l3=0,nwpw_gaunt_lmax
c
c        do m1=-l1,l1
c        do m2=-l2,l2
c        do m3=-l3,l3
c           i = l1**2 + (l1+m1)
c           j = l2**2 + (l2+m2)
c           k = l3**2 + (l3+m3)
c           indx = i
c     >          + j*nwpw_gaunt_sizel2 
c     >          + k*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2
c
c           dbl_mb(nwpw_gaunt_coeff(1)+indx) 
c     >     = nwpw_gaunt3_sub(iscmplx,l1,m1,l2,m2,l3,m3)
c
c        end do
c        end do
c        end do
c
c      end do
c      end do
c      end do
c
c      return
c      end
c
c*     ***********************************
c*     *                                 *
c*     *          nwpw_gaunt3_end        *
c*     *                                 *
c*     ***********************************
c      subroutine nwpw_gaunt3_end()
c      implicit none
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt3_common / nwpw_gaunt_coeff,
c     >                              nwpw_gaunt_lmax,
c     >                              nwpw_gaunt_sizel2,
c     >                              nwpw_gaunt_iscmplx
c
c
c      if (.not.BA_free_heap(nwpw_gaunt_coeff(2)))
c     >   call errquit('nwpw_gaunt3_end:freeing heap memory',0,MA_ERR)
c
c      return
c      end
c
c*     ***********************************
c*     *                                 *
c*     *          nwpw_gaunt3             *
c*     *                                 *
c*     ***********************************
c      real*8 function nwpw_gaunt3(iscmplx,
c     >                           l1,m1,l2,m2,l3,m3)
c      implicit none
c      logical iscmplx
c      integer l1,l2,l3
c      integer m1,m2,m3
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c*     **** nwpw_gaunt common block ****
c      integer nwpw_gaunt_lmax,
c     >        nwpw_gaunt_sizel2,
c     >        nwpw_gaunt_coeff(2)
c      logical nwpw_gaunt_iscmplx
c      common / nwpw_gaunt3_common / nwpw_gaunt_coeff,
c     >                             nwpw_gaunt_lmax,
c     >                             nwpw_gaunt_sizel2,
c     >                             nwpw_gaunt_iscmplx
c
c*     **** local variables ****
c      integer i,j,k,indx
c      real*8  gg
c
c*     **** external functions ****
c      real*8   nwpw_gaunt3_sub
c      external nwpw_gaunt3_sub
c
c      if ((l1.gt.nwpw_gaunt_lmax).or.
c     >    (l2.gt.nwpw_gaunt_lmax).or.
c     >    (l3.gt.nwpw_gaunt_lmax).or.
c     >    (iscmplx.neqv.nwpw_gaunt_iscmplx)) then
c        gg = nwpw_gaunt3_sub(iscmplx,l1,m1,l2,m2,l3,m3)
c      else
c         i = l1**2 + (l1+m1)
c         j = l2**2 + (l2+m2)
c         k = l3**2 + (l3+m3)
c         indx = i
c     >        + j*nwpw_gaunt_sizel2
c     >        + k*nwpw_gaunt_sizel2*nwpw_gaunt_sizel2
c         gg = dbl_mb(nwpw_gaunt_coeff(1)+indx)
c
c      end if
c      nwpw_gaunt3 = gg
c
c      return
c      end


*     ***********************************
*     *                                 *
*     *          nwpw_gaunt3            *
*     *                                 *
*     ***********************************
*
* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt3(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi 
*     |     |
*   = |     | dconjg(dY_lm(theta,phi)/dphi) * dY_l1m1(theta,phi)/dphi * dconjg(Y_l2m2(theta,phi)) * 1/(sin(theta))**2 *  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | dY_lm(theta,phi)/dphi * dconjg(dY_l1m1(theta,phi)/dphi) * Y_l2m2(theta,phi) 1/(sin(theta))**2  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*
* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt3(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | (1/sin(theta))**2 * dT_lm(theta,phi))/dphi * dT_l1m1(theta,phi)/dphi * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
      real*8 function nwpw_gaunt3(iscmplx,
     >                           l1,m1,l2,m2,l3,m3)
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     
*     *** local variables ****
      double precision pi
      double precision piover2,twopi,fourpi
c      parameter(pi = 3.14159265358979323846264338327950288419d0)
c      parameter(piover2= pi*0.5d0)
c      parameter(twopi  = pi*2.0d0)
c      parameter(fourpi = pi*4.0d0)

c     work arrays for integration
      double precision x(100),w(100),coeff

      logical tp,tm
      integer mm(3),newn,n,itmp
      real*8   rtheta_lm,rtheta_lm_div,ytheta_lm,theta_lm_div
      external rtheta_lm,rtheta_lm_div,ytheta_lm,theta_lm_div

      pi = 4.0d0*datan(1.0d0)
      twopi   = 2.0d0*pi
      fourpi  = 4.0d0*pi
      piover2 = pi/2.0d0

      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) 
     > write(*,*)'Invalid parameter in nwpw_gaunt3, negative l'
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2))
     1 write(*,*) 'Invalid parameter in nwpw_gaunt3, m > l'

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gaunt3 = 0.0d0
            return 
         else
           coeff = twopi
         endif
      else
         m1 = -m1
         m2 = -m2
         mm(1) = m1
         mm(2) = m2
         mm(3) = m3
         n = 3
         do while (n>1) 
            newn = 0
            do i=1,n-1
               if (mm(i).lt.mm(i+1)) then
                  itmp    = mm(i) 
                  mm(i)   = mm(i+1)
                  mm(i+1) = itmp
                  newn    = i
               end if
            end do
            n = newn
         end do
         tp = (abs(mm(1)).eq.(abs(mm(2))+abs(mm(3))) )
         tm = (abs(mm(1)).eq.abs((abs(mm(2))-abs(mm(3)))) )

         if ((mm(1).gt.0).and.(mm(2).gt.0).and.(mm(3).gt.0).and.tp) then
            coeff = piover2
         else if ((mm(1).gt.0).and.(mm(2).gt.0)
     >                        .and.(mm(3).eq.0).and.tp) then
            coeff = pi
         else if ((mm(1).gt.0).and.(mm(2).lt.0) 
     >                        .and.(mm(3).lt.0).and.tp) then
            coeff = -piover2
         else if ((mm(1).gt.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.tm) then
            coeff = piover2
         else if ((mm(1).eq.0).and.(mm(2).eq.0)
     >                        .and.(mm(3).eq.0)) then
            coeff = twopi
         else if ((mm(1).eq.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.(tm.or.tp)) then
            coeff = pi
         else
            nwpw_gaunt3 = 0.0d0
            return 
         end if
      endif

c      !**** Check the triangle rule ****
c      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
c         nwpw_gaunt3 = 0.0d0
c         return 
c      endif

c      !**** Check if the integrand is odd function==>integral is zero ****
c      if (mod(l1 + l2 + l3,2) .eq. 1) then
c         nwpw_gaunt3 = 0.0d0
c         return 
c      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gaunt3 = 0.0d0
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL nwpw_gauss_weights(x1,x2,x,w,order)
      nwpw_gaunt3 = 0.0d0
      if (iscmplx) then
         do i = 1, order
            nwpw_gaunt3 = nwpw_gaunt3
     >                 - w(i)*theta_lm_div(l1,m1,x(i)) * (m1)
     >                       *theta_lm_div(l2,m2,x(i)) * (m2)
     >                       *ytheta_lm(l3,m3,x(i))
         end do
      else
         do i = 1, order
            nwpw_gaunt3 = nwpw_gaunt3
     >                 + w(i)*rtheta_lm_div(l1,-m1,x(i)) * (m1)
     >                       *rtheta_lm_div(l2,-m2,x(i)) * (m2)
     >                       *rtheta_lm(l3,m3,x(i))
         end do
      end if

      nwpw_gaunt3 = nwpw_gaunt3*coeff

      return
      end 




c!  *************************************************
c!  
c!     Name    : rtheta2_lm
c!  
c!  
c!     Purpose : calculates rtheta2_lm for a scalar cos_theta
c!               such that
c!                                                           {cos(|m|*phi)   m>0
c!         dT_lm(cos_theta,phi)/dtheta=rtheta2_lm(cos_theta)*{1              m==0
c!                                                           {sin(|m|*phi)   m<0
c!  
c!  
c!  *************************************************
      double precision function rtheta2_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f

      !*** external functions *** 
      double precision rlegendre_lm
      external         rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if ((mod_m+1) .gt. l) then
         rtheta2_lm = 0.0d0
      end if

!     *** find coefficient ***
      if (mod_m.eq.0) then 
         coeff= 0.5d0
      else
         coeff= 1.0d0 
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)
      f     = coeff*rlegendre_lm(l,mod_m+1,cos_theta)

      rtheta2_lm = -f
      return
      end 


c!  *************************************************
c!  
c!     Name    : ytheta2_lm
c!  
c!  
c!     Purpose : calculates ytheta_lm for a scalar cos_theta
c!               such that
c!                                                 
c!         dY_lm(cos_theta,phi)/dtheta=ytheta_lm(cos_theta)*exp(i*m*phi)
c!                                                
c!  
c!  *************************************************
      double precision function ytheta2_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f

      !*** external functions *** 
      double precision legendre_lm
      external         legendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      mod_m = abs(m) 
      if (mod_m+1 .gt. l) then
         ytheta2_lm = 0.0d0
         return
      end if

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 1.0d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif 
      coeff = coeff*(2*l+1)/fourpi
      coeff = dsqrt(coeff)
      f     = coeff*legendre_lm(l,mod_m+1,cos_theta)
      if ((m.lt.0).and.(mod(mod_m+1,2).eq.1)) f=-f

      ytheta2_lm = -f
      return
      end






*     ***************************
*     *				*
*     *   nwpw_factorial_init   *
*     *				*
*     ***************************

      subroutine nwpw_factorial_init()
      implicit none


*     **** Factorial common block ****
      integer FMAX
      parameter (FMAX=30)
      real*8 factorial(0:FMAX)
      common / FACTORIAL_Block / factorial

*     **** local variables ****
      integer n
      real*8  x


      factorial(0) = 1.0d0
      do n=1,FMAX
         x = dble(n)
         factorial(n) = x*factorial(n-1)
      end do

      return
      end

*     ***************************
*     *				*
*     *      nwpw_clebsch	*
*     *				*
*     ***************************
*
* Computes the Clebsch-Gordon coefficients using Racah's formula.
*

      real*8 function nwpw_clebsch(l1,l2,l,m1,m2,m)
      implicit none
      integer l1,l2,l
      integer m1,m2,m


*     **** Factorial common block ****
      integer FMAX
      parameter (FMAX=30)
      real*8 factorial(0:FMAX)
      common / FACTORIAL_Block / factorial
      logical genfac
      save genfac
      data genfac /.true./

*     **** local variables ****
      integer z,zmin,zmax
      real*8 tmp1,tmp2,sign

      if (genfac) then
         call nwpw_factorial_init()
         genfac = .false.
      end if

*     **************************************
*     **** check the triangle condition ****
*     **************************************

*     *** a zero Clebsch-Gordon coefficient ***
      if ( ((m1+m2).ne.m)      .OR.
     >     (abs(m1).gt.abs(l1)).OR.
     >     (abs(m2).gt.abs(l2)).OR.
     >     (abs(m) .gt.abs(l) ).OR.
     >     (l.gt.(l1+l2))      .OR.
     >     (l.lt.abs(l1-l2))   .OR.
     >     (mod((l1+l2+l),2).ne. 0) ) then
         tmp2 = 0.0d0

*     *** a non-zero Clebsch-Gordon coeffient ****
      else
         zmin = 0
         if ((l-l2+m1)     .lt.0) zmin = -l + l2 - m1
         if ((l-l1-m2+zmin).lt.0) zmin = -l + l1 + m2
         zmax = l1+l2-l
         if ((l2+m2-zmax).lt.0) zmax = l2 + m2
         if ((l1-m1-zmax).lt.0) zmax = l1 - m1

         tmp1 = 0.0d0
         do z=zmin,zmax
            sign = 1.0d0
            if (mod(z,2).ne.0) sign = -1.0d0
            tmp1 = tmp1 
     >           + sign
     >             / ( factorial(  (z)        )
     >               * factorial( (l1+l2-l-z) )
     >               * factorial( (l1-m1-z)   )
     >               * factorial( (l2+m2-z)   )
     >               * factorial( (l-l2+m1+z) )
     >               * factorial( (l-l1-m2+z) ))
         end do
         tmp2 = dsqrt( dble(2*l+1) * factorial( (l1+l2-l)   )
     >                           * factorial( (l1-l2+l)   )
     >                           * factorial( (-l1+l2+l)  )
     >                           * factorial( (l1+m1)     )
     >                           * factorial( (l1-m1)     )
     >                           * factorial( (l2+m2)     )
     >                           * factorial( (l2-m2)     )
     >                           * factorial( (l+m)       )
     >                           * factorial( (l-m)       )
     >                           / factorial( (l1+l2+l+1) )
     >               )*tmp1
      end if

      nwpw_clebsch = tmp2
      return
      end



*     ***************************
*     *				*
*     *      nwpw_Gaunt_Clebsh	*
*     *				*
*     ***************************
*
* Computes the Gaunt coefficient
*
*  Gaunt(l,m,l1,m1,l2,m2) =

*      /2pi  / pi
*     |     |
*   = |     | dconjg(Y_lm(theta,phi)) * Y_l1m1(theta,phi) * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | Y_lm(theta,phi) * dconjg(Y_l1m1(theta,phi)) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*   = sqrt( (2*l+1)*(2*l2+1)/(4*pi*(2*l1+1)) ) * Clebsch(l l2 l1; m m2 m1) * Clebsch(l l2 l1; 0 0 0)
*
      real*8 function nwpw_Gaunt_Clebsh(l,m,l1,m1,l2,m2)
      implicit none
      integer l,m
      integer l1,m1
      integer l2,m2

*     **** local variables ****
      real*8 coeff,fourpi
      real*8 tmp,xsign
   
*     **** external functions ***
      real*8   nwpw_clebsch
      external nwpw_clebsch


*     *** a zero Gaunt coefficient ***
      if ( ((m+m2).ne.m1)      .OR.
     >     (abs(m1).gt.abs(l1)).OR.
     >     (abs(m2).gt.abs(l2)).OR.
     >     (abs(m) .gt.abs(l) ).OR.
     >     (l1.gt.(l+l2))      .OR.
     >     (l1.lt.abs(l-l2))   .OR.
     >     (mod((l1+l2+l),2).ne. 0 )) then
         tmp = 0.0d0

*     *** a non-zero Gaunt Coefficient ***
      else
c         if (mod((m*m+m2*m2-3*m-3*m2)/2,2).eq.1) then
c            xsign = -1.0d0
c         else
c            xsign = 1.0d0
c         end if
          xsign = 1.0d0

         fourpi = 16.0d0*datan(1.0d0)
         coeff  = dsqrt( dble((2*l+1)*(2*l2+1))/(fourpi*dble(2*l1+1)))
         tmp    = xsign*coeff*nwpw_clebsch(l,l2,l1, m,m2,m1)
     >                       *nwpw_clebsch(l,l2,l1, 0,0,0)
      end if

      nwpw_Gaunt_Clebsh = tmp
      return
      end

*     **************************************
*     *                                    *
*     *      rSphericalHarmonic3           *
*     *                                    *
*     **************************************

*  Calculates the r*Spherical Harmonic for x,y,z such that
*
*                                                                                 {exp(i*|m|*phi)  m>0
*     rY_lm(cos_theta,phi)=(r**l)*sqrt((l-|m|)!/(l+|m|)!)*rlegendre_lm(cos_theta)*{1               m==0
*                                                                                 {exp(-i*|m|*phi) m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*

      complex*16 function rSphericalHarmonic3(l,m,x,y,z)
      implicit none
      integer l,m
      real*8 x,y,z

*     **** local variables ****
      integer i,mod_m
      real*8 r,cos_theta,phi,coeff
      complex*16 ctmp

*     **** external functions ****
      real*8   rlegendre_lm
      external rlegendre_lm

      mod_m = abs(m)
      if (mod_m .gt. l) 
     >  write(*,*) 
     >  'parameter out of order in function rSphericalHarmonic3'

      r = dsqrt(x*x + y*y + z*z)
      if (r.gt.1.0d-9) then
            cos_theta = z/r
            phi = datan2(y,x)
*        *** find coefficient ***
         if (mod_m.eq.0) then 
            coeff= 1.0d0
            ctmp = dcmplx(1.0d0,0.0d0)
         else if (mod_m.gt.0) then
            coeff= 1.0d0 
            do i=1,2*mod_m
               coeff = coeff/dble(l-mod_m+i)
            end do
            ctmp = dcmplx(dcos(m*phi),dsin(m*phi))
         end if
         coeff = dsqrt(coeff)
         rSphericalHarmonic3 = ctmp*coeff
     >                        *rlegendre_lm(l,mod_m,cos_theta)
     >                        *(r**l)
      else
         if (l.eq.0) then
            ctmp = dcmplx(1.0d0,0.0d0)
         else
            ctmp = dcmplx(0.0d0,0.0d0)
         end if
         rSphericalHarmonic3=ctmp*(r**l)
      end if
      return
      end


*     **************************************
*     *                                    *
*     *      Tesseral3_vector_lm           *
*     *                                    *
*     **************************************

*  Calculates the tesseral harmonic for x,y,z such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine Tesseral3_vector_lm(l,m,nfft3d,x,y,z,Tlm)
      implicit none
      integer l,m
      integer nfft3d
      real*8 x(*),y(*),z(*),Tlm(*)

*     **** local variables ****
      integer i,k,mod_m
      real*8 phi,cos_theta,tmp2,r,twopi,coeff

*     **** external functions
      real*8   rlegendre_lm,rtheta_lm
      external rlegendre_lm,rtheta_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l) 
     >  write(*,*) 
     >  'parameter out of order in function Tesseral3_vector_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then 
         coeff= 0.5d0
      else
         coeff= 1.0d0 
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)


      do k=1,nfft3d
         r = dsqrt(x(k)*x(k) + y(k)*y(k) + z(k)*z(k))
         if (r.gt.1.0d-9) then
            cos_theta = z(k)/r
            phi = datan2(y(k),x(k))

            if (m.lt.0) then
               tmp2 = dsin(mod_m*phi)
            else if (m.gt.0) then
               tmp2 = dcos(mod_m*phi)
            else
               tmp2 = 1.0d0
            end if

            Tlm(k) = tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)
         else
            if (l.eq.0) then
               Tlm(k) = 1.0d0/dsqrt(2*twopi)
            else
               Tlm(k) = 0.0d0
            end if
         end if
      end do
      return
      end

*     **************************************
*     *                                    *
*     *       Tesseral3_rgrid_lm           *
*     *                                    *
*     **************************************

*  Calculates the tesseral harmonic for x,y,z such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine Tesseral3_rgrid_lm(l,m,nfft3d,rgrid,Tlm)
      implicit none
      integer l,m
      integer nfft3d
      real*8 rgrid(3,*),Tlm(*)

*     **** local variables ****
      integer i,k,mod_m
      real*8 phi,cos_theta,tmp2,r,twopi,coeff

*     **** external functions
      real*8   rlegendre_lm,rtheta_lm
      external rlegendre_lm,rtheta_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l)
     >  write(*,*)
     >  'parameter out of order in function Tesseral3_vector_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)


      do k=1,nfft3d
         r = dsqrt(rgrid(1,k)**2 + rgrid(2,k)**2 + rgrid(3,k)**2)
         if (r.gt.1.0d-9) then
            cos_theta = rgrid(3,k)/r
            phi = datan2(rgrid(2,k),rgrid(1,k))

            if (m.lt.0) then
               tmp2 = dsin(mod_m*phi)
            else if (m.gt.0) then
               tmp2 = dcos(mod_m*phi)
            else
               tmp2 = 1.0d0
            end if

            Tlm(k) = tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)
         else
            if (l.eq.0) then
               Tlm(k) = 1.0d0/dsqrt(2*twopi)
            else
               Tlm(k) = 0.0d0
            end if
         end if
      end do
      return
      end


*     **************************************
*     *                                    *
*     *       Tesseral3_rgrid_lm_rl        *
*     *                                    *
*     **************************************

*  Calculates the r**l times the tesseral harmonic for x,y,z such that
*
*                                                         {cos(|m|*phi)   m>0
*      r**l*T_lm(cos_theta,phi)=r**l*rtheta_lm(cos_theta)*{1              m==0
*                                                         {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine Tesseral3_rgrid_lm_rl(l,m,nfft3d,rgrid,rTlm)
      implicit none
      integer l,m
      integer nfft3d
      real*8 rgrid(3,*),rTlm(*)

*     **** local variables ****
      integer i,k,mod_m
      real*8 phi,cos_theta,tmp2,r,twopi,coeff

*     **** external functions
      real*8   rlegendre_lm,rtheta_lm
      external rlegendre_lm,rtheta_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l)
     >  write(*,*)
     >  'parameter out of order in function Tesseral3_vector_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else 
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)

      do k=1,nfft3d
         r = dsqrt(rgrid(1,k)**2 + rgrid(2,k)**2 + rgrid(3,k)**2)
         if (r.gt.1.0d-9) then
            cos_theta = rgrid(3,k)/r
            phi = datan2(rgrid(2,k),rgrid(1,k))

            if (m.lt.0) then
               tmp2 = dsin(mod_m*phi)
            else if (m.gt.0) then
               tmp2 = dcos(mod_m*phi)
            else
               tmp2 = 1.0d0
            end if

            rTlm(k) = (r**l)*tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)
         else
            rTlm(k) = 0.0d0
         end if
      end do
      return
      end


*     **************************************
*     *                                    *
*     *     Tesseral3_rgrid_lm_over_rl     *
*     *                                    *
*     **************************************

*  Calculates the 1/r**(l+1) times the tesseral harmonic for x,y,z such that
*
*                                                                     {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)/r**(l+1)=(1/r**(l+1))*rtheta_lm(cos_theta)*{1              m==0
*                                                                     {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine Tesseral3_rgrid_lm_over_rl(l,m,nfft3d,rgrid,rTlm)
      implicit none
      integer l,m
      integer nfft3d
      real*8 rgrid(3,*),rTlm(*)

*     **** local variables ****
      integer i,k,mod_m
      real*8 phi,cos_theta,tmp2,r,twopi,coeff

*     **** external functions
      real*8   rlegendre_lm,rtheta_lm
      external rlegendre_lm,rtheta_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l)
     >  write(*,*)
     >  'parameter out of order in function Tesseral3_vector_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)


      do k=1,nfft3d
         r = dsqrt(rgrid(1,k)**2 + rgrid(2,k)**2 + rgrid(3,k)**2)
         if (r.gt.1.0d-9) then
            cos_theta = rgrid(3,k)/r
            phi = datan2(rgrid(2,k),rgrid(1,k))

            if (m.lt.0) then
               tmp2 = dsin(mod_m*phi)
            else if (m.gt.0) then
               tmp2 = dcos(mod_m*phi)
            else
               tmp2 = 1.0d0
            end if

            rTlm(k) = tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)
     >                /(r**(l+1))
         else
            rTlm(k) = 0.0d0
         end if
      end do
      return
      end



*     **************************************
*     *                                    *
*     *           Tesseral3_lm             *
*     *                                    *
*     **************************************

*  Calculates the tesseral harmonic for x,y,z such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      double precision function Tesseral3_lm(l,m,x,y,z)
      implicit none
      integer l,m
      real*8 x,y,z

*     **** local variables ****
      integer mod_m
      real*8 phi,cos_theta,tmp2,r

*     **** external functions
      real*8   rtheta_lm
      external rtheta_lm

      r = dsqrt(x*x + y*y + z*z)
      cos_theta = z/r
      phi = datan2(y,x)
      
      mod_m = abs(m)
      if (m.lt.0) then
         tmp2 = dsin(mod_m*phi)
      else if (m.gt.0) then
         tmp2 = dcos(mod_m*phi)
      else
         tmp2 = 1.0d0
      end if

      Tesseral3_lm = rtheta_lm(l,m,cos_theta)*tmp2
      return
      end



*     **************************************
*     *                                    *
*     *       Tesseral3_lm_over_rl         *
*     *                                    *
*     **************************************

*  Calculates the tesseral harmonic divided by r**(l+1) for x,y,z such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*  
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      double precision function Tesseral3_lm_over_rl(l,m,x,y,z)
      implicit none
      integer l,m
      real*8 x,y,z

*     **** local variables ****
      integer mod_m
      real*8 phi,cos_theta,tmp2,r

*     **** external functions
      real*8   rtheta_lm
      external rtheta_lm

      r = dsqrt(x*x + y*y + z*z)
      cos_theta = z/r
      phi = datan2(y,x)

      mod_m = abs(m)
      if (m.lt.0) then
         tmp2 = dsin(mod_m*phi)
      else if (m.gt.0) then
         tmp2 = dcos(mod_m*phi)
      else
         tmp2 = 1.0d0
      end if

      Tesseral3_lm_over_rl = rtheta_lm(l,m,cos_theta)*tmp2/(r**(l+1))
      return
      end


*     **************************************
*     *                                    *
*     *           dTesseral3_lm            *
*     *                                    *
*     **************************************
*
*  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine dTesseral3_lm(l,m,x,y,z,dTx,dTy,dTz)
      implicit none
      integer l,m
      real*8 x,y,z
      real*8 dTx,dTy,dTz

*     **** local variables ****
      integer mod_m
      real*8 phi,cos_theta,tmp2,r,f1,f2

*     **** external functions
      real*8   rtheta_lm_div,drtheta_lm
      external rtheta_lm_div,drtheta_lm

      r = dsqrt(x*x + y*y + z*z)
      cos_theta = z/r
      phi = datan2(y,x)

      mod_m = abs(m)

      if (m.lt.0) then
         tmp2 = dsin(mod_m*phi)
         f2   = mod_m *dcos(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
      else if (m.gt.0) then
         tmp2 = dcos(mod_m*phi)
         f2   = -mod_m*dsin(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
      else
         tmp2 = 1.0d0
         f2   = 0.0d0
      end if
      f1 = drtheta_lm(l,m,cos_theta)*tmp2

      dTx = f1*cos_theta*cos(phi) - f2*sin(phi)
      dTy = f1*cos_theta*sin(phi) + f2*cos(phi)
      dTz = -f1*dsqrt(1.0d0-cos_theta**2)

      return
      end

*     **************************************
*     *                                    *
*     *           dTesseral3_lm_over_rl    *
*     *                                    *
*     **************************************
*
*  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such that
*
*                                                {cos(|m|*phi)   m>0
*      dT_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                                {sin(|m|*phi)   m<0
*
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine dTesseral3_lm_over_rl(l,m,x,y,z,dTx,dTy,dTz)
      implicit none
      integer l,m
      real*8 x,y,z
      real*8 dTx,dTy,dTz

*     **** local variables ****
      integer mod_m,i
      real*8 phi,cos_theta,tmp2,r,f1,f2,rl1,drl1,drTlm,coeff,twopi

*     **** external functions
      real*8   rtheta_lm_div,drtheta_lm,rlegendre_lm
      external rtheta_lm_div,drtheta_lm,rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)
      mod_m = abs(m)
      if (mod_m .gt. l)
     >  write(*,*)
     >  'parameter out of order in function dTesseral3_lm_over_rl'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else 
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      end if
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)

      r = dsqrt(x*x + y*y + z*z)
      if (r.gt.1.0d-9) then
         rl1 = 1.0d0/(r**(l+1))
         drl1 = -dble(l+1)/(r**(l+2))

         cos_theta = z/r
         phi = datan2(y,x)


         if (m.lt.0) then
            tmp2 = dsin(mod_m*phi)
            f2   = mod_m *dcos(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
         else if (m.gt.0) then
            tmp2 = dcos(mod_m*phi)
            f2   = -mod_m*dsin(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
         else
            tmp2 = 1.0d0
            f2   = 0.0d0
         end if
         drTlm = tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)*drl1
         f1 = drtheta_lm(l,m,cos_theta)*tmp2

         dTx = drTlm*x/r + (f1*cos_theta*cos(phi)-f2*sin(phi))*rl1
         dTy = drTlm*y/r + (f1*cos_theta*sin(phi)+f2*cos(phi))*rl1
         dTz = drTlm*z/r - f1*dsqrt(1.0d0-cos_theta**2)*rl1
      else
         dTx = 0.0d0
         dTy = 0.0d0
         dTz = 0.0d0
      end if
      return
      end




*     **************************************
*     *                                    *
*     *      dTesseral3_rgrid_lm_rl        *
*     *                                    *
*     **************************************
*
      subroutine dTesseral3_rgrid_lm_rl(l,m,nfft3d,rgrid,drTlm)
      implicit none
      integer l,m
      integer nfft3d
      real*8 rgrid(3,*)
      real*8 drTlm(3,*)

*     **** local variables ****
      integer mod_m,i,k
      real*8 phi,cos_theta,tmp2,r,f1,f2,Tlm
      real*8 rl,drl,coeff,twopi

*     **** external functions
      real*8   rtheta_lm_div,drtheta_lm,rlegendre_lm
      external rtheta_lm_div,drtheta_lm,rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l)
     >  write(*,*)
     >  'parameter out of order in function dTesseral3_rgrid_lm_rl'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)

      do k=1,nfft3d
         r = dsqrt(rgrid(1,k)**2 + rgrid(2,k)**2 + rgrid(3,k)**2)

         if (r.gt.1.0d-9) then
            rl = (r**l)
            if (l.eq.0) then
               drl = 0.0d0
            else if (l.eq.1) then
               drl = 1.0d0
            else
               drl = dble(l)*(r**(l-1))
            end if

            cos_theta = rgrid(3,k)/r
            phi = datan2(rgrid(2,k),rgrid(1,k))

            if (m.lt.0) then
               tmp2= dsin(mod_m*phi)
               f2  = mod_m*dcos(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
            else if (m.gt.0) then
               tmp2= dcos(mod_m*phi)
               f2  = -mod_m*dsin(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
            else
               tmp2= 1.0d0
               f2  = 0.0d0
            end if
            Tlm = drl*tmp2*coeff*rlegendre_lm(l,mod_m,cos_theta)
            f1  = drtheta_lm(l,m,cos_theta)*tmp2

            drTlm(1,k) = Tlm*rgrid(1,k)/r
     >                 + rl*(f1*cos_theta*cos(phi) - f2*sin(phi))
            drTlm(2,k) = Tlm*rgrid(2,k)/r
     >                 + rl*(f1*cos_theta*sin(phi) + f2*cos(phi))
            drTlm(3,k) = Tlm*rgrid(3,k)/r
     >                 - rl*f1*dsqrt(1.0d0-cos_theta**2)
         else
            drTlm(1,k) = 0.0d0
            drTlm(2,k) = 0.0d0
            drTlm(3,k) = 0.0d0
         end if

      end do

      return
      end




*     **************************************
*     *                                    *
*     *           Tesseral_lm              *
*     *                                    *
*     **************************************

*  Calculates the tesseral harmonic for x,y,z such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
* 
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      double precision function Tesseral_lm(l,m,cos_theta,phi)
      implicit none
      integer l,m
      real*8 cos_theta,phi

*     **** local variables ****
      integer mod_m
      real*8 tmp2,r

*     **** external functions
      real*8   rtheta_lm
      external rtheta_lm

      mod_m = abs(m)
      if (m.lt.0) then
         tmp2 = dsin(mod_m*phi)
      else if (m.gt.0) then
         tmp2 = dcos(mod_m*phi)
      else
         tmp2 = 1.0d0
      end if

      Tesseral_lm = rtheta_lm(l,m,cos_theta)*tmp2
      return
      end


*     **************************************
*     *                                    *
*     *           YSpherical_lm            *
*     *                                    *
*     **************************************

*  Calculates the spherical harmonic for x,y,z such that
*
*      Y_lm(cos_theta,phi)=rtheta_lm(cos_theta)* exp(i*m*phi)
*
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      double complex function YSpherical_lm(l,m,cos_theta,phi)
      implicit none
      integer l,m
      real*8 cos_theta,phi

*     **** external functions
      real*8   ytheta_lm
      external ytheta_lm

      YSpherical_lm = ytheta_lm(l,m,cos_theta)
     >               *dcmplx(dcos(m*phi),dsin(m*phi))
      return
      end





*     **************************************
*     *                                    *
*     *           dTesseral_lm             *
*     *                                    *
*     **************************************
*
*  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such that
*
*                                               {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                               {sin(|m|*phi)   m<0
*
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine dTesseral_lm(l,m,cos_theta,phi,dTx,dTy,dTz)
      implicit none
      integer l,m
      real*8 cos_theta,phi
      real*8 dTx,dTy,dTz

*     **** local variables ****
      integer mod_m
      real*8 tmp2,r,f1,f2

*     **** external functions
      real*8   rtheta_lm_div,drtheta_lm
      external rtheta_lm_div,drtheta_lm

      mod_m = abs(m)
      if (m.lt.0) then
         tmp2 = dsin(mod_m*phi)
         f2   = mod_m *dcos(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
      else if (m.gt.0) then
         tmp2 = dcos(mod_m*phi)
         f2   = -mod_m*dsin(mod_m*phi)*rtheta_lm_div(l,m,cos_theta)
      else
         tmp2 = 1.0d0
         f2   = 0.0d0
      end if
      f1 = drtheta_lm(l,m,cos_theta)*tmp2

      dTx = f1*cos_theta*dcos(phi) - f2*dsin(phi)
      dTy = f1*cos_theta*dsin(phi) + f2*dcos(phi)
      dTz = -f1*dsqrt(1.0d0-cos_theta**2)

      return
      end


*     **************************************
*     *                                    *
*     *           dYspherical_lm           *
*     *                                    *
*     **************************************
*
*  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such that
*
*      Y_lm(cos_theta,phi)=theta_lm(cos_theta)* exp(i*m*phi)
*
*    where   cos_theta = z/r and phi = atan2(y,x)
*
      subroutine dYspherical_lm(l,m,cos_theta,phi,dYx,dYy,dYz)
      implicit none
      integer l,m
      real*8 cos_theta,phi
      complex*16 dYx,dYy,dYz

*     **** local variables ****
      integer mod_m
      complex*16 f1,f2

*     **** external functions
      real*8   dytheta_lm,ytheta_lm_div
      external dytheta_lm,ytheta_lm_div

      !*** derivative with respect to theta ***
      f1 = dytheta_lm(l,m,cos_theta)*dcmplx(dcos(m*phi),dsin(m*phi))

      !*** derivative with respect to phi ***
      if(m.eq.0) then
        f2 = dcmplx(0.0d0,0.0d0)
      else
        f2 = ytheta_lm_div(l,m,cos_theta)
     >       *dcmplx(dcos(m*phi),dsin(m*phi))
     >       *m*dcmplx(0.0d0,1.0d0)
      end if

      dYx = f1*cos_theta*dcos(phi) - f2*dsin(phi)
      dYy = f1*cos_theta*dsin(phi) + f2*dcos(phi)
      dYz = -f1*dsqrt(1.0d0-cos_theta**2)

      return
      end




c!  *************************************************
c!  
c!     Name    : rtheta_lm
c!  
c!  
c!     Purpose : calculates rtheta_lm for a scalar cos_theta
c!               such that
c!                                                  {cos(|m|*phi)   m>0
c!         T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
c!                                                  {sin(|m|*phi)   m<0
c!  
c!  
c!  *************************************************
      double precision function rtheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta
c      double precision theta_lm !*RESULT*

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f

      !*** external functions *** 
      double precision rlegendre_lm
      external         rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l) 
     >write(*,*) 'parameter out of order in function rtheta_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then 
         coeff= 0.5d0
      else
         coeff= 1.0d0 
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)
      f     = coeff*rlegendre_lm(l,mod_m,cos_theta)

      rtheta_lm = f
      return
      end 

c!  *************************************************
c!  
c!     Name    : rtheta_lm_div
c!  
c!  
c!     Purpose : calculates rtheta_lm_div for a scalar cos_theta
c!               such that
c!                                                            {cos(|m|*phi)   m>0
c!     T_lm(cos_theta,phi)/sin_theta=rtheta_lm_div(cos_theta)*{1              m==0
c!                                                            {sin(|m|*phi)   m<0
c!  
c!  *************************************************

      double precision function rtheta_lm_div(l,m,cos_theta)        
      implicit none
      integer l,m                                               
      double precision cos_theta
c      double precision theta_lm !*RESULT*

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f
      
      !*** external functions ***
      double precision rlegendre_lm_div
      external         rlegendre_lm_div
      
      twopi = 8.0d0*datan(1.0d0)

      mod_m = abs(m)
      if (mod_m .gt. l)
     >write(*,*) 'parameter out of order in function rtheta_lm_div'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else 
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)
      f     = coeff*rlegendre_lm_div(l,mod_m,cos_theta)

      rtheta_lm_div = f
      return
      end


c!  *************************************************
c!  
c!     Name    : drtheta_lm
c!  
c!  
c!     Purpose : calculates drtheta_lm for a scalar cos_theta
c!               such that
c!                                                              {cos(|m|*phi)   m>0
c!         dT_lm(cos_theta,phi)/dtheta  = drtheta_lm(cos_theta)*{1              m==0
c!                                                              {sin(|m|*phi)   m<0
c!  
c!  
c!  *************************************************
      double precision function drtheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f,sin_theta

      !*** external functions *** 
      double precision rlegendre_lm
      external         rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)
      sin_theta = dsqrt(1.0d0-cos_theta*cos_theta)

      mod_m = abs(m)
      if (mod_m.gt.l) 
     >write(*,*) 'parameter out of order in function rtheta_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then 
         coeff= 0.5d0
      else
         coeff= 1.0d0 
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)
      if (sin_theta.le.1.0d-9) then
         if (mod_m.eq.1) then
            if (mod(l,2).eq.0) then
               f = 0.5d0*coeff*l*(l+1)
            else
               f = cos_theta*0.5d0*coeff*l*(l+1)
            end if
         else
            f = 0.0d0
         end if
      else
         if (mod_m.lt.l) then
            f = coeff*(mod_m*(cos_theta/sin_theta)
     >                      *rlegendre_lm(l,mod_m,cos_theta)
     >                - rlegendre_lm(l,mod_m+1,cos_theta))
         else
            f = coeff*(mod_m*(cos_theta/sin_theta)
     >                      *rlegendre_lm(l,mod_m,cos_theta))
         end if
      end if
      drtheta_lm = f
      return
      end 


c!  *************************************************
c! 
c!     Name    : drtheta_cos_theta_lm
c! 
c! 
c!     Purpose : calculates drtheta_cos_theta_lm for a scalar cos_theta
c!               such that
c!                                                                             {cos(|m|*phi)   m>0
c!         dT_lm(cos_theta,phi)/d_cos_theta  = drtheta_cos_theta_lm(cos_theta)*{1              m==0
c!                                                                             {sin(|m|*phi)   m<0
c! 
c! 
c!  *************************************************
      double precision function drtheta_cos_theta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f,sin_theta

      !*** external functions ***
      double precision rlegendre_lm
      external         rlegendre_lm

      twopi = 8.0d0*datan(1.0d0)
      sin_theta = dsqrt(1.0d0-cos_theta*cos_theta)

      mod_m = abs(m)
      if (mod_m.gt.l)
     >write(*,*) 'parameter out of order in function rtheta_lm'

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 0.5d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/twopi
      coeff = dsqrt(coeff)
      if (sin_theta.le.1.0d-9) then
         f = 0.0d0
c         if (mod_m.eq.1) then
c            if (mod(l,2).eq.0) then
c               f = 0.5d0*coeff*l*(l+1)
c            else
c               f = cos_theta*0.5d0*coeff*l*(l+1)
c            end if
c         else
c            f = 0.0d0
c         end if
      else
         if (mod_m.lt.l) then
            f = -coeff*(mod_m*(cos_theta/sin_theta**2)
     >                      *rlegendre_lm(l,mod_m,cos_theta)
     >                - rlegendre_lm(l,mod_m+1,cos_theta))
         else
            f = -coeff*(mod_m*(cos_theta/sin_theta**2)
     >                      *rlegendre_lm(l,mod_m,cos_theta))
         end if
      end if
      drtheta_cos_theta_lm = f
      return
      end




c!  *************************************************
c!  
c!     Name    : ytheta_lm
c!  
c!  
c!     Purpose : calculates ytheta_lm for a scalar cos_theta
c!               such that
c!                                                 
c!         Y_lm(cos_theta,phi)=ytheta_lm(cos_theta)*exp(i*m*phi)
c!                                                
c!  
c!  *************************************************
      double precision function ytheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta
c      double precision theta_lm !*RESULT*

      !*** local variables ***
      integer i,mod_m
      double precision coeff,twopi,fourpi,f

      !*** external functions *** 
      double precision legendre_lm
      external         legendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      mod_m = abs(m) 
      if (mod_m .gt. l) then
         write(*,*) 'parameter out of order in function ytheta_lm',l,m
      end if

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 1.0d0
      else
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif 
      coeff = coeff*(2*l+1)/fourpi
      coeff = dsqrt(coeff)
      f     = coeff*legendre_lm(l,mod_m,cos_theta)
      if ((m.lt.0).and.(mod(mod_m,2).eq.1)) f=-f

      ytheta_lm = f
      return
      end


!**************************************************
!
!       Name: dytheta_lm
!
!       Purpose:
!
!       Created:        8/05/2002
!**************************************************
      double precision function dytheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision  cos_theta

      !*** local variables ***
      integer i
      integer mod_m
      double precision  coeff,tmp,fourpi

      !*** external functions ***
      double precision  dlegendre_lm
      external          dlegendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      if( m .gt.l ) then
      call errquit("parameter out of order in function dytheta_lm",0,1)
      end if

      mod_m = abs(m)

      !*** find coefficient ***
      if (mod_m.eq.0) then
        coeff = 1.0d0
      else if(mod_m.gt.0) then
        coeff = 1.0d0
        do i=1,2*mod_m
          coeff = coeff/dble(l-mod_m+i)
        end do
      end if

      coeff = coeff*(2*l+1)/(fourpi)
      coeff = dsqrt(coeff)

      tmp = coeff*dlegendre_lm(l,mod_m,cos_theta)

      if (m.lt.0) then
        if( mod(mod_m,2).eq.1) then
           tmp = -tmp
        end if
      end if

      dytheta_lm = tmp
      return
      end


c  !*************************************************
c  !
c  !   Name    : ytheta_lm_div
c  !
c  !
c  !   Purpose : calculates ytheta_lm/sin(theta)
c  !
c  !
c  !   Created : 10/99
c  !
c  !*************************************************
      double precision function ytheta_lm_div(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,fourpi,f

      !*** external functions ***
      double precision legendre_lm_div
      external         legendre_lm_div

      fourpi = 16.0d0*datan(1.0d0)

      if ( m.gt.l) then
        call errquit(
     >   'parameter out of order in function theta_lm_div',0,0)
      end if

      mod_m = abs(m)

      !*** find coefficient ***
      if (mod_m.eq.0) then
        coeff = 1.0d0
      else if (mod_m.gt.0) then
        coeff = 1.0d0
        do i=1,2*mod_m
          coeff=coeff/dble(l-mod_m+i)
        end do
      end if

      coeff = coeff*(2*l+1)/fourpi
      coeff = dsqrt(coeff)

      f = coeff*legendre_lm_div(l,mod_m,cos_theta)

      if (m.lt.0) then
        if (mod(mod_m,2).eq.1) then
           f = -f
        end if
      end if

      ytheta_lm_div = f
      return
      end




*     ***********************************
*     *                                 *
*     *           rlegendre_lm          *
*     *                                 *
*     ***********************************
*
* Compute the associated legendre polynomial w/o 
* a Condon-Shortley phase.
*
      double precision function rlegendre_lm(l,m,x)
      implicit none
      integer l
      integer m
      double precision x
 
!     *** local variables ***
      integer i,mod_m
      double precision p_mm,p_mp1m,tmp,fact,coeff
 
      tmp = 0.0d0
      mod_m = abs(m)
      if (abs(x).gt.1) write(*,*)'parameter ot of range in legendre_lm'
 
      if (mod_m.gt.l) then
         rlegendre_lm = 0.0d0
         return 
      end if
 
!    *** find P(mod_m,mod_m) for mod_m=0 ***
      p_mm = 1.0d0
 
!    *** find P(mod_m,mod_m) for mod_m > 0 ***
      if (mod_m.gt.0) then
         tmp  = dsqrt((1.0d0-x)*(1.0d0+x))
         fact = 1.0d0
         do i = 1, mod_m
            p_mm = p_mm*fact*tmp
            fact = fact + 2.0d0
         end do
      end if
 
!    *** find P(l,mod_m) ***
      if (mod_m.eq.l) then
         rlegendre_lm = p_mm
      else
      !*** find P(mod_m+1,mod_m) ***
         p_mp1m = x*(2*mod_m + 1)*p_mm
         if (l.eq.mod_m+1) then
            rlegendre_lm = p_mp1m
         else
            do i = mod_m + 2, l
               tmp=(X*(2*i-1)*p_mp1m-(i+mod_m-1)*p_mm)/dble(i-mod_m)
               p_mm = p_mp1m
               p_mp1m = tmp
            end do
            rlegendre_lm = tmp
         end if
      end if
c      LEGENDRE_LM = LEGENDRE_LM*(-1)**MOD_M
 
c should not be needed 
cc!    *** negative m - this routine is only call with negative m from dtheta_lm and ddtheta_lm ***
c      if (m .LT. 0) then
c         coeff = 1.0d0
c         do i=1,2*mod_m
c            coeff = coeff/(l-mod_m+i)
c         END DO
c         rlegendre_lm = rlegendre_lm*coeff
c      end if
      return 
      end 

*     ***********************************
*     *                                 *
*     *           rlegendre_lm_div      *
*     *                                 *
*     ***********************************
*
* Compute the associated legendre polynomial divided by sin(theta) w/o 
* a Condon-Shortley phase.
*
      double precision function rlegendre_lm_div(l,m,x)
      implicit none
      integer l,m
      double precision x

      !*** local variables ***
      integer i
      double precision p_mm,p_mp1m,tmp,fact,f

      !*** check the arguments ***
      if ((m.lt.0)           .or.
     >    (m.gt.l)           .or.
     >    (dabs(x).gt.1.0d0) .or.
     >    (m.eq.0) ) then
        call errquit(
     >        'parameter ot of range in rlegendre_lm_div',0,0)
      end if

      !*** P(1,1,x)/dsqrt(1-x**2) ***
      p_mm = 1.0d0

      !*** P(m,m,x)/dsqrt(1-x**2)            ***
      !***   = (2*m-1)*dsqrt(1-x**2)         ***
      !***     *(P(m-1,m-1,x)/dsqrt(1-x**2)) ***
      fact = 3.0d0
      tmp = dsqrt((1.0d0-x)*(1.0d0+x))

      do i=2,m
        p_mm = p_mm*fact*tmp
        fact = fact + 2.0d0
      end do

      !*** find P(l,m) ***
      if (m.eq.l) then
        f = p_mm
      else
          !*** find P(m+1,m) ***
          p_mp1m = x*(2*m+1)*p_mm
          if (l.eq.(m+1)) then
            f = p_mp1m
          else
            do i=m+2,l
              tmp = (x*(2*i-1)*p_mp1m - (i+m-1)*p_mm)/dble(i-m)
              p_mm = p_mp1m
              p_mp1m = tmp
            end do
            f = tmp
          end if
      end if
      rlegendre_lm_div = f
      return
      end


*     *********************************************
*     *                                           *
*     *               legendre_lm                 *
*     *                                           *
*     *********************************************
*
* Compute the associated legendre polynomial w/
* a Condon-Shortley phase?
*
      DOUBLE PRECISION FUNCTION LEGENDRE_LM(L, M, X)
      INTEGER L
      INTEGER M
      DOUBLE PRECISION X
 
!    *** local variables ***
      INTEGER I
      INTEGER MOD_M
      DOUBLE PRECISION P_MM
      DOUBLE PRECISION P_MP1M
      DOUBLE PRECISION TMP
      DOUBLE PRECISION FACT
      DOUBLE PRECISION COEFF
 
      TMP=0.0d0
      MOD_M = ABS(M)
      IF (ABS(X) .GT. 1) CALL errquit(
     1   'parameter ot of range in legendre_lm',0,0)
 
      IF (MOD_M .GT. L) THEN
         LEGENDRE_LM = 0.0d0
         RETURN 
      ENDIF
 
!    *** find P(mod_m,mod_m) for mod_m=0 ***
      P_MM = 1.0D0
!    *** find P(mod_m,mod_m) for mod_m > 0 ***
      IF (MOD_M .GT. 0) THEN
         TMP = SQRT((1.0D0 - X)*(1.0D0 + X))
         FACT = 1.0D0
         DO I = 1, MOD_M
            P_MM = -P_MM*FACT*TMP
            FACT = FACT + 2.0D0
         END DO
      ENDIF
 
!    *** find P(l,mod_m) ***
      IF (MOD_M .EQ. L) THEN
         LEGENDRE_LM = P_MM
      ELSE
 
      !*** find P(mod_m+1,mod_m) ***
         P_MP1M = X*(2*MOD_M + 1)*P_MM
         IF (L .EQ. MOD_M+1) THEN
            LEGENDRE_LM = P_MP1M
         ELSE
            DO I = MOD_M + 2, L
               TMP=(X*(2*I-1)*P_MP1M-(I+MOD_M-1)*P_MM)/dble(I-MOD_M)
               P_MM = P_MP1M
               P_MP1M = TMP
            END DO
            LEGENDRE_LM = TMP
         ENDIF
      ENDIF
 
c!    *** negative m - this routine is only call with negative m from dtheta_lm and ddtheta_lm ***
      IF (M .LT. 0) THEN
         COEFF = 1.0D0
         DO I = 1, 2*MOD_M
            COEFF = COEFF/dble(L - MOD_M + I)
         END DO
         COEFF = COEFF*(-1)**MOD_M
         LEGENDRE_LM = LEGENDRE_LM*COEFF
      ENDIF
      RETURN 
      END 


c !*************************************************
c !
c !   Name    : dlegendre_lm
c !
c !
c !   Purpose : calculates the derivative of the
c !             associated Legendre with respect to angle theta
c !             (x=cos(theta))
c !
c !
c !   Created : 10/99
c !
c !
c !   Notes: Arfken definition of P(l,m,x) is different from mine by a
c !          factor (-1)**m. Thus all recursive relationships should be modified
c !          accordingly
c !*************************************************
      double precision function dlegendre_lm(l,m,x)
      implicit none
      integer l,m
      double precision  x
c      real*8  dlegendre_lm

      !*** local variables ***
      double precision f1,f2

      !*** external functions ****
      double precision legendre_lm
      external         legendre_lm

      !*** check the arguments ***
      if ((m.gt.l) .or. (dabs(x).gt.1.0d0)) then
        call errquit("parameter ot of range in dlegendre_lm",0,0)
      end if

      if (m.eq.l) then
        f1 = 0.0d0
      else
        f1 = legendre_lm(l,m+1,x)
      END IF

      f2 = legendre_lm(l,m-1,x)

      !*** see Arfken p.669 3rd ed, p. 725 4th ed.                             ***
      !*** sqrt(1-x*x)*dP_l^m/dx = 0.5*P_l^(m+1) - 0.5*(l+m)*(l-m+1)*P_l^(m-1) ***
      !*** Note this recursion relationship doesn't have the problems at       ***
      !*** theta=0,pi that the more traditional recursion relationship in term ***
      !*** of l (i.e. p. 115 of Jahnke and Emde) have                          ***
      !*** (1-x*x)*dP_l^m/dx = (l+1)*x*P_l^m - (l-m+1)*P_(l+1)^m               ***

      dlegendre_lm = 0.5d0*f1 - 0.5d0*(l+m)*(l-m+1)*f2

      return
      end


c  !*************************************************
c  !
c  !   Name    : legendre_lm_div
c  !
c  !   Purpose : calculates associated Legendre
c  !             polynomial divided by sin(theta)
c  !             for a scalar argument
c  !
c  !   Created : 10/99
c  !
c  !*************************************************
      double precision function legendre_lm_div(l,m,x)
      implicit none
      integer l,m
      double precision x
c      real*8 legendre_lm_div !*RESULT*

      !*** local variables ***
      integer i
      double precision p_mm,p_mp1m,tmp,fact,f

      !*** check the arguments ***
      if ((m.lt.0)           .or.
     >    (m.gt.l)           .or.
     >    (dabs(x).gt.1.0d0) .or.
     >    (m.eq.0) ) then
        call errquit(
     >        'parameter ot of range in legendre_lm_div',0,0)
      end if

      !*** P(1,1,x)/sqrt(1-x**2) ***
      p_mm = -1.0d0

      fact = 3.0d0
      tmp = dsqrt((1.0d0-x)*(1.0d0+x))

      do i=2,m
        p_mm = -p_mm*fact*tmp
        fact = fact + 2.0d0
      end do

      !*** find P(l,m) ***
      if (m.eq.l) then
        f = p_mm

      else
          !*** find P(m+1,m) ***
          p_mp1m = x*(2*m+1)*p_mm
          if (l.eq.(m+1)) then
            f = p_mp1m
          else
            do i=m+2,l
              tmp = (x*(2*i-1)*p_mp1m - (i+m-1)*p_mm)/(i-m)
              p_mm = p_mp1m
              p_mp1m = tmp
            end do
            f = tmp
          end if
      end if

      legendre_lm_div = f
      return
      end



c!  *************************************************
c!  
c!     Name    : theta_lm
c!  
c!  
c!     Purpose : calculates theta_lm for a scalar cos_theta
c!               such that
c!  
c!         Y_lm(cos_theta,phi)=theta_lm(cos_theta)*exp(i*m*phi)
c!  
c!  
c!     Created : 1/99
c!  
c!  *************************************************

      double precision function theta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta

      !*** local variables ***
      integer i,mod_m
      double precision coeff,fourpi,f

      !*** external functions *** 
      double precision legendre_lm
      external         legendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      if (m .gt. l) 
     > call errquit(
     >      'parameter out of order in function theta_lm',0,0)

      mod_m = abs(m)

!     *** find coefficient ***
      if (mod_m.eq.0) then
         coeff= 1.0d0
      else if (mod_m.gt.0) then
         coeff= 1.0d0
         do i=1,2*mod_m
            coeff = coeff/dble(l-mod_m+i)
         end do
      endif
      coeff = coeff*(2*l+1)/fourpi
      coeff = dsqrt(coeff)
      f     = coeff*legendre_lm(l,mod_m,cos_theta)
      if (m.lt.0) then
         if (mod(mod_m,2).eq.1) f=-f
      endif

      theta_lm = f
      return
      end 

