*----------------------------------------------------------------------*
      subroutine gl_weights(xa,xb,npnts,xrts,wght)
*----------------------------------------------------------------------*
*
* purpose: given the interval between xa and xb, this routine returns
*          npnts roots xrts(*) and weights wght(*) for
*          gauss-legendre quadrature
*
*  algorithm taken from Press et al. "Numerical Recepies"
*
*----------------------------------------------------------------------*

      implicit none

* input
      real*8, intent(in) ::
     &     xa, xb
      integer, intent(in) ::
     &     npnts
* output
      real*8, intent(out) ::
     &     xrts(npnts), wght(npnts)

* local:
* constants
      integer, parameter ::
     &     ntest = 00
      real*8, parameter ::
     &     eps = 3d-14, pi = 3.141592653589793238d0

* variables
      integer ::
     &     i, j, m
      real*8 ::
     &     xm, xl, z, z1, p1, p2, p3, pp

      m = (npnts+1)/2
      xm = 0.5d0*(xb + xa)
      xl = 0.5d0*(xb - xa)
      ! loop over roots
      do i = 1, m
        ! start guess
        z = cos(pi*(dble(i)-.25d0)/(npnts+0.5d0))
        ! Newton
        newton: do
          p1 = 1.0d0
          p2 = 0.0d0
          do j = 1, npnts
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*dble(j)-1.0d0) * z*p2 - (j-1.0d0)*p3)/dble(j)
          end do
          ! p1 is now the value of the Legendre polynomial at the
          ! current approximation of its root
          ! pp is its derivative at this point:
          pp = npnts*(z*p1-p2)/(z*z-1.0d0)
          z1 = z
          ! Newton step:
          z = z1-p1/pp
          if (abs(z-z1).lt.eps) exit newton
        end do newton
        xrts(i)         = xm-xl*z
        xrts(npnts+1-i) = xm+xl*z
        wght(i)         = 2d0*xl/((1.0d0 - z*z)*pp*pp)
        wght(npnts+1-i) = wght(i)
      end do
      
      if (ntest.ge.5) then
        write(6,*) 'Gauss-Legendre roots and weigths for '
        write(6,*) npnts,' points in [',xa,':',xb,']'
        write(6,*)
        write(6,*) 'roots:'
        call wrtmat(xrts,1,npnts,1,npnts) 
        write(6,*) 'weights:'
        call wrtmat(wght,1,npnts,1,npnts) 
      end if

      return
      end

*----------------------------------------------------------------------*
*----------------------------------------------------------------------
      subroutine s_weights(xa,xb,npnts,xrts,wght)
*----------------------------------------------------------------------*
*
* purpose: given the interval between xa and xb, this routine returns
*          npnts roots xrts(*) and weights wght(*) for the simpson rule
*
*----------------------------------------------------------------------*

      implicit none

* input
      real*8, intent(in) ::
     &     xa, xb
      integer, intent(in) ::
     &     npnts
* output
      real*8, intent(out) ::
     &     xrts(npnts), wght(npnts)

* local:
* constants:
      integer, parameter ::
     &     ntest = 00
* variables:
      real*8 ::
     &     h, delt
      integer ::
     &     ipnt


      if (mod(npnts,2).ne.1) then
        write(6,*) 'Simpson: number of quadrature points must be uneven'
        stop 'Simpson'
      end if

      ! equal differences      
      h = abs(xb - xa)/dble(npnts-1)
      delt = xb-xa

      xrts(1) = xa
      wght(1) = h/3d0
      do ipnt = 2, npnts-1
        xrts(ipnt) = xa + delt*dble(ipnt-1)/dble(npnts-1)
        if (mod(ipnt,2).eq.0) wght(ipnt)= (4d0/3d0)*h
        if (mod(ipnt,2).eq.1) wght(ipnt)= (2d0/3d0)*h
      end do
      xrts(npnts) = xb
      wght(npnts) = h/3d0

      if (ntest.ge.5) then
        write(6,*) 'Simpson roots and weigths for '
        write(6,*) npnts,' points in [',xa,':',xb,']'
        write(6,*)
        write(6,*) 'roots:'
        call wrtmat(xrts,1,npnts,1,npnts) 
        write(6,*) 'weights:'
        call wrtmat(wght,1,npnts,1,npnts)
 
      end if

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------
      subroutine test_quad(xa,xb,npnts,xrts,wght)
*----------------------------------------------------------------------*
*
* purpose: given the interval between xa and xb, this routine tests
*          the quadrature on some functions
*
*----------------------------------------------------------------------*

      implicit none

* input
      real*8, intent(in) ::
     &     xa, xb
      integer, intent(in) ::
     &     npnts
      real*8, intent(in) ::
     &     xrts(npnts), wght(npnts)

* local:
* constants:
      integer, parameter ::
     &     ntest = 00
      real*8, parameter ::
     &     pi = 3.141592653589793238d0
* variables:
      real*8 ::
     &     sum, fw
      integer ::
     &     ipnt

      ! test 1: a constant function
      sum = 0d0
      do ipnt = 1, npnts
        sum = sum + wght(ipnt)*1d0/(xb-xa)
        write(6,'(x,i4,4e15.6)')
     &       ipnt, wght(ipnt), xrts(ipnt), 1d0/(xb-xa), sum
      end do
      write (6,*) ' test 1 = ',sum
      write (6,*) '    ref = ',1d0
      ! test 2: a scaled sin function
      sum = 0d0
      do ipnt = 1, npnts
        fw  = wght(ipnt)*(xb-xa)*pi/2d0*sin(xrts(ipnt)/(xb-xa)*pi/2d0)
        sum = sum + fw
        write(6,'(x,i4,4e15.6)')
     &       ipnt, wght(ipnt), xrts(ipnt), fw , sum
      end do
      write (6,*) ' test 2 = ',sum
      write (6,*) '    ref = ',1d0

      return
      end
c $Id$
