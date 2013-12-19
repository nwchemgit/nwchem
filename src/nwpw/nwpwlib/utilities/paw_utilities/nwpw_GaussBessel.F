c
c $Id$
c

*     ******************************************************
*     *                                                    *
*     *               nwpw_GaussBessel                     *
*     *                                                    *
*     ******************************************************
*
*     Calculates the Gaussian Bessel function, 
*
*                                 /infinity
*    GausBessel(n,l,alpha,R)  =   | k**n exp(-alpha**2 * k**2) j_l(R*k) dk
*                                 /0
*
*     This function uses the SpecialKummer function. Note it is assumed that
*     (n+l) is an even integer.
*
      real*8 function nwpw_GaussBessel(n,l,alpha,R)
      implicit none
      integer n,l
      real*8  alpha,R
 
*     *** local variables ***
      integer i
      real*8 c,pi
c      real*8 c,x,y,pi,c1

*     **** external functions ****
c      integer  nwpw_doublefactorial
c      external nwpw_doublefactorial
      real*8   nwpw_SpecialKummer,ln_gamma
      external nwpw_SpecialKummer,ln_gamma

      if (mod(n+l,2).eq.1) then
         call errquit("nwpw_GaussBessel: n+l is not even",0,0)
      end if

      pi = 4.0d0*datan(1.0d0)
c      x = dble(nwpw_doublefactorial(n+l-1))
c      y = dble(nwpw_doublefactorial(2*l+1))
c      c1 = (dsqrt(pi)/((2.0d0**((n+l)/2+1))*(alpha**(n+l+1)))) * (x/y)

      c = (dsqrt(pi)/((2.0d0**(l+2))*(alpha**(n+l+1)))) 
     >  * dexp(ln_gamma((n+l+1)/2.0d0) - ln_gamma(l+1.5d0))
      
      nwpw_GaussBessel =  c
     >               * (R**l) 
     >               * nwpw_SpecialKummer((n+l)/2,l,-(0.5d0*R/alpha)**2)
      return
      end 

*     ******************************************************
*     *                                                    *
*     *                nwpw_dGaussBessel                   *
*     *                                                    *
*     ******************************************************
*
*     Calculates the derivative of the Gaussian Bessel function wrt R,
*
*                                       /infinity
*   dGausBessel(n,l,alpha,R)  = (d/dR)  | k**n exp(-alpha**2 * k**2) j_l(R*k) dk
*                                       /0
*
*     This function uses the SpecialKummer function. Note it is assumed that
*     (n+l) is an even integer.
*
      real*8 function nwpw_dGaussBessel(n,l,alpha,R)
      implicit none
      integer n,l
      real*8  alpha,R

*     *** local variables ***
      integer i
      real*8 c,x,y,pi

*     **** external functions ****
      real*8   nwpw_SpecialKummer,ln_gamma
      external nwpw_SpecialKummer,ln_gamma

      if (mod(n+l,2).eq.1) then
         call errquit("nwpw_GaussBessel: n+l is not even",0,0)
      end if

      pi = 4.0d0*datan(1.0d0)
c      x = dble(nwpw_doublefactorial(n+l-1))
c      y = dble(nwpw_doublefactorial(2*l+1))
c      c = (dsqrt(pi)/((2.0d0**((n+l)/2+1))*(alpha**(n+l+1)))) * (x/y)

      c = (dsqrt(pi)/((2.0d0**(l+2))*(alpha**(n+l+1)))) 
     >  * dexp(ln_gamma((n+l+1)/2.0d0) - ln_gamma(l+1.5d0))

      nwpw_dGaussBessel =  c
     >    * ( (l*R**(l-1))
     >        *nwpw_SpecialKummer((n+l)/2,l,-(0.5d0*R/alpha)**2)
     >      - (0.5d0*R**(l+1)/alpha**2) 
     >        *dble((n+l)/2+0.5d0)/dble(l+1.5d0)
     >        *nwpw_SpecialKummer((n+l)/2+1,l+1,-(0.5d0*R/alpha)**2)
     >      )
      return
      end




*     ****************************************
*     *                                      *
*     *         nwpw_doublefactorial         *
*     *                                      *
*     ****************************************
*
*     This functions returns the double factorial n!! for n from -1 to 16      
*
      integer function nwpw_doublefactorial(n)
      implicit none
      integer n

      integer ndf(18)
      data ndf /    1,  
     >               1,
     >               1,
     >               2,
     >               3,
     >               8,
     >               15,
     >               48,
     >               105,
     >               384,
     >               945,
     >               3840,
     >               10395,
     >               46080,
     >               135135,
     >               645120,
     >               2027025,
     >               10321920/
      
      nwpw_doublefactorial = ndf(n+2)
      return
      end

