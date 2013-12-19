c
c $Id$
c

*     ******************************************************
*     *                                                    *
*     *             nwpw_SpecialKummer                     *
*     *                                                    *
*     ******************************************************
*
*     Calculates a special case of the Kummer confluent hypergeometric 
*     function, M(n+1/2,l+3/2,z) for z .LE. 0
*
*     This function was created by  Marat Valiev, and  modified by Eric Bylaska.
*     See Abramowitz and Stegun for the formulas used in this function.
*
      real*8 function nwpw_SpecialKummer(n,l,z)
      implicit none
      integer n,l
      real*8  z
 
*     *** local variables ***
      real*8 eps
      parameter (eps=1.0d-16)

      integer i
      real*8 a,b,m1,m3,s

*     **** external functions ****
      real*8   util_gamma,util_gammp
      external util_gamma,util_gammp
 
      nwpw_SpecialKummer = 0.0d0

*    *** cannot handle positive z ***
      if (z.gt.0.0d0) then
         call errquit('nwpw_SpecialKummer:invalid parameter, z>0',0,1)
      end if
 
 
*    *** solution for z==0 ***
      if (z.eq.0.0d0) then
         nwpw_SpecialKummer = 1.0d0
         return 
      end if
 
*     ***** M(a,a+1,z) = a * (-z)**(-a) * igamma(a,-z) = a * (-z)**(-a) * P(a,-z) *Gamma(a)  where z is real and a = (n+0.5)  ****
      if (n.eq.l) then
         nwpw_SpecialKummer = util_gammp(n+0.5d0,(-z))
     >                       *(n+0.5d0) 
     >                       *((-z)**((-n)- 0.5d0))
     >                       *util_gamma(n+0.5d0)
         return 
 
*     ***** M(a,a,z) = exp(z)  where a = (n+0.5)  ****
      else if (n.eq.(l+1)) then
         nwpw_SpecialKummer = dexp(z)
         return 
      end if
 
!     *** do inifinite series for small z
      if (dabs(z).le.1.0d0) then
 
         nwpw_SpecialKummer = 1.0d0
         s = 1.0d0
         a = n + 0.5d0
         b = l + 1.5d0
         do i=1,10000
            s = s*(a+i-1)*z/((b+i-1)*i)
            nwpw_SpecialKummer = nwpw_SpecialKummer + s
            if (dabs(s).lt.eps) return 
         end do
         call errquit("nwpw_SpecialKummer:cannot converge",0,1)
         return 
      end if
 
      if (n.lt.l) then
 
      !*** starting point n=l or b=a+1***
         a = n + 0.5d0
         b = n + 1.5d0
 
      !*** m1 = M(a,b-1) ***
      !*** m2 = M(a,b,z) ***
         m1 = dexp(z)
         nwpw_SpecialKummer = util_gammp(a,(-z))*a/(-z)**a*util_gamma(a)
 
      !**********************************************
      ! using recursion formula
      ! z(a-b)M(a,b+1,z)=b(b-1)M(a,b-1,z)+b(1-b-z)M(a,b,z)
      ! obtain M(1/2,3/2+l  ,z) --> m2
      !        M(1/2,3/2+l-1,z) --> m2
      !**********************************************
         do i=1,l-n
            m3=(b*(b-1.0d0)*m1+b*(1.0d0-b-z)*nwpw_SpecialKummer)
     >         /(z*(a-b))
            b = b + 1
            m1 = nwpw_SpecialKummer
            nwpw_SpecialKummer = m3
         end do
 
      else if (n.gt.(l+1)) then
 
      !*** starting point n=l+1 or b=a ***
         a = l + 1.5d0
         b = l + 1.5d0
 
      !*** m1 = M(a-1,b) ***
      !*** m2 = M(a,a,z) ***
         m1 = util_gammp(a-1.0d0,(-z))*(a-1.0d0)/(-z)**(a-1.0d0)*
     >      util_gamma(a-1.0d0)
         nwpw_SpecialKummer = dexp(z)
 
      !**********************************************
      ! using recursion formula
      ! aM(a+1,b,z)=(b-a)M(a-1,b,z)+(2a-b+z)M(a,b,z)
      ! obtain M(n+1/2-1,3/2,z)   --> m1
      !        M(n+1/2  ,3/2,z)   --> m2
      !**********************************************
         do i=1,n-l-1
            m3 = ((b-a)*m1+(2*a-b+z)*nwpw_SpecialKummer)/a
            m1 = nwpw_SpecialKummer
            nwpw_SpecialKummer = m3
            a = a + 1
         end do
      end if
 
      return
      end 
