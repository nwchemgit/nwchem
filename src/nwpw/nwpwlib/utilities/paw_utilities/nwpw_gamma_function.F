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
      real*8   nwpw_gamma,nwpw_gammp
      external nwpw_gamma,nwpw_gammp
 
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
         nwpw_SpecialKummer = nwpw_gammp(n+0.5d0,(-z))
     >                       *(n+0.5d0) 
     >                       *((-z)**((-n)- 0.5d0))
     >                       *nwpw_gamma(n+0.5d0)
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
         nwpw_SpecialKummer = nwpw_gammp(a,(-z))*a/(-z)**a*nwpw_gamma(a)
 
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
         m1 = nwpw_gammp(a-1.0d0,(-z))*(a-1.0d0)/(-z)**(a-1.0d0)*
     >      nwpw_gamma(a-1.0d0)
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
 
*     ***************************************   
*     *                                     *
*     *          nwpw_ln_gamma              *
*     *                                     *
*     ***************************************   

*        This function computes the Log(Gamma) function
*       Serves as backup if intrinic Gamma function is not available
 
      real*8 function nwpw_ln_gamma(XX)
      implicit none
      real*8 XX
      integer J
      real*8 SER
      real*8 STP
      real*8 TMP
      real*8 X
      DOUBLE PRECISION Y, COF(6)
      SAVE STP, COF
      DATA COF, STP/ 76.18009172947146d0, -86.50532032941677d0, 
     1   24.01409824083091d0, -1.231739572450155d0, 
     2   0.1208650973866179d-2, -.5395239384953d-5, 2.5066282746310005d0
     3   / 
      X = XX
      Y = X
      TMP = X + 5.5d0
      TMP = (X + 0.5d0)*dlog(TMP) - TMP
      SER = 1.000000000190015d0
      DO J = 1, 6
         Y = Y + 1.0d0
         SER = SER + COF(J)/Y
      END DO
      nwpw_ln_gamma = TMP + dlog(STP*SER/X)
      return 
      end 
 
*     **********************************************
*     *                                            *
*     *            nwpw_gamma                      *
*     *                                            *
*     **********************************************
*
*     This function computes the Gamma function
*
      real*8 function nwpw_gamma(X)
      implicit none
      real*8 X
 
      real*8 XX
      real*8 nwpw_ln_gamma
      external nwpw_ln_gamma
 
      XX = X
      nwpw_gamma = dexp(nwpw_ln_gamma(XX))
      return
      end 
 
*     **********************************************
*     *                                            *
*     *           nwpw_gammp                       *
*     *                                            *
*     **********************************************
*
*     This function computes the incomplete gamma function P(a,x) = igamma(a,x)/Gamma(a)
*  
*                             /x
*       P(a,x) = 1/Gamma(a) * | exp(-t) * t**(a-1) dt
*                             /0
*
      real*8 function nwpw_gammp(A, X)
      implicit none
      real*8 a,x
      real*8 gammcf,gamser,gln
 
      if (x.lt. (a+1.0d0)) then
         call nwpw_gser(gamser, a, x, gln)
         nwpw_gammp = gamser
      else
         call nwpw_gcf(gammcf,a,x,gln)
         nwpw_gammp = 1.0d0 - gammcf
      end if
      return
      end 
 
*     **********************************************
*     *                                            *
*     *               nwpw_gcf                     *
*     *                                            *
*     **********************************************
      subroutine nwpw_gcf(GAMMCF, A, X, GLN)
      implicit none
      integer ITMAX
      real*8 A,GAMMCF,GLN,X,EPS,FPMIN
      parameter (ITMAX = 100, EPS = 3.0d-16, FPMIN = 1.0d-30)
      real*8 AN,B,C,D,DEL,H
      integer I

      real*8   nwpw_ln_gamma
      external nwpw_ln_gamma
 
      GLN = nwpw_ln_gamma(A)
      B = X + 1.0d0 - A
      C = 1.0d0/1.0d-30
      D = 1.0d0/B
      H = D
      DO I = 1, 100
         AN = -I*(I - A)
         B = B + 2.0d0
         D = AN*D + B
         IF (dabs(D) .LT. 1.0d-30) D = 1.0d-30
         C = B + AN/C
         IF (dabs(C) .LT. 1.0d-30) C = 1.0d-30
         D = 1.0d0/D
         DEL = D*C
         H = H*DEL
         IF (dabs(DEL-1.0d0) .LT. 3.0d-16) GO TO 1
      END DO
      call errquit('a too large, ITMAX too small in nwpw_gcf',0,0)

    1 continue
      GAMMCF = dexp((-X) + A*dlog(X) - GLN)*H
      return
      end 
 
*     **********************************************
*     *                                            *
*     *             nwpw_gser                      *
*     *                                            *
*     **********************************************
      subroutine nwpw_gser(GAMSER, A, X, GLN)
      implicit none
      real*8 A,X
      real*8 GAMSER,GLN
 
!    *** local variables ***
      integer ITMAX
      parameter (ITMAX = 100)
      real*8 EPS
      parameter (EPS = 3.0d-16)
      integer N
      real*8 AP,DEL,SUM

      real*8   nwpw_ln_gamma
      external nwpw_ln_gamma
 
      GLN = nwpw_ln_gamma(A)
 
      if (X .le. 0.0d0) then
         if (X.lt.0.0d0) call errquit("x < 0 in nwpw_gser",0,1)
         GAMSER = 0.0d0
         return
      end if
 
      AP = A
      SUM = 1.0d0/A
      DEL = SUM
      do N = 1, 100
         AP = AP + 1.0d0
         DEL = DEL*X/AP
         SUM = SUM + DEL
         if (dabs(DEL) .lt. dabs(SUM)*3.0d-16) GO TO 1
 
      end do
 
      call errquit
     >     ("a too large,ITMAX too small in nwpw_gser",0,1)
 
    1 continue
      GAMSER = SUM*dexp((-X) + A*dlog(X) - GLN)
 
      return 
      end 
 
c $Id$
