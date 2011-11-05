c $Id$

!  *******************************************************************
!  
!     Name    : nwpw_hypergeometric1F1(n,l,z)
!  
!  
!     Purpose : Calculates a special case of Kummer confluent
!               hypergeometric function, M(n+1/2,l+3/2,z) for z .LE. 0
!  
!     Created : Marat Valiev
!  
!  ********************************************************************
      real*8 function nwpw_hypergeom1F1 (N, L, Z)
      implicit none
 
      integer N, L
      real*8 Z
 
!    *** local variables ***
      integer I
      real*8 A, B, X
      real*8 M1,M3
      real*8 EPS
      real*8 S
 
      real*8   util_gamma,util_gammp
      external util_gamma,util_gammp
 
!    *** cannot handle positive z ***
      IF (Z .GT. 0) THEN
         PRINT *, 'Invalid parameter'
         STOP 
      ENDIF
 
 
!    *** solution for z==0 ***
      IF (Z .EQ. 0.0D0) THEN
         nwpw_hypergeom1F1 = 1.0D0
         RETURN 
      ENDIF
 
      X = dsqrt(dabs(Z))
 
      IF (N .EQ. L) THEN
 
         nwpw_hypergeom1F1 = util_gammp(N + 0.5D0,(-Z))
     >                     *(N + 0.5)*(-Z)**((-N)- 0.5D0)
     >                     *util_gamma(N + 0.5D0)
 
         RETURN 
 
      ELSE IF (N .EQ. (L+1)) THEN
 
         nwpw_hypergeom1F1 = dexp(Z)
 
         RETURN 
 
      ENDIF
 
!    *** do inifinite series for small z
      IF (dabs(Z) .LE. 1.0d0) THEN
 
         EPS = 1.0D-16
 
         nwpw_hypergeom1F1 = 1.0D0
         S = 1.0D0
 
         A = N + 0.5D0
         B = L + 1.5D0
         DO I = 1, 10000
 
            S = S*(A + I - 1)*Z/((B + I - 1)*I)
 
            nwpw_hypergeom1F1 = nwpw_hypergeom1F1 + S
 
            IF (dabs(S) .LT. EPS) RETURN 
 
         END DO
 
         call errquit("paw:cannot converge in 1f1",0,1)
 
         RETURN 
 
      ENDIF
 
      IF (N .LT. L) THEN
 
      !*** starting point n=l or b=a+1***
         A = N + 0.5D0
         B = N + 1.5D0
 
      !*** M1 = M(a,b-1) ***
      !*** M2 = M(a,b,z) ***
         M1 = dexp(Z)
         nwpw_hypergeom1F1 = util_gammp(A,(-Z))*A/(-Z)**A*util_gamma(A)
 
      !**********************************************
      ! using recursion formula
      ! z(a-b)M(a,b+1,z)=b(b-1)M(a,b-1,z)+b(1-b-z)M(a,b,z)
      ! obtain M(1/2,3/2+l  ,z) --> M2
      !        M(1/2,3/2+l-1,z) --> M2
      !**********************************************
         DO I = 1, L - N
 
            M3 = (B*(B - 1.0D0)*M1 + B*(1.0D0 - B - Z)*nwpw_hypergeom1F1
     1         )/(Z*(A - B))
 
            B = B + 1
            M1 = nwpw_hypergeom1F1
            nwpw_hypergeom1F1 = M3
 
         END DO
 
      ELSE IF (N .GT. (L+1)) THEN
 
      !*** starting point n=l+1 or b=a ***
         A = L + 1.5D0
         B = L + 1.5D0
 
      !*** M1 = M(a-1,b) ***
      !*** M2 = M(a,a,z) ***
         M1 = util_gammp(A - 1.0D0,(-Z))*(A - 1.0D0)/(-Z)**(A - 1.0D0)*
     >      util_gamma(A - 1.0D0)
         nwpw_hypergeom1F1 = dexp(Z)
 
      !**********************************************
      ! using recursion formula
      ! aM(a+1,b,z)=(b-a)M(a-1,b,z)+(2a-b+z)M(a,b,z)
      ! obtain M(n+1/2-1,3/2,z)   --> M1
      !        M(n+1/2  ,3/2,z)   --> M2
      !**********************************************
         DO I = 1, N - L - 1
 
            M3 = ((B - A)*M1 + (2*A - B + Z)*nwpw_hypergeom1F1)/A
 
            M1 = nwpw_hypergeom1F1
            nwpw_hypergeom1F1 = M3
            A = A + 1
 
         END DO
 
      ENDIF
 
      return
      end 
 
 
 
 
 
