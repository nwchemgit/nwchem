!  *******************************************************************
!  
!     Name    : HyperGeometric1F1(n,l,z)
!  
!  
!     Purpose : Calculates a special case of Kummer confluent
!               hypergeometric function, M(n+1/2,l+3/2,z) for z .LE. 0
!  
!     Created : Marat Valiev
!  
!  ********************************************************************
      DOUBLE PRECISION FUNCTION PAW_HYPERGEOM1F1 (N, L, Z)
      implicit none
 
      INTEGER N, L
      DOUBLE PRECISION Z
 
!    *** local variables ***
      INTEGER I
      DOUBLE PRECISION A, B, X
      DOUBLE PRECISION M1,M3
      DOUBLE PRECISION EPS
      DOUBLE PRECISION S
 
      DOUBLE PRECISION PAW_GAMMA,PAW_GAMMP
      EXTERNAL         PAW_GAMMA,PAW_GAMMP
 
!    *** cannot handle positive z ***
      IF (Z .GT. 0) THEN
         PRINT *, 'Invalid parameter'
         STOP 
      ENDIF
 
 
!    *** solution for z==0 ***
      IF (Z .EQ. 0.0D0) THEN
         PAW_HYPERGEOM1F1 = 1.0D0
         RETURN 
      ENDIF
 
      X = dsqrt(dabs(Z))
 
      IF (N .EQ. L) THEN
 
         PAW_HYPERGEOM1F1 = PAW_GAMMP(N + 0.5D0,(-Z))
     >                     *(N + 0.5)*(-Z)**((-N)- 0.5D0)
     >                     *PAW_GAMMA(N + 0.5D0)
 
         RETURN 
 
      ELSE IF (N .EQ. (L+1)) THEN
 
         PAW_HYPERGEOM1F1 = DEXP(Z)
 
         RETURN 
 
      ENDIF
 
!    *** do inifinite series for small z
      IF (dabs(Z) .LE. 1.0d0) THEN
 
         EPS = 1.0D-16
 
         PAW_HYPERGEOM1F1 = 1.0D0
         S = 1.0D0
 
         A = N + 0.5D0
         B = L + 1.5D0
         DO I = 1, 10000
 
            S = S*(A + I - 1)*Z/((B + I - 1)*I)
 
            PAW_HYPERGEOM1F1 = PAW_HYPERGEOM1F1 + S
 
            IF (dabs(S) .LT. EPS) RETURN 
 
         END DO
 
      CALL errquit("paw:cannot converge in 1f1",0,1)
 
         RETURN 
 
      ENDIF
 
      IF (N .LT. L) THEN
 
      !*** starting point n=l or b=a+1***
         A = N + 0.5D0
         B = N + 1.5D0
 
      !*** M1 = M(a,b-1) ***
      !*** M2 = M(a,b,z) ***
         M1 = dexp(Z)
         PAW_HYPERGEOM1F1 = PAW_GAMMP(A,(-Z))*A/(-Z)**A*PAW_GAMMA(A)
 
      !**********************************************
      ! using recursion formula
      ! z(a-b)M(a,b+1,z)=b(b-1)M(a,b-1,z)+b(1-b-z)M(a,b,z)
      ! obtain M(1/2,3/2+l  ,z) --> M2
      !        M(1/2,3/2+l-1,z) --> M2
      !**********************************************
         DO I = 1, L - N
 
            M3 = (B*(B - 1.0D0)*M1 + B*(1.0D0 - B - Z)*PAW_HYPERGEOM1F1
     1         )/(Z*(A - B))
 
            B = B + 1
            M1 = PAW_HYPERGEOM1F1
            PAW_HYPERGEOM1F1 = M3
 
         END DO
 
      ELSE IF (N .GT. (L+1)) THEN
 
      !*** starting point n=l+1 or b=a ***
         A = L + 1.5D0
         B = L + 1.5D0
 
      !*** M1 = M(a-1,b) ***
      !*** M2 = M(a,a,z) ***
         M1 = PAW_GAMMP(A - 1.0D0,(-Z))*(A - 1.0D0)/(-Z)**(A - 1.0D0)*
     1      PAW_GAMMA(A - 1.0D0)
         PAW_HYPERGEOM1F1 = dexp(Z)
 
      !**********************************************
      ! using recursion formula
      ! aM(a+1,b,z)=(b-a)M(a-1,b,z)+(2a-b+z)M(a,b,z)
      ! obtain M(n+1/2-1,3/2,z)   --> M1
      !        M(n+1/2  ,3/2,z)   --> M2
      !**********************************************
         DO I = 1, N - L - 1
 
            M3 = ((B - A)*M1 + (2*A - B + Z)*PAW_HYPERGEOM1F1)/A
 
            M1 = PAW_HYPERGEOM1F1
            PAW_HYPERGEOM1F1 = M3
            A = A + 1
 
         END DO
 
      ENDIF
 
      return
      end 
 
 
 
 
 
c $Id$
