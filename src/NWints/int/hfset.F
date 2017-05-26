      Subroutine hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &                 Bxyz,Bprims,Bcoefs,NPB,NCB,
     &                 GENCON,alpha,ipair,ES,NPP)
! $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical GENCON

      Parameter (PI=3.1415926535898d0)
      Parameter (EXPLIM=100.d0)
#include "apiP.fh"

!--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

!--> Exponents, Pair Index & Prefactors for 2-ctr Overlap Distributions

      Dimension alpha(2,(NPA*NPB)),ipair(2,(NPA*NPB)),ES(3,(NPA*NPB))
!
! Compute the prefactors of the overlap distributions formed by the product of 
! two primitive Gaussians. These prefactors are defined as
!
!   ES = ESx * ESy * ESz, where
!
!
!          /  PI   \ 1/2      /   a b     2 \
!   ESx = | ------- |     EXP| - -----  Rx   |
!          \ a + b /          \  a + b      /
!
!
! N.B. 1) Overlap distributions with prefactors less than a given tolerance
!         are removed from the list. This shortened list is of length "NPP".
!      2) For segmented contractions, the product of contraction coefficients 
!         is also incorporated in the prefactor.
!
!******************************************************************************
 
      ! Python: 1.5*math.log(math.pi) = 1.7170948287741004
      parameter (LOGPIx15=1.7170948287741d0)
      parameter (N3OVER2=-1.5d0)
      double precision eps_small
      parameter (eps_small=1d-32)
      double precision const
      if(val_int_acc.lt.eps_small)  then
         const = LOGPIx15 - log(eps_small)
      else
         const = LOGPIx15 - log(val_int_acc)
      endif
      

      Rx2 = (Axyz(1) - Bxyz(1))**2
      Ry2 = (Axyz(2) - Bxyz(2))**2
      Rz2 = (Axyz(3) - Bxyz(3))**2

#if 0

      ! ORIGINAL CODE BELOW

      m2 = 0
      do mpa = 1,NPA
        a = Aprims(mpa)
        do mpb = 1,NPB
          b = Bprims(mpb)
         
          abi = 1.0d0/(a+b)
          beta = a*b*abi
         
          if ( ( beta*(Rx2+Ry2+Rz2) ) .lt.
     &         ( const + N3OVER2 * log(a+b) ) ) then
         
            m2 = m2 + 1
            
            alpha(1,m2) = a
            alpha(2,m2) = b
            
            ipair(1,m2) = mpa
            ipair(2,m2) = mpb
            
            s = sqrt(PI*abi)

            ESx = s*exp(-min(EXPLIM,beta*Rx2))
            ESy = s*exp(-min(EXPLIM,beta*Ry2))
            ESz = s*exp(-min(EXPLIM,beta*Rz2))
            
            ! GENCON is an easy branch to predict but should
            ! still be hoisted out of the loops.
            if( GENCON )then
              ES(1,m2) = ESx
              ES(2,m2) = ESy
              ES(3,m2) = ESz
            else ! GENCON
              ES(1,m2) = ESx*(Acoefs(mpa,1)*Bcoefs(mpb,1))
              ES(2,m2) = ESy
              ES(3,m2) = ESz
            end if ! GENCON
         
          end if ! threshold
        
        enddo ! mpb
      enddo ! mpa
      NPP = m2

#else

      m2 = 0
      do mpa = 1,NPA
        a = Aprims(mpa)
        do mpb = 1,NPB
          b = Bprims(mpb)
         
          abi = 1.0d0/(a+b)
          beta = a*b*abi
          if ( ( beta*(Rx2+Ry2+Rz2) ) .lt.
     &         ( const + N3OVER2 * log(a+b) ) ) then
         
            m2 = m2 + 1
            
            alpha(1,m2) = a
            alpha(2,m2) = b
            
            ipair(1,m2) = mpa
            ipair(2,m2) = mpb
            
          end if ! threshold
        enddo ! mpb
      enddo ! mpa
      NPP = m2

      if (GENCON) then
        do m2 = 1, NPP
       
          a = alpha(1,m2)
          b = alpha(2,m2)
       
          abi = 1.0d0/(a+b)
          beta = a*b*abi
          
          s = sqrt(PI*abi)
          
          ES(1,m2) = s*exp(-min(EXPLIM,beta*Rx2))
          ES(2,m2) = s*exp(-min(EXPLIM,beta*Ry2))
          ES(3,m2) = s*exp(-min(EXPLIM,beta*Rz2))
           
        enddo ! m2
      else ! GENCON
        do m2 = 1, NPP
       
          a = alpha(1,m2)
          b = alpha(2,m2)
       
          abi = 1.0d0/(a+b)
          beta = a*b*abi
          
          mpa = ipair(1,m2)
          mpb = ipair(2,m2)

          s = sqrt(PI*abi)
       
          c = Acoefs(mpa,1)
          d = Bcoefs(mpb,1)

          ESx = s*exp(-min(EXPLIM,beta*Rx2))
          ESy = s*exp(-min(EXPLIM,beta*Ry2))
          ESz = s*exp(-min(EXPLIM,beta*Rz2))
          
          !ES(1,m2) = ESx*(Acoefs(mpa,1)*Bcoefs(mpb,1))
          ES(1,m2) = ESx*c*d
          ES(2,m2) = ESy
          ES(3,m2) = ESz
           
        enddo ! m2
      endif ! GENCON

#endif

!      write(6,*)'-----start------ pair matrix '
!      do i=1,NPP
!        write(6,*)i,ipair(1,i),ipair(2,i)
!      enddo
!      write(6,*)'----- end ------ pair matrix '
      end
