      Subroutine hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &                 Bxyz,Bprims,Bcoefs,NPB,NCB,
     &                 GENCON,alpha,ipair,ES,NPP)
! $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical GENCON

      Parameter (PI=3.1415926535898D0)
      Parameter (EXPLIM=100.D0)
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

      m1 = 0
      do mpa = 1,NPA
      do mpb = 1,NPB
       m1 = m1 + 1

       ipair(1,m1) = mpa
       ipair(2,m1) = mpb

      enddo
      enddo

      Rx2 = (Axyz(1) - Bxyz(1))**2
      Ry2 = (Axyz(2) - Bxyz(2))**2
      Rz2 = (Axyz(3) - Bxyz(3))**2

      m2 = 0
      do m1 = 1,(NPA*NPB)

       mpa = ipair(1,m1)
       mpb = ipair(2,m1)

       a = Aprims(mpa)
       b = Bprims(mpb)

       abi = 1/(a+b)
       beta = a*b*abi

       s = sqrt(PI*abi)

       ESx = s*exp(-min(EXPLIM,beta*Rx2))
       ESy = s*exp(-min(EXPLIM,beta*Ry2))
       ESz = s*exp(-min(EXPLIM,beta*Rz2))

       if( (ESx*ESy*ESz).gt.val_int_acc )then

        m2 = m2 + 1

        alpha(1,m2) = a
        alpha(2,m2) = b

        ipair(1,m2) = mpa
        ipair(2,m2) = mpb

        if( GENCON )then
         ES(1,m2) = ESx
         ES(2,m2) = ESy
         ES(3,m2) = ESz
        else ! GENCON
         ES(1,m2) = ESx*(Acoefs(mpa,1)*Bcoefs(mpb,1))
         ES(2,m2) = ESy
         ES(3,m2) = ESz
        end if ! GENCON

       end if ! val_int_acc

      enddo ! m1

      NPP = m2

!      write(6,*)'-----start------ pair matrix '
!      do i=1,NPP
!        write(6,*)i,ipair(1,i),ipair(2,i)
!      enddo
!      write(6,*)'----- end ------ pair matrix '
      end
