      Subroutine hf1set3(Axyz,Aprims,Acoef,NPA,
     &                   Bxyz,Bprims,Bcoef,NPB,
     &                   Cxyz,Cprims,Ccoef,NPC,
     &                   alpha,NABC)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Parameter (PI=3.1415926535898D0)
      Parameter (EXPLIM=100.D0)
#include "apiP.fh"

      Dimension Axyz(3),Aprims(NPA),Acoef(NPA)
      Dimension Bxyz(3),Bprims(NPB),Bcoef(NPB)
      Dimension Cxyz(3),Cprims(NPC),Ccoef(NPC)

      Dimension alpha(4,NPA*NPB*NPC)
c
c The prefactor of the charge distribution, formed by the product of 
c three Gaussians, is defined as
c
c         /    PI     \ 3/2     /   a b   __2 \      /   (a + b) c   __2 \
c   ES = | ----------- |    EXP| - -----  AB   | EXP| - -----------  PC   |
c         \ a + b + c /         \  a + b      /      \   a + b + c       /
c
c N.B. 1) This prefactor is returned as the 4th index of the exponents array.
c         (i.e., "alpha(4,m)") 
c      2) Charge distributions, whose prefactor is less than a given tolerance,
c         are removed from the list. The shortened list is of length "NABC".
c      3) The product of contraction coefficients is also incorporated in
c         the prefactor. For a general contraction implementation, this 
c         should not be done at this point.
c
c******************************************************************************

      m = 0
      do 10 ipa = 1,NPA
      do 10 ipb = 1,NPB
      do 10 ipc = 1,NPC
       m = m + 1

       alpha(1,m) = Aprims(ipa)
       alpha(2,m) = Bprims(ipb)
       alpha(3,m) = Cprims(ipc)

       alpha(4,m) = Acoef(ipa)*Bcoef(ipb)*Ccoef(ipc)

   10 continue

      Ax = Axyz(1)
      Ay = Axyz(2)
      Az = Axyz(3)

      Bx = Bxyz(1)
      By = Bxyz(2)
      Bz = Bxyz(3)

      Cx = Cxyz(1)
      Cy = Cxyz(2)
      Cz = Cxyz(3)

      RABsq = (Ax-Bx)**2 + (Ay-By)**2 + (Az-Bz)**2

      m2 = 0
      do 20 m1 = 1,NPA*NPB*NPC

       A = alpha(1,m1)
       B = alpha(2,m1)
       C = alpha(3,m1)

       AB  = A + B
       ABI = 1/AB

       ES = exp(-min(EXPLIM,A*B*ABI*RABsq))

       Px = ABI*(A*Ax + B*Bx)
       Py = ABI*(A*Ay + B*By)
       Pz = ABI*(A*Az + B*Bz)

       ABCI = 1/(A + B + C)

       RPCsq = (Px-Cx)**2 + (Py-Cy)**2 + (Pz-Cz)**2

       ES = ES*exp(-min(EXPLIM,(AB*C*ABCI)*RPCsq))

       s = sqrt(PI*ABCI)

       ES = s*s*s*ES

       if( ES .gt. val_int_acc )then

        m2 = m2 + 1

        alpha(1,m2) = A
        alpha(2,m2) = B
        alpha(3,m2) = C

        alpha(4,m2) = ES*alpha(4,m1)

       end if

   20 continue
      NABC = m2

      end
