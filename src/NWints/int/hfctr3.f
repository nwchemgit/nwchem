      Subroutine hfctr3(Axyz,Bxyz,Cxyz,alpha,G,NABC)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension Axyz(3),Bxyz(3),Cxyz(3)

      Dimension alpha(4,NABC),G(3,NABC)
c
c Define the center of the charge distribution formed by the product of 
c three Gaussians.
c
c******************************************************************************

      Ax = Axyz(1)
      Ay = Axyz(2)
      Az = Axyz(3)

      Bx = Bxyz(1)
      By = Bxyz(2)
      Bz = Bxyz(3)

      Cx = Cxyz(1)
      Cy = Cxyz(2)
      Cz = Cxyz(3)

      do 10 m = 1,NABC

       A = alpha(1,m)
       B = alpha(2,m)
       C = alpha(3,m)

       ABCI = 1/(A + B + C)

       G(1,m) = ABCI*(A*Ax + B*Bx + C*Cx)
       G(2,m) = ABCI*(A*Ay + B*By + C*Cy)
       G(3,m) = ABCI*(A*Az + B*Bz + C*Cz)

   10 continue

      end
