c
c $Id: shake_bonddiff.F 21176 2011-10-10 06:35:49Z d3y133 $
c

*     ********************************************
*     *                                          *
*     *           shake_bonddiff                 *
*     *                                          *
*     ********************************************
*
*   This routine handles the bond difference constraint.
*
*   sigma = |Ra-Rb|^2 - |Rc-Rd|^2 - gamma
*
*   Entry - indx(4) - indexes of the 4 atoms a,b,c,d
*           tol - tolerance of lagrange multiplier iteration
*           maxit - number iterations to solve for Lagrange multiplier
*           gamma - constraint value
*           mass(4) - masses of the atoms a,b,c,d
*           r2(3,*) - partially updated positions
*           r1(3,*) - previous positions
*           x - lagrange multiplier constraint
*
      subroutine shake_bonddiff(indx,
     >                       tol,maxit,
     >                       gamma,mass,
     >                       r2,r1,x)
      implicit none
      integer indx(4)
      real*8  tol
      integer maxit
      real*8  gamma,mass(4)
      real*8  r2(3,*), r1(3,*)
      real*8  x

      logical done
      integer it
      real*8 Rab(3),Rcd(3),tRab(3),tRcd(3),mabinv,mcdinv
      real*8 A,B,C,xold

      Rab(1) = r1(1,indx(1)) - r1(1,indx(2))
      Rab(2) = r1(2,indx(1)) - r1(2,indx(2))
      Rab(3) = r1(3,indx(1)) - r1(3,indx(2))

      Rcd(1) = r1(1,indx(3)) - r1(1,indx(4))
      Rcd(2) = r1(2,indx(3)) - r1(2,indx(4))
      Rcd(3) = r1(3,indx(3)) - r1(3,indx(4))

      tRab(1) = r2(1,indx(1)) - r2(1,indx(2))
      tRab(2) = r2(2,indx(1)) - r2(2,indx(2))
      tRab(3) = r2(3,indx(1)) - r2(3,indx(2))

      tRcd(1) = r2(1,indx(3)) - r2(1,indx(4))
      tRcd(2) = r2(2,indx(3)) - r2(2,indx(4))
      tRcd(3) = r2(3,indx(3)) - r2(3,indx(4))

      mabinv = 1.0d0/mass(1) + 1.0d0/mass(2)
      mcdinv = 1.0d0/mass(3) + 1.0d0/mass(4)

      A = (tRab(1)**2 + tRab(2)**2 + tRab(3)**2)
     >  - (tRcd(1)**2 + tRcd(2)**2 + tRcd(3)**2) - gamma
   
      B = 2.0d0*mabinv*(tRab(1)*Rab(1)+tRab(2)*Rab(2)+tRab(3)*Rab(3))
     >  - 2.0d0*mcdinv*(tRcd(1)*Rcd(1)+tRcd(2)*Rcd(2)+tRcd(3)*Rcd(3))
   
      C = (mabinv**2)*(Rab(1)**2+Rab(2)**2+Rab(3)**2)
     >  - (mabinv**2)*(Rab(1)**2+Rab(2)**2+Rab(3)**2)

      x = -C
      it = 0
      done = .false.
      do while (.not. done)
         xold = x
         x    = -C + (1.0d0-B)*x + A*x*x
         it = it + 1
         done = (dabs(x-xold).lt.tol).or.(it.gt.maxit)
      end do

      r2(1,indx(1)) = r2(1,indx(1)) + x*Rab(1)/mass(1)
      r2(2,indx(1)) = r2(2,indx(1)) + x*Rab(2)/mass(1)
      r2(3,indx(1)) = r2(3,indx(1)) + x*Rab(3)/mass(1)

      r2(1,indx(2)) = r2(1,indx(2)) - x*Rab(1)/mass(2)
      r2(2,indx(2)) = r2(2,indx(2)) - x*Rab(2)/mass(2)
      r2(3,indx(2)) = r2(3,indx(2)) - x*Rab(3)/mass(2)

      r2(1,indx(3)) = r2(1,indx(3)) + x*Rcd(1)/mass(3)
      r2(2,indx(3)) = r2(2,indx(3)) + x*Rcd(2)/mass(3)
      r2(3,indx(3)) = r2(3,indx(3)) + x*Rcd(3)/mass(3)

      r2(1,indx(4)) = r2(1,indx(4)) - x*Rcd(1)/mass(4)
      r2(2,indx(4)) = r2(2,indx(4)) - x*Rcd(2)/mass(4)
      r2(3,indx(4)) = r2(3,indx(4)) - x*Rcd(3)/mass(4)

      return
      end
