c
c $Id$
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
*   Entry - na - number of atoms (na==3) or (na==4)
*           indx(*) - indexes of the na atoms a,b,c,d
*           tol - tolerance of lagrange multiplier iteration
*           maxit - number iterations to solve for Lagrange multiplier
*           gamma - constraint value
*           mass(*) - masses of the atoms a,b,c,d
*           r2(3,*) - partially updated positions
*           r1(3,*) - previous positions
*           x - lagrange multiplier constraint
*
      subroutine shake_bonddiff(na,indx,
     >                       tol,maxit,
     >                       gamma,mass,
     >                       r2,r1,x)
      implicit none
      integer na
      integer indx(4)
      real*8  tol
      integer maxit
      real*8  gamma,mass(4)
      real*8  r2(3,*), r1(3,*)
      real*8  x

      logical done
      integer it
      real*8 Rab(3),Rcd(3),tRab(3),tRcd(3),mabinv,mcdinv
      real*8 Rac(3),Rbc(3),tRac(3),tRbc(3),macinv,mbcinv
      real*8 mainv,mbinv,mcinv
      real*8 A,B,C,xold,dsc,y1,y2

      if (na.eq.3) then
         Rac(1) = r1(1,indx(1)) - r1(1,indx(3))
         Rac(2) = r1(2,indx(1)) - r1(2,indx(3))
         Rac(3) = r1(3,indx(1)) - r1(3,indx(3))
         call lattice_min_difference(Rac(1),Rac(2),Rac(3))

         Rbc(1) = r1(1,indx(2)) - r1(1,indx(3))
         Rbc(2) = r1(2,indx(2)) - r1(2,indx(3))
         Rbc(3) = r1(3,indx(2)) - r1(3,indx(3))
         call lattice_min_difference(Rbc(1),Rbc(2),Rbc(3))

         tRac(1) = r2(1,indx(1)) - r2(1,indx(3))
         tRac(2) = r2(2,indx(1)) - r2(2,indx(3))
         tRac(3) = r2(3,indx(1)) - r2(3,indx(3))
         call lattice_min_difference(tRac(1),tRac(2),tRac(3))

         tRbc(1) = r2(1,indx(2)) - r2(1,indx(3))
         tRbc(2) = r2(2,indx(2)) - r2(2,indx(3))
         tRbc(3) = r2(3,indx(2)) - r2(3,indx(3))
         call lattice_min_difference(tRbc(1),tRbc(2),tRbc(3))

         macinv = 1.0d0/mass(1) + 1.0d0/mass(3)
         mbcinv = 1.0d0/mass(2) + 1.0d0/mass(3)
         mainv  = 1.0d0/mass(1)
         mbinv  = 1.0d0/mass(2)
         mcinv  = 1.0d0/mass(3)

         A = (macinv**2 - mcinv**2)*(Rac(1)**2+Rac(2)**2+Rac(3)**2)
     >     + (mcinv**2 - mbcinv**2)*(Rbc(1)**2+Rbc(2)**2+Rbc(3)**2)
     >     + 2.0d0*mcinv*(mbinv+mainv)
     >       *(Rac(1)*Rbc(1)+Rac(2)*Rbc(2)+Rac(3)*Rbc(3))

         B = 2.0d0*macinv*(tRac(1)*Rac(1)+tRac(2)*Rac(2)+tRac(3)*Rac(3))
     >     - 2.0d0*mcinv *(tRac(1)*Rbc(1)+tRac(2)*Rbc(2)+tRac(3)*Rbc(3))
     >     + 2.0d0*mbcinv*(tRbc(1)*Rbc(1)+tRbc(2)*Rbc(2)+tRbc(3)*Rbc(3))
     >     - 2.0d0*mcinv *(tRbc(1)*Rac(1)+tRbc(2)*Rac(2)+tRbc(3)*Rac(3))

         C = (tRac(1)**2 + tRac(2)**2 + tRac(3)**2)
     >     - (tRbc(1)**2 + tRbc(2)**2 + tRbc(3)**2) - gamma


         if (dabs(A).gt.1.0d-12) then
            x = -C
            it = 0
            done = .false.
            do while (.not. done)
               xold = x
               x    = -C + (1.0d0-B)*x - A*x*x
               it = it + 1
               done = (dabs(x-xold).lt.tol).or.(it.gt.maxit)
            end do

            dsc = B*B-4.0d0*A*C
            if (dsc.ge.0.0d0) then
               y1 = (-B + dsqrt(dsc))/(2.0d0*A)
               y2 = (-B - dsqrt(dsc))/(2.0d0*A)
               if (dabs(x-y1).lt.1.0d-1) then
                  x = y1
               else if (dabs(x-y2).lt.1.0d-1) then
                  x = y2
               end if
            end if
         else
            x = -C/B
            it = 0
         end if

         r2(1,indx(1)) = r2(1,indx(1)) + x*Rac(1)*mainv
         r2(2,indx(1)) = r2(2,indx(1)) + x*Rac(2)*mainv
         r2(3,indx(1)) = r2(3,indx(1)) + x*Rac(3)*mainv

         r2(1,indx(2)) = r2(1,indx(2)) - x*Rbc(1)*mbinv
         r2(2,indx(2)) = r2(2,indx(2)) - x*Rbc(2)*mbinv
         r2(3,indx(2)) = r2(3,indx(2)) - x*Rbc(3)*mbinv

         r2(1,indx(3)) = r2(1,indx(3)) + x*(Rbc(1)-Rac(1))*mcinv
         r2(2,indx(3)) = r2(2,indx(3)) + x*(Rbc(2)-Rac(2))*mcinv
         r2(3,indx(3)) = r2(3,indx(3)) + x*(Rbc(3)-Rac(3))*mcinv

c         tRac(1) = r2(1,indx(1)) - r2(1,indx(3))
c         tRac(2) = r2(2,indx(1)) - r2(2,indx(3))
c         tRac(3) = r2(3,indx(1)) - r2(3,indx(3))
c
c         tRbc(1) = r2(1,indx(2)) - r2(1,indx(3))
c         tRbc(2) = r2(2,indx(2)) - r2(2,indx(3))
c         tRbc(3) = r2(3,indx(2)) - r2(3,indx(3))
c
c         C = (tRac(1)**2 + tRac(2)**2 + tRac(3)**2)
c     >     - (tRbc(1)**2 + tRbc(2)**2 + tRbc(3)**2) - gamma
c
c         write(*,*) "x,y1,y2,constraint=",x,y1,y2,C

      else
         Rab(1) = r1(1,indx(1)) - r1(1,indx(2))
         Rab(2) = r1(2,indx(1)) - r1(2,indx(2))
         Rab(3) = r1(3,indx(1)) - r1(3,indx(2))
         call lattice_min_difference(Rab(1),Rab(2),Rab(3))

         Rcd(1) = r1(1,indx(3)) - r1(1,indx(4))
         Rcd(2) = r1(2,indx(3)) - r1(2,indx(4))
         Rcd(3) = r1(3,indx(3)) - r1(3,indx(4))
         call lattice_min_difference(Rcd(1),Rcd(2),Rcd(3))

         tRab(1) = r2(1,indx(1)) - r2(1,indx(2))
         tRab(2) = r2(2,indx(1)) - r2(2,indx(2))
         tRab(3) = r2(3,indx(1)) - r2(3,indx(2))
         call lattice_min_difference(tRab(1),tRab(2),tRab(3))

         tRcd(1) = r2(1,indx(3)) - r2(1,indx(4))
         tRcd(2) = r2(2,indx(3)) - r2(2,indx(4))
         tRcd(3) = r2(3,indx(3)) - r2(3,indx(4))
         call lattice_min_difference(tRcd(1),tRcd(2),tRcd(3))

         mabinv = 1.0d0/mass(1) + 1.0d0/mass(2)
         mcdinv = 1.0d0/mass(3) + 1.0d0/mass(4)

         A = (mabinv**2)*(Rab(1)**2+Rab(2)**2+Rab(3)**2)
     >     - (mcdinv**2)*(Rcd(1)**2+Rcd(2)**2+Rcd(3)**2)
   
         B = 2.0d0*mabinv*(tRab(1)*Rab(1)+tRab(2)*Rab(2)+tRab(3)*Rab(3))
     >     + 2.0d0*mcdinv*(tRcd(1)*Rcd(1)+tRcd(2)*Rcd(2)+tRcd(3)*Rcd(3))

         C = (tRab(1)**2 + tRab(2)**2 + tRab(3)**2)
     >     - (tRcd(1)**2 + tRcd(2)**2 + tRcd(3)**2) - gamma

         if (dabs(A).gt.1.0d-12) then
            x = -C
            it = 0
            done = .false.
            do while (.not. done)
               xold = x
               x    = -C + (1.0d0-B)*x - A*x*x
               it = it + 1
               done = (dabs(x-xold).lt.tol).or.(it.gt.maxit)
            end do

            dsc = B*B-4.0d0*A*C
            if (dsc.ge.0.0d0) then
               y1 = (-B + dsqrt(dsc))/(2.0d0*A)
               y2 = (-B - dsqrt(dsc))/(2.0d0*A)
               if (dabs(x-y1).lt.1.0d-1) then
                  x = y1
               else if (dabs(x-y2).lt.1.0d-1) then
                  x = y2
               end if
            end if
         else
            x = -C/B
            it = 0
         end if

         r2(1,indx(1)) = r2(1,indx(1)) + x*Rab(1)/mass(1)
         r2(2,indx(1)) = r2(2,indx(1)) + x*Rab(2)/mass(1)
         r2(3,indx(1)) = r2(3,indx(1)) + x*Rab(3)/mass(1)

         r2(1,indx(2)) = r2(1,indx(2)) - x*Rab(1)/mass(2)
         r2(2,indx(2)) = r2(2,indx(2)) - x*Rab(2)/mass(2)
         r2(3,indx(2)) = r2(3,indx(2)) - x*Rab(3)/mass(2)

         r2(1,indx(3)) = r2(1,indx(3)) - x*Rcd(1)/mass(3)
         r2(2,indx(3)) = r2(2,indx(3)) - x*Rcd(2)/mass(3)
         r2(3,indx(3)) = r2(3,indx(3)) - x*Rcd(3)/mass(3)

         r2(1,indx(4)) = r2(1,indx(4)) + x*Rcd(1)/mass(4)
         r2(2,indx(4)) = r2(2,indx(4)) + x*Rcd(2)/mass(4)
         r2(3,indx(4)) = r2(3,indx(4)) + x*Rcd(3)/mass(4)

      end if


      return
      end
