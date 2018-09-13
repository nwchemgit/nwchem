

*     **************************************
*     *                                    *
*     *      metadynamics_sphereoverlap    *
*     *                                    *
*     **************************************

*  This routine computes the intersection volume of two spheres
* of radii Ra and Rb seperated by a distance d.  It also calculates its
* derivatives wrt to Ra, Rb and d.
*
      subroutine metadynamics_sphereoverlap(Ra,Rb,d,s,dsdRa,dsdRb,dsdd)
      implicit none
      real*8 Ra,Rb,d
      real*8 s,dsdRa,dsdRb,dsdd

*     **** local variables ****
      real*8 pi


      if (d.le.(Ra+Rb)) then
         pi = 4.0d0*datan(1.0d0)
         s = pi*( (Ra+Rb-d)**2 
     >        *(d**2 + 2*d*(Ra+Rb) + 6.0d0*Ra*Rb - 3.0d0*(Ra**2+Rb**2)))
     >        /(12.0d0*d)

         dsdRa = pi*(Ra*(Ra+Rb-d)*(Rb-Ra+d))/d
         dsdRb = pi*(Rb*(Ra+Rb-d)*(Ra-Rb+d))/d
         dsdd  = pi*(  d**4 + Ra**4 + Rb**4 
     >               - 2.0d0*(Ra**2*Rb**2 + d**2*(Ra**2+Rb**2)))
     >              /(4.0d0*d**2)
      else
         s     = 0.0d0
         dsdRa = 0.0d0
         dsdRb = 0.0d0
         dsdd  = 0.0d0
      end if
      return
      end

*     *******************************************
*     *                                         *
*     *      metadynamics_coordsphereoverlap    *
*     *                                         *
*     *******************************************

*  This routine computes the intersection volume of two spheres centered
*  at R1 and R2 with volumes defined by V1=V0+dV0*N1 and  V2=V0+dV0*N2, where N1 and N2 
*  are the coordination numbers at each center.  It also calculates the 
*  derivatives wrt to R1, R2, N1, and N2.
*
*  Entry - V0: volume of base sphere
*        - dV0: defined so that volume of sphere of coorindation 1 is V0+dV0
*          N1,N2: coordinations numbers
*          R1,R2: atom positions
*  Exit - V: intersection volume
*         dVdN1,dVdN2: partial derivatives of V wrt N1 and N2
*         dVdR1,dVdR2: partial derivatives of V wrt R1 and R2

      subroutine metadynamics_coordsphereoverlap(V0,dV0,N1,N2,R1,R2,
     >                                           V,
     >                                           dVdN1,dVdN2,
     >                                           dVdR1,dVdR2)
      implicit none
      real*8 V0,dV0
      real*8 N1,N2
      real*8 R1(3),R2(3)
      real*8 V,dVdN1,dVdN2
      real*8 dVdR1(3),dVdR2(3)

*     **** local variables ****
      real*8 onethird
      parameter (onethird=1.0d0/3.0d0)
      real*8 dx,dy,dz,d,fourpi,Ra,Rb,dRadN1,dRbdN2,dVdRa,dVdRb,dVdd

      fourpi = 16.0d0*datan(1.0d0)
      Ra = (3.0d0*(V0+dV0*N1)/fourpi)**(onethird)
      Rb = (3.0d0*(V0+dV0*N2)/fourpi)**(onethird)
      dRadN1 = dV0/(((9.0d0*fourpi)*(V0+dV0*N1)**2)**(onethird))
      dRbdN2 = dV0/(((9.0d0*fourpi)*(V0+dV0*N2)**2)**(onethird))
      dx = R2(1) - R1(1)
      dy = R2(2) - R1(2)
      dz = R2(3) - R1(3)
      call lattice_min_difference(dx,dy,dz)
      d = dsqrt(dx**2 + dy**2 + dz**2)

      call metadynamics_sphereoverlap(Ra,Rb,d,V,dVdRa,dVdRb,dVdd)
      dVdN1 = dVdRa*dRadN1
      dVdN2 = dVdRb*dRbdN2
      dVdR1(1) = -dVdd*dx/d
      dVdR1(2) = -dVdd*dy/d
      dVdR1(3) = -dVdd*dz/d
      dVdR2(1) =  dVdd*dx/d
      dVdR2(2) =  dVdd*dy/d
      dVdR2(3) =  dVdd*dz/d

      return
      end


*     ****************************************************
*     *                                                  *
*     *           metadynamics_coordnumber               *
*     *                                                  *
*     ****************************************************

      real*8 function metadynamics_coordnumber(sprik,n,m,r0,
     >                                         ii,nw,indxw,rion)
      implicit none
      logical sprik
      real*8 n,m,r0
      integer ii,nw,indxw(*),rion(3,*)

*     **** local variables ****
      integer j
      real*8  f,dx,dy,dz,x1,y1,z1,x2,y2,z2,r

      f = 0.0d0
      x1 = rion(1,ii)
      y1 = rion(2,ii)
      z1 = rion(3,ii)
      do j=1,nw
         x2 = rion(1,indxw(j))
         y2 = rion(2,indxw(j))
         z2 = rion(3,indxw(j))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r = dsqrt(dx**2 + dy**2 + dz**2)
         if (sprik) then
            f = f + 1.d0/(1.d0+dexp(n*(r-r0)))
         else
            f = f + (1.0d0-(r/r0)**n)/(1.0d0-(r/r0)**m)
         end if
      end do

      metadynamics_coordnumber = f
      return
      end

*     ****************************************************
*     *                                                  *
*     *         metadynamics_coordnumber_force           *
*     *                                                  *
*     ****************************************************
      subroutine metadynamics_coordnumber_force(sprik,n,m,r0,
     >                                          ii,nw,indxw,dv,
     >                                          rion,fion)
      implicit none
      logical sprik
      real*8  n,m,r0
      integer ii,nw,indxw(*)
      real*8  dv
      real*8  rion(3,*),fion(3,*)

*     *** local variables ****
      integer j
      real*8  df,dx,dy,dz,x1,y1,z1,x2,y2,z2,r,rn,rm

      x1 = rion(1,ii)
      y1 = rion(2,ii)
      z1 = rion(3,ii)
      do j=1,nw
         x2 = rion(1,indxw(j))
         y2 = rion(2,indxw(j))
         z2 = rion(3,indxw(j))
         dx = x1-x2
         dy = y1-y2
         dz = z1-z2
         call lattice_min_difference(dx,dy,dz)
         r = dsqrt(dx**2 + dy**2 + dz**2)
         if (sprik) then
            rn = 1.d0+dexp(n*(r-r0))
            df = -n*(rn-1.d0)/(rn*rn)
         else
            rn = (1.0d0-(r/r0)**n)
            rm = (1.0d0-(r/r0)**m)
            df = (-n*rm/r0*(r/r0)**(n-1)+m*rn/r0*(r/r0)**(m-1))/(rm**2)
         end if
         fion(1,ii) = fion(1,ii) - (dx/r)*df*dv
         fion(2,ii) = fion(2,ii) - (dy/r)*df*dv
         fion(3,ii) = fion(3,ii) - (dz/r)*df*dv
         fion(1,indxw(j)) = fion(1,indxw(j)) + (dx/r)*df*dv
         fion(2,indxw(j)) = fion(2,indxw(j)) + (dy/r)*df*dv
         fion(3,indxw(j)) = fion(3,indxw(j)) + (dz/r)*df*dv
      end do

      return
      end


*     ****************************************************
*     *                                                  *
*     *         metadynamics_coordspherediff             *
*     *                                                  *
*     ****************************************************
      real*8 function metadynamics_coordspherediff(V0,dV0,
     >                                        sprik,n,m,r0,
     >                                        n1,indx1,
     >                                        nw,indxw,
     >                                        rion)
      implicit none
      real*8 V0,dV0
      logical sprik
      real*8  n,m,r0
      integer n1,indx1(*)
      integer nw,indxw(*)
      real*8  rion(3,*)

*     **** local variables ***
      real*8  onethird,Vs
      parameter (onethird=1.0d0/3.0d0,Vs=1.0d0)
      integer i,isgn
      real*8  f,vscal,C1,C2,V,dVdC1,dVdC2,dVdR1(3),dVdR2(3)

*     **** external functions ****
      real*8   metadynamics_coordnumber
      external metadynamics_coordnumber


      f = 0.0d0
      if (n1.gt.1) then
         vscal = (3.0d0/(16.0d0*datan(1.0d0)))**(onethird)
         isgn = 1
         C1 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                 indx1(1),nw,indxw,rion)
         C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                 indx1(2),nw,indxw,rion)

         call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(1)),
     >                                       rion(1,indx1(2)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
         f = vscal*((Vs+V)**onethird)*isgn
         isgn = -1*isgn
         if (n1.gt.2) then
            do i=2,n1-2
               C1 = C2
               C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                    indx1(i+1),nw,indxw,rion)
               call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(i)),
     >                                       rion(1,indx1(i+1)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
               f    = f + 2*vscal*((Vs+V)**onethird)*isgn
               isgn = -1*isgn
            end do
            C1 = C2
            C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                    indx1(n1),nw,indxw,rion)
            call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(n1-1)),
     >                                       rion(1,indx1(n1)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
            f = f + vscal*((Vs+V)**onethird)*isgn
         end if
      end if

      metadynamics_coordspherediff = f
      return
      end

*     ****************************************************
*     *                                                  *
*     *       metadynamics_coordspherediff_force         *
*     *                                                  *
*     ****************************************************

      subroutine metadynamics_coordspherediff_force(dv,V0,dV0,
     >                                        sprik,n,m,r0,
     >                                        n1,indx1,
     >                                        nw,indxw,
     >                                        rion,
     >                                        fion)
      implicit none
      real*8 dv
      real*8 V0,dV0
      logical sprik
      real*8  n,m,r0
      integer n1,indx1(*)
      integer nw,indxw(*)
      real*8  rion(3,*)
      real*8  fion(3,*)
 
*     **** local variables ***
      real*8  onethird,mtwothird,Vs
      parameter (onethird=1.0d0/3.0d0,mtwothird=-2.0d0/3.0d0,Vs=1.0d0)
      integer i,isgn
      real*8  f,vscal,dfidV,C1,C2,V,dVdC1,dVdC2,dVdR1(3),dVdR2(3)

*     **** external functions ****
      real*8   metadynamics_coordnumber
      external metadynamics_coordnumber

      if (n1.gt.1) then
         vscal = (3.0d0/(16.0d0*datan(1.0d0)))**(onethird)
         isgn = 1
         C1 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                 indx1(1),nw,indxw,rion)
         C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                 indx1(2),nw,indxw,rion)

         call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(1)),
     >                                       rion(1,indx1(2)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
         dfidV = vscal*onethird*(Vs+V)**mtwothird
         fion(1,indx1(1)) = fion(1,indx1(1)) - dv*dfidV*dVdR1(1)*isgn
         fion(2,indx1(1)) = fion(2,indx1(1)) - dv*dfidV*dVdR1(2)*isgn
         fion(3,indx1(1)) = fion(3,indx1(1)) - dv*dfidV*dVdR1(3)*isgn
         fion(1,indx1(2)) = fion(1,indx1(2)) - dv*dfidV*dVdR2(1)*isgn
         fion(2,indx1(2)) = fion(2,indx1(2)) - dv*dfidV*dVdR2(2)*isgn
         fion(3,indx1(2)) = fion(3,indx1(2)) - dv*dfidV*dVdR2(3)*isgn
         call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                                  indx1(1),nw,indxw,
     >                                  dv*dfidV*dVdC1*isgn,
     >                                  rion,fion)
         call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                                  indx1(2),nw,indxw,
     >                                  dv*dfidV*dVdC2*isgn,
     >                                  rion,fion)
         isgn = -1*isgn
         if (n1.gt.2) then
            do i=2,n1-2
               C1 = C2
               C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                      indx1(i+1),nw,indxw,rion)
               call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(i)),
     >                                       rion(1,indx1(i+1)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
               dfidV = vscal*onethird*(Vs+V)**mtwothird
               fion(1,indx1(n1-1)) = fion(1,indx1(n1-1))
     >                             - dv*dfidV*dVdR1(1)*2.0d0*isgn
               fion(2,indx1(n1-1)) = fion(2,indx1(n1-1))
     >                             - dv*dfidV*dVdR1(2)*2.0d0*isgn
               fion(3,indx1(n1-1)) = fion(3,indx1(n1-1))
     >                             - dv*dfidV*dVdR1(3)*2.0d0*isgn
               fion(1,indx1(n1))   = fion(1,indx1(n1))
     >                             - dv*dfidV*dVdR2(1)*2.0d0*isgn
               fion(2,indx1(n1))   = fion(2,indx1(n1))
     >                             - dv*dfidV*dVdR2(2)*2.0d0*isgn
               fion(3,indx1(n1))   = fion(3,indx1(n1))
     >                             - dv*dfidV*dVdR2(3)*2.0d0*isgn
               call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                         indx1(i),nw,indxw,
     >                         dv*dfidV*dVdC1*2.0d0*isgn,
     >                         rion,fion)
               call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                         indx1(i+1),nw,indxw,
     >                         dv*dfidV*dVdC2*2.0d0*isgn,
     >                         rion,fion)
               isgn = -1*isgn
            end do
            C1 = C2
            C2 = metadynamics_coordnumber(sprik,n,m,r0,
     >                                    indx1(n1),nw,indxw,rion)
            call metadynamics_coordsphereoverlap(V0,dV0,C1,C2,
     >                                       rion(1,indx1(n1-1)),
     >                                       rion(1,indx1(n1)),
     >                                       V,
     >                                       dVdC1,dVdC2,
     >                                       dVdR1,dVdR2)
            dfidV = vscal*onethird*(Vs+V)**mtwothird
            fion(1,indx1(n1-1)) = fion(1,indx1(n1-1))
     >                          - dv*dfidV*dVdR1(1)*isgn
            fion(2,indx1(n1-1)) = fion(2,indx1(n1-1))
     >                          - dv*dfidV*dVdR1(2)*isgn
            fion(3,indx1(n1-1)) = fion(3,indx1(n1-1))
     >                          - dv*dfidV*dVdR1(3)*isgn
            fion(1,indx1(n1))   = fion(1,indx1(n1))
     >                          - dv*dfidV*dVdR2(1)*isgn
            fion(2,indx1(n1))   = fion(2,indx1(n1))
     >                          - dv*dfidV*dVdR2(2)*isgn
            fion(3,indx1(n1))   = fion(3,indx1(n1))
     >                          - dv*dfidV*dVdR2(3)*isgn
            call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                               indx1(n1-1),nw,indxw,
     >                               dv*dfidV*dVdC1*isgn,
     >                               rion,fion)
            call metadynamics_coordnumber_force(sprik,n,m,r0,
     >                              indx1(n1),nw,indxw,
     >                              dv*dfidV*dVdC2*isgn,
     >                              rion,fion)
         end if

      end if
      return
      end
c $Id$
