      subroutine qhop_newpc(xr,xacc,xhydr,dist,alpha,xnew)
c
c $Id$
c     The new coordinates of the proton, after its eventual transfer,
c     are determined. 
c     
c     description of arguments
c     
c     in: xr(3)    coord. of the atom R to which the acceptor atom is bounded
c         xacc(3)  coord. of the acceptor atom A
c         xhydr(3) coord of the transferrable proton H
c         dist     equilibrium distance xacc-proton
c         alpha    R-A-H angle (in radians)
c     
c     out: xnew(3)  coord. of thr proton after transfer
c      
c     In order to avoid energy drifts, the new position
c     is given by two condition: the distance from the acceptor 
c     atom dist is constrained to be the equilibrium distance for the bond, 
c     and the angle R-A-H  is set to its equilibrium value too.
c     The points satisfying this condition lie on a circle
c     of radius dist*cos(pi-alpha). The point on this cirle minimizing
c     the distance with the actual H position is chosen to be the new 
c     position of H after the eventual transfer.
c     The problem can be analitically solved in a new coordinate system,
c     with the z-axis parallel to the R-A axis. 
c     
c     The coordinate of H are first trandformed, via an affine transformation
c     in the new referement, then the new H position is determined in this
c     new system and finally it is retransformed in the actual referement 
c     system.
c     
c     The matrix notation of the rotation is:
c     
c            cos(theta)*cos(phi)    -sin(phi)           cos(phi)*sin(theta)
c           -sin(phi)*cos(theta)     cos(phi)          -sin(phi)*sin(theta)    
c               -sin(theta)              0                  cos(theta)
c
      implicit none
c
c     declaration of arguments
c     
      real*8 xr(3),xacc(3),xhydr(3),dist,alpha,xnew(3)
c
c     declaration of arguments
c     
      real*8 xr(3),xacc(3),xhydr(3),dist,alpha,xnew(3)
c
c     declaration of local variables
c
      real*8 costh,sinth,cosphi,sinphi,theta,phi
      real*8 trasl,radius,factor
      real*8 xhydr2(3),xnew2(2),xdr(3),xdr1prime,xdr2prime,norm
c     
      parameter(tol=1E-10)
      parameter(pi=3.14159)
c      

      do i=1,3
         xdr(i)=xacc(i)-xr(i)
      enddo          

      norm=sqrt(xdr(1)**2+xdr(2)**2+xdr(3)**2)
      trasl=-dist*cos(alpha)+norm
      radius=dist*sin(pi-alpha)
  
      if(xdr(1).lt.tol.and.xdr(2).lt.tol)then
c   
c     no rotation needed
c            
         factor=sqrt(xhydr2(1)**2+xhydr2(2)**2)
         if(factor.gt.tol)then
            xnew(1)=xhydr(1)*(radius/factor)
            xnew(2)=xhydr(2)*(radius/factor)
         else
c
c     H projection on the circle plane coincides with the center of 
c     the circle, all points are equally distant
c     
         xnew(1)=radius
         xnew(2)=0.
         endif
      return
      else
         phi=atan(-(xdr(2)/(xdr(1))))
      endif

      sinphi=sin(phi)
      cosphi=cos(phi)

      xdr1prime=xdr(1)*cos(phi)-xdr(2)*sin(phi)
      if((xdr(3)).lt.tol)then
         theta=pi/2.
      else
         theta=atan(xdr1prime/xdr(3))
      endif
      xdr2prime=sin(theta)*xdr1prime+cos(theta)*xdr(3)
      if(xdr2prime.lt.0)theta=theta+pi
      
      sinth=sin(theta)
      costh=cos(theta)

c     
c     trasforming the coordinates of the hydrogen in the new system
c     
      xhydr2(1)=cosphi*costh*(xhydr(1)-xr(1))-costh*sinphi*(xhydr(2)
     $     -xr(2))-sinth*(xhydr(3)-xr(3))
      xhydr2(2)=sinphi*(xhydr(1)-xr(1))+cosphi*(xhydr(2)-xr(2))
      xhydr2(3)=sinth*cosphi*(xhydr(1)-xr(1))-sinth*sinphi*(xhydr(2)
     $     -xr(2))+costh*(xhydr(3)-xr(3))
         
c
c     calculating the new hydrogen coordinates in transformed system
c      
      factor=sqrt(xhydr2(1)**2+xhydr2(2)**2)
      if(factor.gt.tol)then
         xnew2(1)=xhydr2(1)*(radius/factor)
         xnew2(2)=xhydr2(2)*(radius/factor)
      else
         xnew2(1)=radius
         xnew2(2)=0
      endif
c     
c     Retrasforming the coordinates into the original system 
c     
      xnew(1)=costh*cosphi*xnew2(1)+sinphi*xnew2(2)+sinth*cosphi
     $     *trasl+xr(1)
      xnew(2)=-costh*sinphi*xnew2(1)+cosphi*xnew2(2)-sinth*sinphi
     $     *trasl+xr(2)
      xnew(3)=-sinth*xnew2(1)+costh*trasl+xr(3)
      return
      end
