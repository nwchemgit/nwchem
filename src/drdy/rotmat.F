C $Id$
      Subroutine drdy_rotmat (x,xs,amat,rr,rrs,rri,natom,ndim)
C
C  Calculate transformation matrix, amat, from rotated orientation,
C     vector xs, to original orientation, vector x.  
C
      Implicit none
      Logical done
      Integer natom,ndim
      Double Precision x(ndim),xs(ndim),amat(3,3),
     * rr(3,3),rrs(3,3),rri(3,3)
C
      Integer i,iat,j,k,koff
      Double Precision dotp,rr1mag,rr2mag,sum,zz
      Double Precision eps,one
      Data eps/1.d-4/, one/1.0d00/
C
C  form two orthogonal vectors for each configuration (Gaussian 
C     standard, xs, and original input, x).
C     1 - vector from atom 1 to 2
C     2 - vector from midpoint between atoms 1 and 2 to third atom
C         the third atom should not be collinear with 1 and 2
C         if all atoms are collinar, then use special procedure
      rr1mag = 0.0d00
C First vector is between atoms 1 and 2
      do k = 1,3
         rr(1,k) = x(k+3)-x(k)
         rrs(1,k) = xs(k+3)-xs(k)
         rr1mag = rr1mag + rr(1,k)*rr(1,k)
      enddo
      rr1mag = sqrt(rr1mag)
C Loop over other atoms until a second vector is formed that is not
C  collinear with vector 1.  (The dot product is -1 or 1 for collinear.)
      dotp=1.0d00
      iat = 2
      do while (1.0-abs(dotp).lt.eps .and. iat.lt.natom)
         iat = iat + 1
         koff = 3*(iat-1)
         dotp = 0.0d00
         rr2mag = 0.0d00
         do k = 1,3
            rr(2,k) = 0.5d00*(x(k)+x(k+3))-x(k+koff)
            rrs(2,k) = 0.5d00*(xs(k)+xs(k+3))-xs(k+koff)
            dotp = dotp + rr(1,k)*rr(2,k)
            rr2mag = rr2mag + rr(2,k)*rr(2,k)
         enddo
         rr2mag = sqrt(rr2mag)
         dotp = dotp/(rr1mag*rr2mag)
      enddo
C
      done = .false.
      if (1.0d00-abs(dotp).lt.eps) then
C all atoms are collinear, form a second vector that is the cross
C  product of rr(1,k) and rrs(1,k)
         rr(2,1) = rr(1,2)*rrs(1,3) - rr(1,3)*rrs(1,2)
         rr(2,2) = rr(1,3)*rrs(1,1) - rr(1,1)*rrs(1,3)
         rr(2,3) = rr(1,1)*rrs(1,2) - rr(1,2)*rrs(1,1)
         rrs(2,1) = rr(2,1)
         rrs(2,2) = rr(2,2)
         rrs(2,3) = rr(2,3)
C if the magnitude of the cross product is sufficiently small,
C  then no rotation is needed (we're done)
         zz = rr(2,1)*rr(2,1)+rr(2,2)*rr(2,2)+rr(2,3)*rr(2,3)
         done = abs(zz).lt.eps
      endif
      if (done) then
C rr(1,k) and rrs(1,k) are collinear, get sign of their dot product
         dotp = 0.0d00
         do k = 1,3
            dotp = dotp + rr(1,k)*rrs(1,k)
         enddo
         dotp = sign(one,dotp)
         do i = 1,3
            do j = 1,3
               amat(i,j) = 0.0d00
            enddo
            amat(i,i) = dotp
         enddo
      else
C  form third orthogonal vector for each configuration (Gaussian 
C     standard and original input) from cross product
         rr(3,1) = rr(1,2)*rr(2,3) - rr(1,3)*rr(2,2)
         rr(3,2) = rr(1,3)*rr(2,1) - rr(1,1)*rr(2,3)
         rr(3,3) = rr(1,1)*rr(2,2) - rr(1,2)*rr(2,1)
         rrs(3,1) = rrs(1,2)*rrs(2,3) - rrs(1,3)*rrs(2,2)
         rrs(3,2) = rrs(1,3)*rrs(2,1) - rrs(1,1)*rrs(2,3)
         rrs(3,3) = rrs(1,1)*rrs(2,2) - rrs(1,2)*rrs(2,1)
C  Invert rr matrix
C    Get det(rr)
         zz =      rr(1,1)*(rr(2,2)*rr(3,3)-rr(2,3)*rr(3,2))
         zz = zz - rr(1,2)*(rr(2,1)*rr(3,3)-rr(2,3)*rr(3,1))
         zz = zz + rr(1,3)*(rr(2,1)*rr(3,2)-rr(2,2)*rr(3,1))
C
         rri(1,1) =  (rr(2,2)*rr(3,3)-rr(2,3)*rr(3,2))/zz
         rri(2,1) = -(rr(2,1)*rr(3,3)-rr(2,3)*rr(3,1))/zz
         rri(3,1) =  (rr(2,1)*rr(3,2)-rr(2,2)*rr(3,1))/zz
         rri(1,2) = -(rr(1,2)*rr(3,3)-rr(1,3)*rr(3,2))/zz
         rri(2,2) =  (rr(1,1)*rr(3,3)-rr(1,3)*rr(3,1))/zz
         rri(3,2) = -(rr(1,1)*rr(3,2)-rr(1,2)*rr(3,1))/zz
         rri(1,3) =  (rr(1,2)*rr(2,3)-rr(1,3)*rr(2,2))/zz
         rri(2,3) = -(rr(1,1)*rr(2,3)-rr(1,3)*rr(2,1))/zz
         rri(3,3) =  (rr(1,1)*rr(2,2)-rr(1,2)*rr(2,1))/zz
C
         do i = 1,3
            do j = 1,3
               sum = 0.0d00
               do k = 1,3
                  sum = sum + rri(i,k)*rrs(k,j)
               enddo
               amat(i,j) = sum
            enddo
         enddo
      endif
      return
      end
