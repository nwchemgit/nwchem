C $Id: rotx.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_rotx (x1,x2,scr,natom,ndim)
C
C  Rotate orientation x2 so principle moments of interia align with those
C  of orientation x1
C
      Implicit None ! Double Precision (A-H,O-Z)
      Integer ndim
      Integer natom
      Double Precision x1(ndim),x2(ndim),scr(ndim)
      Integer ndima
      Parameter (ndima=3)
      Logical lcol1,lcol2,ldone
      Double Precision
     &      amat(ndima,ndima),rr1(ndima,ndima),rr2(ndima,ndima),
     &      rri(ndima,ndima)
      Integer i, j, k, ii
      Double Precision dot12, r2mag, zz, sum, dotp
      Double Precision eps
      Data eps/1.e-6/
C
C      write (6,*) ' entered rotx, call moiten'
      ldone = .false.
C Construct and diagonalize moment of inertia tensors
      call drdy_moiten(x1,amat,scr,rri,rr1,rri(1,2),rri(1,3),lcol1,
     *      natom,ndim,ndima)
      call drdy_moiten(x2,amat,scr,rri,rr2,rri(1,2),rri(1,3),lcol2,
     *      natom,ndim,ndima)
C      write (6,*) ' rr1(k,1)='
C      write (6,*) (rr1(k,1),k=1,3)
C      write (6,*) ' rr2(k,1)='
C      write (6,*) (rr2(k,1),k=1,3)
      if (lcol1.or.lcol2) then
C Either one or both of the orientations are collinear, just need to align
C  the two collinear vectors, or align the collinear vector and the 
C  priniciple axis with the smallest eigenvalue
C  check if rr1(k,1) and rr2(k,1) are already collinear (dot product is -1 or 1)
        dot12 = 0.0d00
C         write (6,*) ' rr1,rr2,dot12='
        do k = 1,3
          dot12 = dot12 + rr1(k,1)*rr2(k,1)
C            write (6,*) rr1(k,1),rr2(k,1),dot12
        enddo
        ldone = 1.0d00-abs(dot12).lt.eps
C         write (6,*) ' collinear orientation, dot12,ldone='
C         write (6,*) dot12,ldone
        if (ldone) then
          do i = 1,3
            do j = 1,3
              amat(i,j) = 0.0d00
            enddo
            amat(i,i) = dot12
          enddo
        else
C form the cross product of rr1(k,1) and rr2(k,1)
          rr1(1,2) = rr1(2,1)*rr2(3,1) - rr1(3,1)*rr2(2,1)
          rr1(2,2) = rr1(3,1)*rr2(1,1) - rr1(1,1)*rr2(3,1)
          rr1(3,2) = rr1(1,1)*rr2(2,1) - rr1(2,1)*rr2(1,1)
          r2mag = 0.0d00
          do k = 1,3
            r2mag = r2mag + rr1(k,2)*rr1(k,2)
          enddo
          r2mag = sqrt(r2mag)
          do k = 1,3
            rr1(k,2) = rr1(k,2)/r2mag
            rr2(k,2) = rr1(k,2)
          enddo
C  form third orthogonal vector for each orientation from cross product
C    of the other two vectors for that orientation
          rr1(1,3) = rr1(2,1)*rr1(3,2) - rr1(3,1)*rr1(2,2)
          rr1(2,3) = rr1(3,1)*rr1(1,2) - rr1(1,1)*rr1(3,2)
          rr1(3,3) = rr1(1,1)*rr1(2,2) - rr1(2,1)*rr1(1,2)
          rr2(1,3) = rr2(2,1)*rr2(3,2) - rr2(3,1)*rr2(2,2)
          rr2(2,3) = rr2(3,1)*rr2(1,2) - rr2(1,1)*rr2(3,2)
          rr2(3,3) = rr2(1,1)*rr2(2,2) - rr2(2,1)*rr2(1,2)
        endif
      endif
      if(.not.ldone) then
C  Invert rr2 matrix
C    Get det(rr)
C         write (6,*) ' rr1='
C         write (6,*) ((rr1(i,j),j=1,3),i=1,3)
C         write (6,*) ' rr2='
C         write (6,*) ((rr2(i,j),j=1,3),i=1,3)
        zz =      rr1(1,1)*(rr1(2,2)*rr1(3,3)-rr1(3,2)*rr1(2,3))
        zz = zz - rr1(2,1)*(rr1(1,2)*rr1(3,3)-rr1(3,2)*rr1(1,3))
        zz = zz + rr1(3,1)*(rr1(1,2)*rr1(2,3)-rr1(2,2)*rr1(1,3))
C
        rri(1,1) =  (rr1(2,2)*rr1(3,3)-rr1(3,2)*rr1(2,3))/zz
        rri(1,2) = -(rr1(1,2)*rr1(3,3)-rr1(3,2)*rr1(1,3))/zz
        rri(1,3) =  (rr1(1,2)*rr1(2,3)-rr1(2,2)*rr1(1,3))/zz
        rri(2,1) = -(rr1(2,1)*rr1(3,3)-rr1(3,1)*rr1(2,3))/zz
        rri(2,2) =  (rr1(1,1)*rr1(3,3)-rr1(3,1)*rr1(1,3))/zz
        rri(2,3) = -(rr1(1,1)*rr1(2,3)-rr1(2,1)*rr1(1,3))/zz
        rri(3,1) =  (rr1(2,1)*rr1(3,2)-rr1(3,1)*rr1(2,2))/zz
        rri(3,2) = -(rr1(1,1)*rr1(3,2)-rr1(3,1)*rr1(1,2))/zz
        rri(3,3) =  (rr1(1,1)*rr1(2,2)-rr1(2,1)*rr1(1,2))/zz
C         write (6,*) ' rri='
C         write (6,*) ((rri(i,j),j=1,3),i=1,3)
C
        do i = 1,3
          do j = 1,3
            sum = 0.0d00
            do k = 1,3
              sum = sum + rr2(i,k)*rri(k,j)
            enddo
            amat(i,j) = sum
          enddo
        enddo
      endif
C      write (6,*) ' amat='
C      write (6,*) ((amat(j,i),i=1,3),j=1,3)
C
C Rotate orientation x2
      do i = 1,natom
        ii = ndima*(i-1)
        do j = 1,ndima
          sum=0.0d00
          do k = 1,ndima
            sum = sum + amat(k,j)*x2(ii+k)
          enddo
          scr(j) = sum
        enddo
        do j = 1,ndima
          x2(ii+j) = scr(j)
        enddo
      enddo
C Calculate dot product of x1 and x2
      dotp = 0.0d00
C      write (6,*) ' x1,x2,dotp='
      do i = 1,3*natom
        dotp = dotp + x2(i)*x1(i)
C         write (6,*) x1(i),x2(i),dotp
      enddo
C      write (6,*) ' dotp=',dotp
      if (dotp.lt.0.0) then
        do i = 1,ndim
          x2(i) = -x2(i)
        enddo
      endif
      return
      end

