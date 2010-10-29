C $Id$
      Subroutine drdy_rotx (x1,x2,scr,natom,ndim)
C
C  Rotate orientation x2 so principle moments of interia align with those
C  of orientation x1
C
      Implicit none
      Integer natom,ndim
C
      Double Precision x1(ndim),x2(ndim),scr(ndim)
C      
      Integer i,ia,ia1,ia2,ii,ii1,ii2,iv,j,k
      Logical lcol1,lcol2,ldone
      Double Precision amat(3,3),rr1(3,3),rr2(3,3),rri(3,3)
      Double Precision dot12,dotp,dotp1,dotp2,dotp1mx,dotp2mx,r2mag,sgn,
     * sum,xnorm1,xnorm2,xnorm11,xnorm21,zz
      Double Precision eps,one
      Data eps/1.d-6/,one/1.0d00/
      Integer fu6
      Data fu6/6/
C
C      write(fu6,*) ' entered rotx, call drdy_moiten'
      ldone = .false.
C Construct and diagonalize moment of inertia tensors
      call drdy_moiten(x1,amat,scr,rri,rr1,rri(1,2),rri(1,3),lcol1,
     *   natom,ndim)
      call drdy_moiten(x2,amat,scr,rri,rr2,rri(1,2),rri(1,3),lcol2,
     *   natom,ndim)
C      write(fu6,*) ' rr1(k,1)='
C      write(fu6,6000) (rr1(k,1),k=1,3)
C      write(fu6,*) ' rr2(k,1)='
C      write(fu6,6000) (rr2(k,1),k=1,3)
      if (lcol1.or.lcol2) then
C Either one or both of the orientations are collinear, just need to align
C  the two collinear vectors, or align the collinear vector and the 
C  priniciple axis with the smallest eigenvalue
C  check if rr1(k,1) and rr2(k,1) are already collinear (dot product is -1 or 1)
         dot12 = 0.0d00
C         write(fu6,*) ' rr1,rr2,dot12='
         do k = 1,3
            dot12 = dot12 + rr1(k,1)*rr2(k,1)
C            write(fu6,6000) rr1(k,1),rr2(k,1),dot12
         enddo
         ldone = 1.0d00-abs(dot12).lt.eps
C         write(fu6,*) ' collinear orientation, dot12,ldone='
C         write(fu6,*) dot12,ldone
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
      else
C Set signs of eigenvectors so that dot product of eigenvector with vector to 
C  one of the atoms has the same sign for both orientations.  Use atom that
C  give the largest magnitude dot products.
         do iv = 1,3
C            write (fu6,*) ' get sign for eignvector',iv
            ii = 0
            dotp1mx = 0.0d00
            dotp2mx = 0.0d00
            sgn = 0.0
C            write (fu6,*) ' ia,dotp1,dotp2='
            do ia = 1,natom
               dotp1 = 0.0d00
               dotp2 = 0.0d00
               xnorm1 = 0.0d00
               xnorm2 = 0.0d00
               do k = 1,3
                  dotp1 = dotp1 + rr1(k,iv)*x1(ii+k)
                  dotp2 = dotp2 + rr2(k,iv)*x2(ii+k)
                  xnorm1 = xnorm1 + x1(ii+k)*x1(ii+k)
                  xnorm2 = xnorm2 + x2(ii+k)*x2(ii+k)
               enddo
               if (xnorm1.gt.0.0d00) dotp1 = dotp1/sqrt(xnorm1)
               if (xnorm2.gt.0.0d00) dotp2 = dotp2/sqrt(xnorm2)
               if (abs(dotp1).gt.dotp1mx.and.abs(dotp2).gt.dotp2mx) then
                  dotp1mx = abs(dotp1)
                  dotp2mx = abs(dotp2)
                  sgn = sign(one,dotp1*dotp2)
               endif
C               write (fu6,*) ia,dotp1,dotp2
               ii = ii + 3
            enddo
            if (dotp1mx.lt.eps .or. dotp2mx.lt.eps) then
C  Planar arrangement, get cross product of vectors for two atoms 
C   and compare with eigenvector
C               write (fu6,*) ' Planar'
               ia = 0
               xnorm1 = 0.0d00
               xnorm2 = 0.0d00
C  Get first atom (not at origin) for constructing cross product
               do while (ia.lt.natom .and. 
     *                  (xnorm1.le.eps .or. xnorm2.le.eps))
                  ia = ia + 1
                  ii = 3*(ia-1)
                  xnorm1 = 0.0d00
                  xnorm2 = 0.0d00
                  do k = 1,3
                     xnorm1 = xnorm1 + x1(ii+k)*x1(ii+k)
                     xnorm2 = xnorm2 + x2(ii+k)*x2(ii+k)
                  enddo
                  xnorm1 = sqrt(xnorm1) 
                  xnorm2 = sqrt(xnorm2)
C                  write (fu6,*) ' ia,xnorm1,xnorm2='
C                  write (fu6,*) ia,xnorm1,xnorm2
               enddo
               ia1 = ia
               ia2 = ia
               ii1 = ii
               xnorm11 = xnorm1
               xnorm21 = xnorm2
               xnorm1 = 0.0d00
               xnorm2 = 0.0d00
C  Get second atom (not at origin and not collinear with first)
               dotp1mx = 0.0d00
               dotp2mx = 0.0d00
               do while (ia.lt.natom)
                  ia = ia + 1
                  ii = 3*(ia-1)
                  xnorm1 = 0.0d00
                  xnorm2 = 0.0d00
                  do k = 1,3
                     xnorm1 = xnorm1 + x1(ii+k)*x1(ii+k)
                     xnorm2 = xnorm2 + x2(ii+k)*x2(ii+k)
                  enddo
                  xnorm1 = sqrt(xnorm1) 
                  xnorm2 = sqrt(xnorm2)
C                  write (fu6,*) ' ia,xnorm1,xnorm2='
C                  write (fu6,*) ia,xnorm1,xnorm2
                  if(xnorm1.gt.eps .and. xnorm2.gt.eps) then
                     dotp1 = 0.0d00
                     dotp2 = 0.0d00
                     do k = 1,3
                        dotp1 = dotp1 + x1(ii1+k)*x1(ii+k)
                        dotp2 = dotp2 + x2(ii1+k)*x2(ii+k)
                     enddo
                     dotp1 = dotp1/(xnorm1*xnorm11)
                     dotp2 = dotp2/(xnorm2*xnorm21)
C                     write (fu6,*) ' ia,dotp1,dotp2='
C                     write (fu6,*) ia,dotp1,dotp2
                     if (1.0d00-abs(dotp1).gt.dotp1mx .and. 
     *                1.0d00-abs(dotp2).gt.dotp2mx) then
                        dotp1mx = abs(dotp1)
                        dotp2mx = abs(dotp2)
                        ia2 = ia
                     endif
                  endif
               enddo
               if (ia2.le.ia1) then
                  write (fu6,*) ' problem in rotx setting sign on',
     *             ' eigenvector',iv
                  sgn = 1.0d00
               else
C form the cross product of vectors for atoms ia1 and ia2 and get
C  dot product of this vector with eigenvector
                  ii1 = 3*(ia1-1)
                  ii2 = 3*(ia2-1)
                  scr(1) = x1(ii1+2)*x1(ii2+3) - x1(ii1+3)*x1(ii2+2)
                  scr(2) = x1(ii1+3)*x1(ii2+1) - x1(ii1+1)*x1(ii2+3)
                  scr(3) = x1(ii1+1)*x1(ii2+2) - x1(ii1+2)*x1(ii2+1)
                  dotp1 = 0.0d00
                  do k = 1,3
                     dotp1 = dotp1 + scr(k)*rr1(k,iv)
                  enddo
                  scr(1) = x2(ii1+2)*x2(ii2+3) - x2(ii1+3)*x2(ii2+2)
                  scr(2) = x2(ii1+3)*x2(ii2+1) - x2(ii1+1)*x2(ii2+3)
                  scr(3) = x2(ii1+1)*x2(ii2+2) - x2(ii1+2)*x2(ii2+1)
                  dotp2 = 0.0d00
                  do k = 1,3
                     dotp2 = dotp2 + scr(k)*rr2(k,iv)
                  enddo
                  sgn = sign(one,dotp1*dotp2)
               endif
            endif
C            write (fu6,*) ' sgn=', sgn
            if (sgn.lt.0.0d00) then
               do k = 1,3
                  rr2(k,iv) = -rr2(k,iv)
               enddo
            endif
         enddo
      endif
      if(.not.ldone) then
C  Invert rr2 matrix
C    Get det(rr)
C         write(fu6,*) ' rr1='
C         write(fu6,*) ((rr1(i,j),j=1,3),i=1,3)
C         write(fu6,*) ' rr2='
C         write(fu6,*) ((rr2(i,j),j=1,3),i=1,3)
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
C         write(fu6,*) ' rri='
C         write(fu6,*) ((rri(i,j),j=1,3),i=1,3)
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
C      write(fu6,*) ' amat='
C      write(fu6,*) ((amat(j,i),i=1,3),j=1,3)
C
C Rotate orientation x2
      do i = 1,natom
         ii = 3*(i-1)
         do j = 1,3
            sum=0.0
            do k = 1,3
               sum = sum + amat(k,j)*x2(ii+k)
            enddo
            scr(j) = sum
         enddo
         do j = 1,3
            x2(ii+j) = scr(j)
         enddo
      enddo
C Calculate dot product of x1 and x2
      dotp = 0.0d00
C      write(fu6,*) ' x1,x2,dotp='
      do i = 1,3*natom
         dotp = dotp + x2(i)*x1(i)
C         write(fu6,*) x1(i),x2(i),dotp
      enddo
C      write(fu6,*) ' dotp=',dotp
      if (dotp.lt.0.0d00) then
         do i = 1,3*natom
            x2(i) = -x2(i)
         enddo
      endif
      return
6000  Format(1x,1p3e13.5)
      end
