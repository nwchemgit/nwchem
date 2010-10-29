      subroutine hnd_diag(vec,eig,ndim,order,righthand)
c
c $Id$
c
      implicit none
c
      integer ndim, i, j, jj
      double precision vec(ndim,ndim), eig(ndim)
      double precision test, t, zero, one, conv, xx
      logical order, righthand
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data conv   /1.0d-10/
c
c     ----- diagonalize the matrix
c
      call util_jacobinosort(ndim,vec,ndim,eig)
c
      if (righthand) then
c
c     ----- check for right handedness, correct if not -----
c
         test = vec(1,3)*( vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2) )
     1        + vec(2,3)*( vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2) )
     2        + vec(3,3)*( vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2) )
         if (test.le.zero) then
            if (abs(eig(1)-eig(2)).gt.conv) then
               if (abs(eig(2)-eig(3)).gt.conv) then
                 do i=1,3
                    vec(i,3) = - vec(i,3)
                 enddo
               else
                 t = eig(2)
                 eig(2) = eig(3)
                 eig(3) = t
                 do i=1,3
                    t = vec(i,2)
                    vec(i,2) = vec(i,3)
                    vec(i,3) = t
                 enddo
               endif
            else
               t = eig(1)
               eig(1) = eig(2)
               eig(2) = t
               do i=1,3
                  t = vec(i,1)
                  vec(i,1) = vec(i,2)
                  vec(i,2) = t
               enddo
            endif
         endif
      endif
c
      if (order) then
c
c     ----- put the principal components of the tensor and -----
c           the corresponding eigenvectors in increasing
c           order according to the absolute value of the
c           eigenvalue ( the principal component ) .
c
         do 50 i=1,3
            jj=i
            do 30 j=i,3
               if( abs(eig(j)).gt. abs(eig(jj))) jj=j
   30       continue
            if(jj.eq.i) go to 50
c
            xx=eig(jj)
            eig(jj)=eig(i)
            eig(i)=xx
            do 40 j=1,3
               xx=vec(j,jj)
               vec(j,jj)=vec(j,i)
               vec(j,i)=xx
   40       continue
c
   50    continue
      endif
c
      return
      end
