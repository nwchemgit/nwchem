      subroutine hnd_diag(vec,eig,ndim,order,righthand)
c
c $Id: hnd_diag.F,v 1.2 2004-07-22 19:37:23 bert Exp $
c
      implicit none
c
      double precision vec(ndim,ndim), eig(ndim)
      double precision test, t, zero, one, conv
      integer ndim, i
      logical order, righthand
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data conv   /1.0d-10/
c
c     ----- diagonalize the matrix
c
      if (order) then
c     ----- put the principal components of the tensor and -----
c           the corresponding eigenvectors in increasing
c           order according to the absolute value of the
c           eigenvalue ( the principal component ) .
         call util_jacobi(ndim,vec,ndim,eig)
      else
         call util_jacobinosort(ndim,vec,ndim,eig)
      endif
c
      if (righthand) then
c
c     ----- check for right handedness, correct if not -----
c
         test = vec(1,3)*( vec(2,1)*vec(3,2) - vec(3,1)*vec(2,2) )
     1        + vec(2,3)*( vec(3,1)*vec(1,2) - vec(1,1)*vec(3,2) )
     2        + vec(3,3)*( vec(1,1)*vec(2,2) - vec(2,1)*vec(1,2) )
         if(test.gt.zero) return
         if( abs(eig(1)-eig(2)).gt.conv) go to 160
         t = eig(1)
         eig(1) = eig(2)
         eig(2) = t
         do 150 i=1,3
            t = vec(i,1)
            vec(i,1) = vec(i,2)
            vec(i,2) = t
  150    continue
         return
  160    if( abs(eig(2)-eig(3)).gt.conv) go to 180
         t = eig(2)
         eig(2) = eig(3)
         eig(3) = t
         do 170 i=1,3
            t = vec(i,2)
            vec(i,2) = vec(i,3)
            vec(i,3) = t
  170    continue
         return
  180    do 190 i=1,3
            vec(i,3) = - vec(i,3)
  190    continue
      endif
c
      return
      end
