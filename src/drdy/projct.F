C $Id$
      Subroutine drdy_projct(ierr)
C
C Calculates projected force constant matrix, 
C  this routine assumes that coordinates, derivatives, and hessian 
C  are mass-weighted and that the center of mass of the coordinates
C  is at the origin
C
      Implicit none
      Integer ierr
C
#include "drdyP.fh"
C
      Integer i,ia,ib,ic,ii,indx,ip,iscr(3),j,ja,jb,jc,jend,jj,jndx,jp,
     * jrnk,k,kp,l
      Double Precision det,dmag,pgrad,prot,ptrans,rot(3,3),sum,zmom
C
      Double Precision eps,tens(3,3,3)
      Data eps/1.0d-9/
      Data tens / 5*0.0d00,-1.0d00,0.0d00,1.0d00,3*0.0d00,1.0d00,
     *            3*0.0d00,-1.0d00,3*0.0d00,-1.0d00,0.0d00, 1.0d00,
     *            5*0.0d00 /
C
C Normalize mass-weighted gradient and put it into scr1
      dmag=0.d0
      ierr = 0
      do i = 1,n3
         dmag = dmag + dx(i)*dx(i)
      enddo
      dmag = sqrt(dmag)
      do i = 1,n3
         scr1(i) = dx(i)/dmag
      enddo
      if (icode(5).eq.3) then
C Compute P matrix for linear case
C  Assume that the molecule lies along the z-axis
         zmom = 0.0d00
         do i = 1, natom
            l = 3*i
            zmom = zmom+x(l)**2
         enddo
         do i = 1,n3
            do j = 1,n3
               proj(i,j) = 0.0d00
            enddo
         enddo
         do ip = 1, natom
            indx = 3*(ip-1)
            do jp = 1, ip
               jndx = 3*(jp-1)
C  translational projector
               ptrans = amass(indx+1)*amass(jndx+1)
C  rotational projector
               prot = x(indx+3)*x(jndx+3)/zmom
C  gradient projector
               pgrad = scr1(indx+3)*scr1(jndx+3)
C  x,x component
               proj(indx+1,jndx+1) = prot + ptrans
C  y,y component
               proj(indx+2,jndx+2) = prot + ptrans
C  z,z component
               proj(indx+3,jndx+3) = pgrad + ptrans
            enddo
         enddo
      else
C Compute P matrix for nonlinear case
C  First, compute inertia tensor and invert (non-linear case)
         do i = 1, 3
            do j = 1, 3
               rot(i,j) = 0.0d00
            enddo
         enddo
         do i = 1, natom
            l = 3*(i-1)+1
            rot(1,1) = rot(1,1)+x(l+1)**2+x(l+2)**2
            rot(1,2) = rot(1,2)-x(l)*x(l+1)
            rot(1,3) = rot(1,3)-x(l)*x(l+2)
            rot(2,2) = rot(2,2)+x(l)**2+x(l+2)**2
            rot(2,3) = rot(2,3)-x(l+1)*x(l+2)
            rot(3,3) = rot(3,3)+x(l)**2+x(l+1)**2
         enddo
         rot(2,1) = rot(1,2)
         rot(3,1) = rot(1,3)
         rot(3,2) = rot(2,3)
         call drdy_mxlneq (rot,3,3,det,jrnk,eps,iscr,0)             
         if (jrnk.lt.3) then
            write (fu6,*) 
     *       ' problem with drdy_mxlneq in subroutine drdy_project'
            ierr = 1
            return
         endif
C
C Compute p matrix
         do ip = 1, natom
            indx = 3*(ip-1)
            do jp = 1, ip
               jndx = 3*(jp-1)
               do ic = 1, 3
                  ii = indx+ic
                  jend = 3
                  if (jp.eq.ip) jend = ic
                  do jc = 1, jend
                     jj = jndx+jc
                     sum = 0.0d00
C  rotational projector
                     do ia = 1,3
                        do ib = 1,3
                           do ja = 1,3
                              do jb = 1,3
                                 sum = sum+tens(ia,ib,ic)*tens(ja,jb,jc)
     *                            *rot(ia,ja)*x(indx+ib)*x(jndx+jb)
                              enddo
                           enddo
                        enddo
                     enddo
C  gradient projector
                     proj(ii,jj) = sum + scr1(ii)*scr1(jj)
                  enddo
                  jj = jndx + ic
C  translational projector
                  proj(ii,jj) = proj(ii,jj) + amass(ii)*amass(jj)
               enddo
            enddo
         enddo
      endif
C
      do i = 1, n3
         do j = 1, i
            proj(i,j) = -proj(i,j)
            proj(j,i) = proj(i,j)
         enddo
         proj(i,i) = 1.0d00+proj(i,i)
      enddo
C
C post and premultiply f by p
      indx = 0
      do i = 1, n3
         do j = 1, i
            sum = 0.0d00
            do k = 1, n3
               do kp = 1,n3
                  sum = sum + proj(i,k)*f(k,kp)*proj(kp,j)
               enddo
            enddo 
            indx = indx + 1
            hess(indx) = sum
         enddo
      enddo
      return
      end
