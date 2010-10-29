*
* $Id$
*

      subroutine paw_dirac_x(n2ft3d,ng,ispin,dn,xcp,xce,x)
      implicit none

      integer ng
      integer n2ft3d
      integer ispin
      real*8 dn(n2ft3d,2)
      real*8 xcp(n2ft3d,2)
      real*8 xce(n2ft3d,2)
      real*8 x(n2ft3d)

      integer k
      double precision coeff1,coeff2,n_tmp
      parameter (coeff1=0.984745021842697d0)
      parameter (coeff2=0.738558766382022d0)

      do k=1,ng
         x(k)=(dn(k,1)+dn(k,2)+1.0d-30)**(1.0d0/3)
         xcp(k,1) = -coeff1*x(k)
         xcp(k,2) = xcp(k,1)
         xce(k,1) = -coeff2*x(k)
      end do  

      return
      end
