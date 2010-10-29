*
* $Id$
*


      subroutine Check_Real(nx,ny,nz,nfft3d,A)
      implicit none
      integer nx,ny,nz,nfft3d
      complex*16 A(nfft3d)


*     **** local variables ****
      integer i,j,k,kr
      integer k1,k2,k3
      integer nxh,nyh,nzh
      integer indx,indx_r

      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

*     **** k1=0,k2=0, k3=-(nzh-1),nzh-1 ****
      k1 = 0
      k2 = 0
      do k3=1,(nzh-1)
        indx   = k3     *(nxh+1)*ny + 1
        indx_r = (nz-k3)*(nxh+1)*ny + 1
        if (A(indx).ne.dconjg(A(indx_r))) then
            write(*,*) "NR:",k1,k2,k3,A(indx),A(indx_r),
     >                  (A(indx)-A(indx_r))
        end if  
      end do
      
*     **** k1=0, k2=-(nyh-1),nyh-1, k3=-(nzh-1),nzh-1 ****
      k1 = 0
      do k3=(-nzh+1),(nzh-1)
      do k2=1,(nyh-1)
         k  =  k3
         kr = -k3
         if (k.lt.0)  k  = k  + nz
         if (kr.lt.0) kr = kr + nz
         indx   = k *(nxh+1)*ny +  k2    *(nxh+1) + 1
         indx_r = kr*(nxh+1)*ny + (ny-k2)*(nxh+1) + 1
         if (A(indx).ne.dconjg(A(indx_r))) then
            write(*,*) "NR:",k1,k2,k3,A(indx),A(indx_r),
     >                  (A(indx)-dconjg(A(indx_r)))
        end if  
      end do
      end do

      return
      end
