*
* $Id$
*

      subroutine lattice_min_difference(x,y,z)
      implicit none
      real*8 x,y,z

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

      real*8 ub(3,3)
      common / lattice_block2 / ub

*     *** local variables ****
      real*8 c1,c2,c3

      c1 = x*ub(1,1) + y*ub(2,1) + z*ub(3,1)
      c2 = x*ub(1,2) + y*ub(2,2) + z*ub(3,2)
      c3 = x*ub(1,3) + y*ub(2,3) + z*ub(3,3)
      c1 = c1 - DNINT(c1)
      c2 = c2 - DNINT(c2)
      c3 = c3 - DNINT(c3)
      x = ua(1,1)*c1 + ua(1,2)*c2 + ua(1,3)*c3
      y = ua(2,1)*c1 + ua(2,2)*c2 + ua(2,3)*c3
      z = ua(3,1)*c1 + ua(3,2)*c2 + ua(3,3)*c3

      return
      end

      subroutine lattice_frac_to_r1(n,f1,r1)
      implicit none
      integer n
      real*8 f1(3,*),r1(3,*)

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

*     **** local variables ***
      integer i

      do i=1,n
         r1(1,i) = ua(1,1)*f1(1,i) + ua(1,2)*f1(2,i) + ua(1,3)*f1(3,i)
         r1(2,i) = ua(2,1)*f1(1,i) + ua(2,2)*f1(2,i) + ua(2,3)*f1(3,i)
         r1(3,i) = ua(3,1)*f1(1,i) + ua(3,2)*f1(2,i) + ua(3,3)*f1(3,i)
      end do
      return
      end

      subroutine lattice_r1_to_frac(n,r1,f1)
      implicit none
      integer n
      real*8 r1(3,*),f1(3,*)

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

*     **** local variables ***
      integer i

      do i=1,n
         f1(1,i) = r1(1,i)*ub(1,1) + r1(2,i)*ub(2,1) + r1(3,i)*ub(3,1)
         f1(2,i) = r1(1,i)*ub(1,2) + r1(2,i)*ub(2,2) + r1(3,i)*ub(3,2)
         f1(3,i) = r1(1,i)*ub(1,3) + r1(2,i)*ub(2,3) + r1(3,i)*ub(3,3)
      end do
      return
      end


      subroutine lattice_center_xyz_to_ijk(x,y,z,c1,c2,c3)
      implicit none
      real*8 x,y,z
      integer c1,c2,c3

      integer np1,np2,np3
      real*8  f1,f2,f3

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      f1 = x*ub(1,1) + y*ub(2,1) + z*ub(3,1)
      f2 = x*ub(1,2) + y*ub(2,2) + z*ub(3,2)
      f3 = x*ub(1,3) + y*ub(2,3) + z*ub(3,3)
      c1 = (f1-0.5d0+0.5d0/dble(np1))*np1
      c2 = (f2-0.5d0+0.5d0/dble(np2))*np2 
      c3 = (f3-0.5d0+0.5d0/dble(np3))*np3 
      return
      end 

      subroutine lattice_center0_xyz_to_ijk(x,y,z,c1,c2,c3)
      implicit none
      real*8 x,y,z
      integer c1,c2,c3

      integer np1,np2,np3
      real*8  f1,f2,f3

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      f1 = x*ub(1,1) + y*ub(2,1) + z*ub(3,1)
      f2 = x*ub(1,2) + y*ub(2,2) + z*ub(3,2)
      f3 = x*ub(1,3) + y*ub(2,3) + z*ub(3,3)
      c1 = (f1+0.5d0/dble(np1))*np1
      c2 = (f2+0.5d0/dble(np2))*np2 
      c3 = (f3+0.5d0/dble(np3))*np3 
      return
      end 


      subroutine lattice_fragcell(n1,r1)
      implicit none
      real*8  rcm(3)
      integer n1
      real*8  r1(3,*)

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

*     **** local variables ****
      integer i,j
      real*8 fcm

      do i=2,n1
         do j=1,3
            fcm = (r1(1,i)-r1(1,1))*ub(1,j) 
     >          + (r1(2,i)-r1(2,1))*ub(2,j) 
     >          + (r1(3,i)-r1(3,1))*ub(3,j)
            do while (fcm.gt.0.5)
               r1(1,i) = r1(1,i) - ua(1,j)
               r1(2,i) = r1(2,i) - ua(2,j)
               r1(3,i) = r1(3,i) - ua(3,j)
               fcm = (r1(1,i)-r1(1,1))*ub(1,j) 
     >             + (r1(2,i)-r1(2,1))*ub(2,j) 
     >             + (r1(3,i)-r1(3,1))*ub(3,j)
            end do
            do while (fcm.le.(-0.5))
               r1(1,i) = r1(1,i) + ua(1,j)
               r1(2,i) = r1(2,i) + ua(2,j)
               r1(3,i) = r1(3,i) + ua(3,j)
               fcm = (r1(1,i)-r1(1,1))*ub(1,j) 
     >             + (r1(2,i)-r1(2,1))*ub(2,j) 
     >             + (r1(3,i)-r1(3,1))*ub(3,j)
            end do
         end do
      end do
      return
      end


      subroutine lattice_incell1_frag(rcm,n1,r1)
      implicit none 
      real*8  rcm(3)
      integer n1
      real*8  r1(3,*)

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

*     **** local variables ****
      integer i
      real*8 fcm(3)

      call lattice_r1_to_frac(1,rcm,fcm)
      
      do while (fcm(1).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,1)
         rcm(2) = rcm(2) - ua(2,1)
         rcm(3) = rcm(3) - ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,1)
           r1(2,i) = r1(2,i) - ua(2,1)
           r1(3,i) = r1(3,i) - ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do
      do while (fcm(1).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,1)
         rcm(2) = rcm(2) + ua(2,1)
         rcm(3) = rcm(3) + ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,1)
           r1(2,i) = r1(2,i) + ua(2,1)
           r1(3,i) = r1(3,i) + ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do

      do while (fcm(2).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,2)
         rcm(2) = rcm(2) - ua(2,2)
         rcm(3) = rcm(3) - ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,2)
           r1(2,i) = r1(2,i) - ua(2,2)
           r1(3,i) = r1(3,i) - ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do
      do while (fcm(2).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,2)
         rcm(2) = rcm(2) + ua(2,2)
         rcm(3) = rcm(3) + ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,2)
           r1(2,i) = r1(2,i) + ua(2,2)
           r1(3,i) = r1(3,i) + ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do

      do while (fcm(3).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,3)
         rcm(2) = rcm(2) - ua(2,3)
         rcm(3) = rcm(3) - ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,3)
           r1(2,i) = r1(2,i) - ua(2,3)
           r1(3,i) = r1(3,i) - ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      do while (fcm(3).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,3)
         rcm(2) = rcm(2) + ua(2,3)
         rcm(3) = rcm(3) + ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,3)
           r1(2,i) = r1(2,i) + ua(2,3)
           r1(3,i) = r1(3,i) + ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      return
      end

      subroutine lattice_incell2_frag(rcm,n1,r1,r2)
      implicit none 
      real*8  rcm(3)
      integer n1
      real*8  r1(3,*)
      real*8  r2(3,*)

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

*     **** local variables ****
      integer i
      real*8 fcm(3)

      call lattice_r1_to_frac(1,rcm,fcm)
      
      do while (fcm(1).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,1)
         rcm(2) = rcm(2) - ua(2,1)
         rcm(3) = rcm(3) - ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,1)
           r1(2,i) = r1(2,i) - ua(2,1)
           r1(3,i) = r1(3,i) - ua(3,1)

           r2(1,i) = r2(1,i) - ua(1,1)
           r2(2,i) = r2(2,i) - ua(2,1)
           r2(3,i) = r2(3,i) - ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do
      do while (fcm(1).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,1)
         rcm(2) = rcm(2) + ua(2,1)
         rcm(3) = rcm(3) + ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,1)
           r1(2,i) = r1(2,i) + ua(2,1)
           r1(3,i) = r1(3,i) + ua(3,1)

           r2(1,i) = r2(1,i) + ua(1,1)
           r2(2,i) = r2(2,i) + ua(2,1)
           r2(3,i) = r2(3,i) + ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do

      do while (fcm(2).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,2)
         rcm(2) = rcm(2) - ua(2,2)
         rcm(3) = rcm(3) - ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,2)
           r1(2,i) = r1(2,i) - ua(2,2)
           r1(3,i) = r1(3,i) - ua(3,2)

           r2(1,i) = r2(1,i) - ua(1,2)
           r2(2,i) = r2(2,i) - ua(2,2)
           r2(3,i) = r2(3,i) - ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do
      do while (fcm(2).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,2)
         rcm(2) = rcm(2) + ua(2,2)
         rcm(3) = rcm(3) + ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,2)
           r1(2,i) = r1(2,i) + ua(2,2)
           r1(3,i) = r1(3,i) + ua(3,2)

           r2(1,i) = r2(1,i) + ua(1,2)
           r2(2,i) = r2(2,i) + ua(2,2)
           r2(3,i) = r2(3,i) + ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do

      do while (fcm(3).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,3)
         rcm(2) = rcm(2) - ua(2,3)
         rcm(3) = rcm(3) - ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,3)
           r1(2,i) = r1(2,i) - ua(2,3)
           r1(3,i) = r1(3,i) - ua(3,3)

           r2(1,i) = r2(1,i) - ua(1,3)
           r2(2,i) = r2(2,i) - ua(2,3)
           r2(3,i) = r2(3,i) - ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      do while (fcm(3).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,3)
         rcm(2) = rcm(2) + ua(2,3)
         rcm(3) = rcm(3) + ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,3)
           r1(2,i) = r1(2,i) + ua(2,3)
           r1(3,i) = r1(3,i) + ua(3,3)

           r2(1,i) = r2(1,i) + ua(1,3)
           r2(2,i) = r2(2,i) + ua(2,3)
           r2(3,i) = r2(3,i) + ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      return
      end

      subroutine lattice_incell3_frag(rcm,n1,r1,r2,r3)
      implicit none 
      real*8  rcm(3)
      integer n1
      real*8  r1(3,*)
      real*8  r2(3,*)
      real*8  r3(3,*)

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 ua(3,3),unitg(3,3)
      common / lattice_block / ua,unitg,ecut,wcut,omega

*     **** common block ****
      real*8 ub(3,3)
      common / lattice_block2 / ub

*     **** local variables ****
      integer i
      real*8 fcm(3)

      call lattice_r1_to_frac(1,rcm,fcm)
      
      do while (fcm(1).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,1)
         rcm(2) = rcm(2) - ua(2,1)
         rcm(3) = rcm(3) - ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,1)
           r1(2,i) = r1(2,i) - ua(2,1)
           r1(3,i) = r1(3,i) - ua(3,1)

           r2(1,i) = r2(1,i) - ua(1,1)
           r2(2,i) = r2(2,i) - ua(2,1)
           r2(3,i) = r2(3,i) - ua(3,1)

           r3(1,i) = r3(1,i) - ua(1,1)
           r3(2,i) = r3(2,i) - ua(2,1)
           r3(3,i) = r3(3,i) - ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do
      do while (fcm(1).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,1)
         rcm(2) = rcm(2) + ua(2,1)
         rcm(3) = rcm(3) + ua(3,1)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,1)
           r1(2,i) = r1(2,i) + ua(2,1)
           r1(3,i) = r1(3,i) + ua(3,1)

           r2(1,i) = r2(1,i) + ua(1,1)
           r2(2,i) = r2(2,i) + ua(2,1)
           r2(3,i) = r2(3,i) + ua(3,1)

           r3(1,i) = r3(1,i) + ua(1,1)
           r3(2,i) = r3(2,i) + ua(2,1)
           r3(3,i) = r3(3,i) + ua(3,1)
         end do
         fcm(1) = rcm(1)*ub(1,1) + rcm(2)*ub(2,1) + rcm(3)*ub(3,1)
      end do

      do while (fcm(2).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,2)
         rcm(2) = rcm(2) - ua(2,2)
         rcm(3) = rcm(3) - ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,2)
           r1(2,i) = r1(2,i) - ua(2,2)
           r1(3,i) = r1(3,i) - ua(3,2)

           r2(1,i) = r2(1,i) - ua(1,2)
           r2(2,i) = r2(2,i) - ua(2,2)
           r2(3,i) = r2(3,i) - ua(3,2)

           r3(1,i) = r3(1,i) - ua(1,2)
           r3(2,i) = r3(2,i) - ua(2,2)
           r3(3,i) = r3(3,i) - ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do
      do while (fcm(2).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,2)
         rcm(2) = rcm(2) + ua(2,2)
         rcm(3) = rcm(3) + ua(3,2)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,2)
           r1(2,i) = r1(2,i) + ua(2,2)
           r1(3,i) = r1(3,i) + ua(3,2)

           r2(1,i) = r2(1,i) + ua(1,2)
           r2(2,i) = r2(2,i) + ua(2,2)
           r2(3,i) = r2(3,i) + ua(3,2)

           r3(1,i) = r3(1,i) + ua(1,2)
           r3(2,i) = r3(2,i) + ua(2,2)
           r3(3,i) = r3(3,i) + ua(3,2)
         end do
         fcm(2) = rcm(1)*ub(1,2) + rcm(2)*ub(2,2) + rcm(3)*ub(3,2)
      end do

      do while (fcm(3).gt.(0.5d0))
         rcm(1) = rcm(1) - ua(1,3)
         rcm(2) = rcm(2) - ua(2,3)
         rcm(3) = rcm(3) - ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) - ua(1,3)
           r1(2,i) = r1(2,i) - ua(2,3)
           r1(3,i) = r1(3,i) - ua(3,3)

           r2(1,i) = r2(1,i) - ua(1,3)
           r2(2,i) = r2(2,i) - ua(2,3)
           r2(3,i) = r2(3,i) - ua(3,3)

           r3(1,i) = r3(1,i) - ua(1,3)
           r3(2,i) = r3(2,i) - ua(2,3)
           r3(3,i) = r3(3,i) - ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      do while (fcm(3).le.(-0.5d0))
         rcm(1) = rcm(1) + ua(1,3)
         rcm(2) = rcm(2) + ua(2,3)
         rcm(3) = rcm(3) + ua(3,3)
         do i=1,n1
           r1(1,i) = r1(1,i) + ua(1,3)
           r1(2,i) = r1(2,i) + ua(2,3)
           r1(3,i) = r1(3,i) + ua(3,3)

           r2(1,i) = r2(1,i) + ua(1,3)
           r2(2,i) = r2(2,i) + ua(2,3)
           r2(3,i) = r2(3,i) + ua(3,3)

           r3(1,i) = r3(1,i) + ua(1,3)
           r3(2,i) = r3(2,i) + ua(2,3)
           r3(3,i) = r3(3,i) + ua(3,3)
         end do
         fcm(3) = rcm(1)*ub(1,3) + rcm(2)*ub(2,3) + rcm(3)*ub(3,3)
      end do
      return
      end


      real*8 function lattice_wcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_wcut = wcut
      return
      end

      real*8 function lattice_ecut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_ecut = ecut
      return
      end

      real*8 function lattice_ggcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_ggcut = 2.0d0*ecut
      return
      end

      real*8 function lattice_wggcut()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_wggcut = 2.0d0*wcut
      return
      end



      real*8 function lattice_omega()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_omega = omega
      return
      end

      real*8 function lattice_unita(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_unita = unita(i,j)
      return
      end


      real*8 function lattice_unitg(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      lattice_unitg = unitg(i,j)
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_unita_small             *
*     *                                          *
*     ********************************************
      real*8 function lattice_unita_small(i,j)
      implicit none
      integer i,j

*     **** common block ****
      logical has_small
      real*8  omega_small
      real*8  unita_small(3,3),unitg_small(3,3)
      common /lattice_small_block/ unita_small,unitg_small,
     >                             omega_small,has_small

      lattice_unita_small = unita_small(i,j)
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_unita_frozen_small      *
*     *                                          *
*     ********************************************
      real*8 function lattice_unita_frozen_small(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8  omega_frozen_small
      real*8  unita_frozen_small(3,3),unitg_frozen_small(3,3)
      common /lattice_small_frozen_block/ unita_frozen_small,
     >                             unitg_frozen_small,
     >                             omega_frozen_small

      lattice_unita_frozen_small = unita_frozen_small(i,j)
      return
      end


*     ********************************************
*     *                                          *
*     *          lattice_unitg_small             *
*     *                                          *
*     ********************************************
      real*8 function lattice_unitg_small(i,j)
      implicit none
      integer i,j

*     **** common block ****
      logical has_small
      real*8  omega_small
      real*8  unita_small(3,3),unitg_small(3,3)
      common /lattice_small_block/ unita_small,unitg_small,
     >                             omega_small,has_small

      lattice_unitg_small = unitg_small(i,j)
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_unitg_frozen_small      *
*     *                                          *
*     ********************************************
      real*8 function lattice_unitg_frozen_small(i,j)
      implicit none
      integer i,j

*     **** common block ****
      real*8  omega_frozen_small
      real*8  unita_frozen_small(3,3),unitg_frozen_small(3,3)
      common /lattice_small_frozen_block/ unita_frozen_small,
     >                             unitg_frozen_small,
     >                             omega_frozen_small

      lattice_unitg_frozen_small = unitg_frozen_small(i,j)
      return
      end



*     ********************************************
*     *                                          *
*     *          lattice_omega_small             *
*     *                                          *
*     ********************************************
      real*8 function lattice_omega_small(i,j)
      implicit none
      integer i,j

*     **** common block ****
      logical has_small
      real*8  omega_small
      real*8  unita_small(3,3),unitg_small(3,3)
      common /lattice_small_block/ unita_small,unitg_small,
     >                             omega_small,has_small

      lattice_omega_small = omega_small
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_omega_frozen_small      *
*     *                                          *
*     ********************************************
      real*8 function lattice_omega_frozen_small()
      implicit none
      integer i,j

*     **** common block ****
      real*8  omega_frozen_small
      real*8  unita_frozen_small(3,3),unitg_frozen_small(3,3)
      common /lattice_small_frozen_block/ unita_frozen_small,
     >                             unitg_frozen_small,
     >                             omega_frozen_small

      lattice_omega_frozen_small = omega_frozen_small
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_has_small               *
*     *                                          *
*     ********************************************
      logical function lattice_has_small()
      implicit none

*     **** common block ****
      logical has_small
      real*8  omega_small
      real*8  unita_small(3,3),unitg_small(3,3)
      common /lattice_small_block/ unita_small,unitg_small,
     >                             omega_small,has_small

      lattice_has_small = has_small
      return
      end





*     ********************************************
*     *                                          *
*     *          lattice_unitg_frozen            *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_unitg_frozen(i,j)
      implicit none
      integer i,j

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_unitg_frozen = unitg_frozen(i,j)
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_unita_frozen            *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_unita_frozen(i,j)
      implicit none
      integer i,j

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_unita_frozen = unita_frozen(i,j)
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_omega_frozen            *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_omega_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_omega_frozen = omega_frozen
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_wcut_frozen             *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_wcut_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_wcut_frozen = wcut_frozen
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_ecut_frozen             *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_ecut_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_ecut_frozen = ecut_frozen
      return
      end

*     ********************************************
*     *                                          *
*     *          lattice_ggcut_frozen            *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_ggcut_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_ggcut_frozen = 2.0d0*ecut_frozen
      return
      end

*     ********************************************
*     *                                          *
*     *         lattice_wggcut_frozen            *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      real*8 function lattice_wggcut_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_wggcut_frozen = 2.0d0*wcut_frozen
      return
      end

*     ********************************************
*     *                                          *
*     *         lattice_frozen                   *
*     *                                          *
*     ********************************************

*     frozen lattice structure - used to determine masking arrays

      logical function lattice_frozen()
      implicit none

*     **** common blocks ****
      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      lattice_frozen = frozen
      return
      end




*     *********************************************
*     *                                           *
*     *               lattice_init                *
*     *                                           *
*     *********************************************

      subroutine lattice_init()
      implicit none

*     **** common block ****
      real*8 ecut,wcut,omega
      real*8 unita(3,3),unitg(3,3)
      common / lattice_block / unita,unitg,ecut,wcut,omega

      real*8 ub(3,3)
      common / lattice_block2 / ub

      logical frozen
      real*8 ecut_frozen,wcut_frozen,omega_frozen
      real*8 unita_frozen(3,3),unitg_frozen(3,3)
      common / lattice_froze / unita_frozen,unitg_frozen,
     >                         ecut_frozen,wcut_frozen,
     >                         omega_frozen,frozen

      logical has_small
      real*8  omega_small
      real*8  unita_small(3,3),unitg_small(3,3)
      common /lattice_small_block/ unita_small,unitg_small,
     >                             omega_small,has_small

      real*8  omega_frozen_small
      real*8  unita_frozen_small(3,3),unitg_frozen_small(3,3)
      common /lattice_small_frozen_block/ unita_frozen_small,
     >                             unitg_frozen_small,
     >                             omega_frozen_small


*     **** local variables ****
      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8  gx,gy,gz,gg
      real*8  gg1,gg2,gg3
      real*8  ecut0,wcut0

*     **** external functions ****
      logical  control_frozen,control_has_ngrid_small
      integer  control_ngrid,control_ngrid_small
      real*8   control_unita,control_ecut,control_wcut
      real*8   control_unita_frozen
      external control_frozen,control_has_ngrid_small
      external control_ngrid,control_ngrid_small
      external control_unita,control_ecut,control_wcut
      external control_unita_frozen
        
      ecut0 = control_ecut()
      wcut0 = control_wcut()

*     **** define lattice ****
      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)
      call get_cube(unita,unitg,omega)

*     **** define frozen lattice - note this is only different from unita if nwpw:frozen_lattice is set ****
      frozen = control_frozen()
      unita_frozen(1,1) = control_unita_frozen(1,1)
      unita_frozen(2,1) = control_unita_frozen(2,1)
      unita_frozen(3,1) = control_unita_frozen(3,1)
      unita_frozen(1,2) = control_unita_frozen(1,2)
      unita_frozen(2,2) = control_unita_frozen(2,2)
      unita_frozen(3,2) = control_unita_frozen(3,2)
      unita_frozen(1,3) = control_unita_frozen(1,3)
      unita_frozen(2,3) = control_unita_frozen(2,3)
      unita_frozen(3,3) = control_unita_frozen(3,3)
      call get_cube(unita_frozen,unitg_frozen,omega_frozen)

*     **** define ub ****
      call dcopy(9,unitg,1,ub,1)
      call dscal(9,(1.0d0/(8.0d0*datan(1.0d0))),ub,1)
       


*     *** set the ecut variable using the frozen lattice***
c     call D3dB_nx(1,nx)
c     call D3dB_ny(1,ny)
c     call D3dB_nz(1,nz)
      nx = control_ngrid(1)
      ny = control_ngrid(2)
      nz = control_ngrid(3)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      gx = unitg_frozen(1,1)*dble(nxh)
      gy = unitg_frozen(2,1)*dble(nxh)
      gz = unitg_frozen(3,1)*dble(nxh)
      gg1 = gx*gx + gy*gy + gz*gz

      gx = unitg_frozen(1,2)*dble(nyh)
      gy = unitg_frozen(2,2)*dble(nyh)
      gz = unitg_frozen(3,2)*dble(nyh)
      gg2 = gx*gx + gy*gy + gz*gz

      gx = unitg_frozen(1,3)*dble(nzh)
      gy = unitg_frozen(2,3)*dble(nzh)
      gz = unitg_frozen(3,3)*dble(nzh)
      gg3 = gx*gx + gy*gy + gz*gz

      gg = gg1
      if (gg2.lt.gg) gg=gg2
      if (gg3.lt.gg) gg=gg3

      ecut = 0.5d0*gg
      if (ecut0.lt.ecut) then
         ecut = ecut0
      end if

      wcut = ecut
      if (wcut0.lt.wcut) then
         wcut = wcut0
      end if

*     **** set ecut,wcut for frozen lattice ****
      wcut_frozen = wcut
      ecut_frozen = ecut


*     **** small cell ****
      has_small = control_has_ngrid_small()
      if (has_small) then
         gg1 = dble(control_ngrid_small(1))/dble(nx)
         gg2 = dble(control_ngrid_small(2))/dble(ny)
         gg3 = dble(control_ngrid_small(3))/dble(nz)
         unita_small(1,1) = unita(1,1)*gg1
         unita_small(2,1) = unita(2,1)*gg1
         unita_small(3,1) = unita(3,1)*gg1
         unita_small(1,2) = unita(1,2)*gg2
         unita_small(2,2) = unita(2,2)*gg2
         unita_small(3,2) = unita(3,2)*gg2
         unita_small(1,3) = unita(1,3)*gg3
         unita_small(2,3) = unita(2,3)*gg3
         unita_small(3,3) = unita(3,3)*gg3
         call get_cube(unita_small,unitg_small,omega_small)
         unita_frozen_small(1,1) = unita_frozen(1,1)*gg1
         unita_frozen_small(2,1) = unita_frozen(2,1)*gg1
         unita_frozen_small(3,1) = unita_frozen(3,1)*gg1
         unita_frozen_small(1,2) = unita_frozen(1,2)*gg2
         unita_frozen_small(2,2) = unita_frozen(2,2)*gg2
         unita_frozen_small(3,2) = unita_frozen(3,2)*gg2
         unita_frozen_small(1,3) = unita_frozen(1,3)*gg3
         unita_frozen_small(2,3) = unita_frozen(2,3)*gg3
         unita_frozen_small(3,3) = unita_frozen(3,3)*gg3
         call get_cube(unita_frozen_small,
     >                 unitg_frozen_small,
     >                 omega_frozen_small)
      end if

      return
      end


*     *******************************
*     *                             *
*     *         lattice_p_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          D3dB_nx --- number of grid points in direction 1
*          D3dB_ny --- number of grid points in direction 2
*          D3dB_nz --- number of grid points in direction 2
*          lattice_unita -- primitive lattice vectors in real space
*
*     Exit -
*          xs  --- coordinates of grid points (Rx,Ry,Rz)
*
*
      subroutine lattice_p_grid(xs)
      implicit none
      real*8 xs(*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer idx,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),dk1,dk2,dk3,twopi


*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      twopi = 8.0d0*datan(1.0d0)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      call dcopy(9,0.0d0,0,a,1)
      a(1,1) = twopi/dble(np1)
      a(2,2) = twopi/dble(np2)
      a(3,3) = twopi/dble(np3)

      call dcopy(6*n2ft3d,0.0d0,0,xs,1)

*     **** grid points in coordination space ****
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call D3dB_ijktoindex2p(1,i+1,j+1,k+1,idx,p)
               if (p .eq. taskid) then
                dk1=dble(k1)
                dk2=dble(k2)
                dk3=dble(k3)
                xs(idx)         =dcos(a(1,1)*dk1+a(1,2)*dk2+a(1,3)*dk3)
                xs(idx+n2ft3d)  =dsin(a(1,1)*dk1+a(1,2)*dk2+a(1,3)*dk3)
                xs(idx+2*n2ft3d)=dcos(a(2,1)*dk1+a(2,2)*dk2+a(2,3)*dk3)
                xs(idx+3*n2ft3d)=dsin(a(2,1)*dk1+a(2,2)*dk2+a(2,3)*dk3)
                xs(idx+4*n2ft3d)=dcos(a(3,1)*dk1+a(3,2)*dk2+a(3,3)*dk3)
                xs(idx+5*n2ft3d)=dsin(a(3,1)*dk1+a(3,2)*dk2+a(3,3)*dk3)
               end if
          end do
        end do
      end do

      return
      end




*     *******************************
*     *                             *
*     *         lattice_r_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          D3dB_nx --- number of grid points in direction 1
*          D3dB_ny --- number of grid points in direction 2
*          D3dB_nz --- number of grid points in direction 2
*          lattice_unita -- primitive lattice vectors in real space
*
*     Exit -
*          r  --- coordinates of grid points (Rx,Ry,Rz)
*
*
      subroutine lattice_r_grid(r)
      implicit none
      real*8 r(3,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid,tid,nthreads
      integer index,k1,k2,k3,it
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),dk1,dk2,dk3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/dble(np1)
         a(i,2) = lattice_unita(i,2)/dble(np2)
         a(i,3) = lattice_unita(i,3)/dble(np3)
      end do

      !call dcopy(3*n2ft3d,0.0d0,0,r,1)
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,r)

*     **** grid points in coordination space ****
      it = 0
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
               if ((p.eq.taskid).and.(it.eq.tid)) then
c                 index = (q-1)*(np1+2)*np2
c    >                  + j    *(np1+2)
c    >                  + i+1
                dk1=dble(k1)
                dk2=dble(k2)
                dk3=dble(k3)
                r(1,index) = a(1,1)*dk1 + a(1,2)*dk2 + a(1,3)*dk3
                r(2,index) = a(2,1)*dk1 + a(2,2)*dk2 + a(2,3)*dk3
                r(3,index) = a(3,1)*dk1 + a(3,2)*dk2 + a(3,3)*dk3

c*               **** reverse y and z ****
c                r(1,index) = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
c                r(2,index) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
c                r(3,index) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

               end if
               it = mod(it+1,nthreads)
          end do
        end do
      end do

      return
      end

*     *******************************
*     *                             *
*     *         lattice_r_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          D3dB_nx --- number of grid points in direction 1
*          D3dB_ny --- number of grid points in direction 2
*          D3dB_nz --- number of grid points in direction 2
*          lattice_unita -- primitive lattice vectors in real space
*
*     Exit -
*          r  --- coordinates of grid points (Rx,Ry,Rz)
*
*
      subroutine c_lattice_r_grid(r)
      implicit none
      real*8 r(3,*)

*     **** local variables ****
      integer n2ft3d
      integer i,j,k,p,taskid,tid,nthreads
      integer index,k1,k2,k3,it
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),dk1,dk2,dk3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


*     **** constants ****
      call Parallel3d_taskid_i(taskid)
      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()
      call C3dB_n2ft3d(1,n2ft3d)
      call C3dB_nx(1,np1)
      call C3dB_ny(1,np2)
      call C3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/dble(np1)
         a(i,2) = lattice_unita(i,2)/dble(np2)
         a(i,3) = lattice_unita(i,3)/dble(np3)
      end do

      !call dcopy(3*n2ft3d,0.0d0,0,r,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,r)

*     **** grid points in coordination space ****
      it = 0
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               !call D3dB_ktoqp(1,k+1,q,p)
               call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
               if ((p.eq.taskid).and.(it.eq.tid)) then
c                 index = (q-1)*(np1+2)*np2
c    >                  + j    *(np1+2)
c    >                  + i+1
                dk1=dble(k1)
                dk2=dble(k2)
                dk3=dble(k3)
                r(1,index) = a(1,1)*dk1 + a(1,2)*dk2 + a(1,3)*dk3
                r(2,index) = a(2,1)*dk1 + a(2,2)*dk2 + a(2,3)*dk3
                r(3,index) = a(3,1)*dk1 + a(3,2)*dk2 + a(3,3)*dk3

c*               **** reverse y and z ****
c                r(1,index) = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
c                r(2,index) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
c                r(3,index) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

               end if
               it = mod(it+1,nthreads)
          end do
        end do
      end do

      return
      end



      subroutine lattice_r_grid_sym(r)
      implicit none
      real*8 r(3,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),dk1,dk2,dk3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/dble(np1)
         a(i,2) = lattice_unita(i,2)/dble(np2)
         a(i,3) = lattice_unita(i,3)/dble(np3)
      end do

      call dcopy(3*n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3+1, nph3-1
        do k2 = -nph2+1, nph2-1
          do k1 = -nph1+1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
c                 index = (q-1)*(np1+2)*np2
c    >                  + j    *(np1+2)
c    >                  + i+1
                dk1=dble(k1)
                dk2=dble(k2)
                dk3=dble(k3)
                r(1,index) = a(1,1)*dk1 + a(1,2)*dk2 + a(1,3)*dk3
                r(2,index) = a(2,1)*dk1 + a(2,2)*dk2 + a(2,3)*dk3
                r(3,index) = a(3,1)*dk1 + a(3,2)*dk2 + a(3,3)*dk3

c*               **** reverse y and z ****
c                r(1,index) = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
c                r(2,index) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
c                r(3,index) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

               end if
          end do
        end do
      end do

      return
      end


      subroutine lattice_r_grid_sym0(r)
      implicit none
      real*8 r(*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),dk1,dk2,dk3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/dble(np1)
         a(i,2) = lattice_unita(i,2)/dble(np2)
         a(i,3) = lattice_unita(i,3)/dble(np3)
      end do

      call dcopy(3*n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3+1, nph3-1
        do k2 = -nph2+1, nph2-1
          do k1 = -nph1+1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
c                 index = (q-1)*(np1+2)*np2
c    >                  + j    *(np1+2)
c    >                  + i+1
                dk1=dble(k1)
                dk2=dble(k2)
                dk3=dble(k3)
                r(index)          = a(1,1)*dk1 + a(1,2)*dk2 + a(1,3)*dk3
                r(index+  n2ft3d) = a(2,1)*dk1 + a(2,2)*dk2 + a(2,3)*dk3
                r(index+2*n2ft3d) = a(3,1)*dk1 + a(3,2)*dk2 + a(3,3)*dk3
               end if
          end do
        end do
      end do

      return
      end


*     *******************************
*     *                             *
*     *         lattice_i_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          D3dB_nx --- number of grid points in direction 1
*          D3dB_ny --- number of grid points in direction 2
*          D3dB_nz --- number of grid points in direction 2
*
*     Entry - nb
*     Exit -
*          ijk --- coordinates of grid points (k1,k2,k3)
*
*
      subroutine lattice_i_grid(nb,ijk)
      implicit none
      integer nb
      integer ijk(4,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3

*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      call D3dB_nfft3d(nb,nfft3d)

      n2ft3d = 2*nfft3d
      call D3dB_nx(nb,np1)
      call D3dB_ny(nb,np2)
      call D3dB_nz(nb,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

      call icopy(4*n2ft3d,0,0,ijk,1)

*     **** grid points  ****
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call D3dB_ijktoindex2p(nb,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
                  ijk(1,index) = k1
                  ijk(2,index) = k2
                  ijk(3,index) = k3
                  ijk(4,index) = 1
               end if
          end do
        end do
      end do

      return
      end 

*     *******************************
*     *                             *
*     *       c_lattice_i_grid      *
*     *                             *
*     *******************************
*
*     This routine computes coordinates of grid points in
*     the unit cell
*
*     Uses -
*          Parallel_taskid --- processor number
*          C3dB_nx --- number of grid points in direction 1
*          C3dB_ny --- number of grid points in direction 2
*          C3dB_nz --- number of grid points in direction 2
*
*     Entry - nb
*     Exit -
*          ijk --- coordinates of grid points (k1,k2,k3)
*
*

      subroutine c_lattice_i_grid(nb,ijk)
      implicit none
      integer nb
      integer ijk(4,*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3

*     **** constants ****
      call Parallel3d_taskid_i(taskid)
      call C3dB_nfft3d(nb,nfft3d)

      call C3dB_nx(nb,np1)
      call C3dB_ny(nb,np2)
      call C3dB_nz(nb,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

      call icopy(4*nfft3d,0,0,ijk,1)

*     **** grid points  ****
      do k3 = -nph3, nph3-1
        do k2 = -nph2, nph2-1
          do k1 = -nph1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               call C3dB_ijktoindex2p(nb,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
                  ijk(1,index) = k1
                  ijk(2,index) = k2
                  ijk(3,index) = k3
                  ijk(4,index) = 1
               end if
          end do
        end do
      end do

      return
      end





      subroutine lattice_mask_sym(r)
      implicit none
      real*8 r(*)

*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3


*     **** constants ****
      call Parallel2d_taskid_i(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2


      call dcopy(n2ft3d,0.0d0,0,r,1)

*     **** grid points in coordination space ****
      do k3 = -nph3+1, nph3-1
        do k2 = -nph2+1, nph2-1
          do k1 = -nph1+1, nph1-1

               i = k1 + nph1
               j = k2 + nph2
               k = k3 + nph3

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
c                 index = (q-1)*(np1+2)*np2
c    >                  + j    *(np1+2)
c    >                  + i+1

                r(index) =  1.0d0
               end if
          end do
        end do
      end do

      return
      end




      subroutine get_cube(unita,unitg,volume)

******************************************************************************
*                                                                            *
*     This routine computes primitive vectors both in coordination           *
*     space and in reciporocal space and the volume of primitive cell.       *
*                                                                            *
*     Inputs:                                                                *
*             type --- type of cube (1=SC, 2=FCC, 3=BCC, 4=linear)           *
*             unit --- lattice constants                                     *
*                                                                            *
*     Outputs:                                                               *
*             volume --- volume of primitive cell                            *
*             unita  --- primitive vectors in coordination space             *
*             unitg  --- primitive vectors in reciprocal space               *
*                                                                            *
*     Library:  DSCAL from BLAS                                              *
*                                                                            *
*     Last modification:  7/03/93  by R. Kawai                               *
*                                                                            *
******************************************************************************

      implicit none

*     ------------------
*     argument variables
*     ------------------
      double precision unita(3,3), unitg(3,3)
      double precision volume

*     ---------------
*     local variables
*     ---------------
      double precision twopi

      twopi = 8.0d0*datan(1.0d0)


*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
*     ---------------------
*     volume of a unit cell
*     ---------------------
      volume = unita(1,1)*unitg(1,1)
     >       + unita(2,1)*unitg(2,1)
     >       + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/volume,unitg,1)

      volume=dabs(volume)
      return
      end

*     *******************************
*     *                             *
*     *     lattice_abc_abg         *
*     *                             *
*     *******************************
*
*     This routine computes a,b,c,alpha,beta,gamma.
*
      subroutine lattice_abc_abg(a,b,c,alpha,beta,gamma)
      implicit none
      real*8 a,b,c
      real*8 alpha,beta,gamma

*     *** local variables ****
      real*8 d2,pi

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      a = dsqrt(lattice_unita(1,1)**2 
     >        + lattice_unita(2,1)**2 
     >        + lattice_unita(3,1)**2)
      b = dsqrt(lattice_unita(1,2)**2 
     >        + lattice_unita(2,2)**2 
     >        + lattice_unita(3,2)**2)
      c = dsqrt(lattice_unita(1,3)**2 
     >        + lattice_unita(2,3)**2 
     >        + lattice_unita(3,3)**2)
 
      d2 = (lattice_unita(1,2)-lattice_unita(1,3))**2 
     >   + (lattice_unita(2,2)-lattice_unita(2,3))**2 
     >   + (lattice_unita(3,2)-lattice_unita(3,3))**2
      alpha = (b*b + c*c - d2)/(2.0d0*b*c)
      alpha = dacos(alpha)*180.0d0/pi
 
      d2 = (lattice_unita(1,3)-lattice_unita(1,1))**2 
     >   + (lattice_unita(2,3)-lattice_unita(2,1))**2 
     >   + (lattice_unita(3,3)-lattice_unita(3,1))**2
      beta = (c*c + a*a - d2)/(2.0d0*c*a)
      beta = dacos(beta)*180.0d0/pi
 
      d2 = (lattice_unita(1,1)-lattice_unita(1,2))**2 
     >   + (lattice_unita(2,1)-lattice_unita(2,2))**2 
     >   + (lattice_unita(3,1)-lattice_unita(3,2))**2
      gamma = (a*a + b*b - d2)/(2.0d0*a*b)
      gamma = dacos(gamma)*180.0d0/pi

      return
      end



*     *************************************
*     *                                   *
*     *     lattice_small_abc_abg         *
*     *                                   *
*     *************************************
*
*     This routine computes a,b,c,alpha,beta,gamma.
*
      subroutine lattice_small_abc_abg(a,b,c,alpha,beta,gamma)
      implicit none
      real*8 a,b,c
      real*8 alpha,beta,gamma

*     *** local variables ****
      real*8 d2,pi

*     **** external functions ****
      real*8   lattice_unita_small
      external lattice_unita_small

*     **** determine a,b,c,alpha,beta,gmma ***
      pi = 4.0d0*datan(1.0d0)
      a = dsqrt(lattice_unita_small(1,1)**2
     >        + lattice_unita_small(2,1)**2
     >        + lattice_unita_small(3,1)**2)
      b = dsqrt(lattice_unita_small(1,2)**2
     >        + lattice_unita_small(2,2)**2
     >        + lattice_unita_small(3,2)**2)
      c = dsqrt(lattice_unita_small(1,3)**2
     >        + lattice_unita_small(2,3)**2
     >        + lattice_unita_small(3,3)**2)

      d2 = (lattice_unita_small(1,2)-lattice_unita_small(1,3))**2
     >   + (lattice_unita_small(2,2)-lattice_unita_small(2,3))**2
     >   + (lattice_unita_small(3,2)-lattice_unita_small(3,3))**2
      alpha = (b*b + c*c - d2)/(2.0d0*b*c)
      alpha = dacos(alpha)*180.0d0/pi

      d2 = (lattice_unita_small(1,3)-lattice_unita_small(1,1))**2
     >   + (lattice_unita_small(2,3)-lattice_unita_small(2,1))**2
     >   + (lattice_unita_small(3,3)-lattice_unita_small(3,1))**2
      beta = (c*c + a*a - d2)/(2.0d0*c*a)
      beta = dacos(beta)*180.0d0/pi

      d2 = (lattice_unita_small(1,1)-lattice_unita_small(1,2))**2
     >   + (lattice_unita_small(2,1)-lattice_unita_small(2,2))**2
     >   + (lattice_unita_small(3,1)-lattice_unita_small(3,2))**2
      gamma = (a*a + b*b - d2)/(2.0d0*a*b)
      gamma = dacos(gamma)*180.0d0/pi

      return
      end




