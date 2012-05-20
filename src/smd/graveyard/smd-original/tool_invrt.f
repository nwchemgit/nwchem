      SUBROUTINE tool_invrt(latt,rlatt,det)
      
      implicit none

      real*8 latt,rlatt,det

      dimension latt(3,3),rlatt(3,3)
      
      rlatt(1,1)=latt(2,2)*latt(3,3)-latt(3,2)*latt(2,3)
      rlatt(2,1)=latt(3,1)*latt(2,3)-latt(2,1)*latt(3,3)
      rlatt(3,1)=latt(2,1)*latt(3,2)-latt(3,1)*latt(2,2)
      rlatt(1,2)=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      rlatt(2,2)=latt(1,1)*latt(3,3)-latt(3,1)*latt(1,3)
      rlatt(3,2)=latt(3,1)*latt(1,2)-latt(1,1)*latt(3,2)
      rlatt(1,3)=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)
      rlatt(2,3)=latt(2,1)*latt(1,3)-latt(1,1)*latt(2,3)
      rlatt(3,3)=latt(1,1)*latt(2,2)-latt(2,1)*latt(1,2)
      
      det=latt(1,1)*rlatt(1,1)+latt(1,2)*rlatt(2,1)+latt(1,3)*rlatt(3,1)
      if(abs(det).gt.0.d0)det=1.d0/det
      
      rlatt(1,1)=det*rlatt(1,1)
      rlatt(2,1)=det*rlatt(2,1)
      rlatt(3,1)=det*rlatt(3,1)
      rlatt(1,2)=det*rlatt(1,2)
      rlatt(2,2)=det*rlatt(2,2)
      rlatt(3,2)=det*rlatt(3,2)
      rlatt(1,3)=det*rlatt(1,3)
      rlatt(2,3)=det*rlatt(2,3)
      rlatt(3,3)=det*rlatt(3,3)

      return

      end
c $Id$
