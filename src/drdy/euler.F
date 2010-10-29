C $Id$
      Subroutine drdy_euler (n3tm,n3,x,dx,dels)
      Implicit None
      Integer n3tm,n3
      Double Precision x(n3tm),dx(n3tm),dels
C
      Integer i
      Double Precision dxmag
C
      dxmag = 0.0d00
      do i = 1,n3
         dxmag = dxmag + dx(i)*dx(i)
      enddo
      dxmag = sqrt(dxmag)
      do i = 1,n3
         x(i) = x(i) - dels*dx(i)/dxmag
      enddo
      return
      end
