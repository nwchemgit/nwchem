C $Id: euler.F,v 1.2 2001-01-08 15:51:59 windus Exp $
      Subroutine drdy_euler (n3tm,n3,x,dx,dels)
      Implicit None
      Integer n3tm, n3
      Double Precision x(n3tm),dx(n3tm), dels
*
      Double Precision dxmag
      Integer i
*
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
