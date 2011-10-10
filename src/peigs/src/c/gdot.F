      double precision function gdot(n, buf1, nx, buf2, ny)
c     
c kahan's careful sum applied to dot products
c
      integer n, nx, ny
      double precision buf1(*), buf2(*)
c     
      integer i
      double precision s_1, c_1, y_1, t_1, f_1
c     
      s_1 = 0.d0
      c_1 = 0.d0
      do i = 1, n
         y_1 = c_1 + buf1(i)*buf2(i)
         t_1 = s_1 + y_1
         f_1 = 0.d0
         if ( y_1 .lt. 0.d0 .and. s_1 .lt. 0.d0 )
     $        f_1 = (0.46d0*t_1 - t_1)+t_1
         if ( y_1 .gt. 0.d0 .and. s_1 .gt. 0.d0 )
     $        f_1 = (0.46d0*t_1 - t_1)+t_1
         c_1 = ((s_1 - f_1) - ( t_1 - f_1)) + y_1
         s_1 = t_1
      enddo
      s_1 = s_1 + c_1
      gdot = s_1
      return
      end
c $Id$
