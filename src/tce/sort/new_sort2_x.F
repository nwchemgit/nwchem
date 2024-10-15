!
! These are just matrix transposes.  The code should be much simpler.
!
      subroutine tce_sort_2_x(u,s,a,b,i,j,f)
      implicit none
      integer, intent(in) :: a,b,i,j
      double precision, intent(out) :: s(a*b)
      double precision, intent(in) :: u(a*b)
      double precision, intent(in) :: f
      integer :: ia,ib,j1,j2
      ia = 0
      if (i.eq.2) then ! transpose
       do concurrent (j1 = 1:a, j2 = 1:b)
         ia = j2+b*(j1-1)
         ib = j1+a*(j2-1)
         s(ib) = u(ia) * f
       end do
      else ! i.eq.1 copy
       do concurrent (ia = 1:a*b)
         s(ia) = u(ia) * f
       end do
      endif
      end
