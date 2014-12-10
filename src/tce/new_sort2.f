      subroutine tce_sort_2(u,s,a,b,i,j,f)
      implicit none
      integer a,b
      integer i,j
      integer j1,j2
      double precision s(b,a)
      double precision u(a,b)
      double precision f
      if (j.eq.1) then ! copy
      do j1 = 1,a
       do j2 = 1,b
        s(j2,j1) = u(j2,j1) * f
       enddo
      enddo
      else ! transpose
      do j1 = 1,a
       do j2 = 1,b
        s(j2,j1) = u(j1,j2) * f
       enddo
      enddo
      endif
      return
      end

      subroutine tce_sortacc_2(u,s,a,b,i,j,f)
      implicit none
      integer a,b
      integer i,j
      integer j1,j2
      double precision s(b,a)
      double precision u(a,b)
      double precision f
      if (j.eq.1) then ! copy
      do j1 = 1,a
       do j2 = 1,b
        s(j2,j1) = s(j2,j1) + u(j2,j1) * f
       enddo
      enddo
      else ! transpose
      do j1 = 1,a
       do j2 = 1,b
        s(j2,j1) = s(j2,j1) + u(j1,j2) * f
       enddo
      enddo
      endif
      return
      end
