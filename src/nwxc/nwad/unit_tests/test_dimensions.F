C> \ingroup nwad_tests
C> @{
C>
C> \brief Test the dimensions of the arrays
C>
C> For every function value there are O(N) gradients, O(N**2) 2nd derivatives
C> and O(N**3) 3rd derivatives. However these quantities are not square arrays
C> but triangles. So we better check that we have the dimensions right. For
C> up to 2nd order that is trivial as we use that all the time, but for
C> 3rd order this is new.
C>
      program test_dimensions
      implicit none
      integer nd1, nd2, nd3, n, i1, i2, i3, nk1, nk2, nk3, num
      parameter (num = 10)
c
      write(*,'(2x,"n",8x,"Counted",7x,2x,7x,"Computed")')
      do n = 1, num
        nk1 = 0
        do i1 = 1, n
          nk1 = nk1 + 1
        enddo
        nd1 = n
        nk2 = 0
        do i1 = 1, n
          do i2 = 1, i1
            nk2 = nk2 + 1
          enddo
        enddo
        nd2 = n*(n+1)/2
        nk3 = 0
        do i1 = 1, n
          do i2 = 1, i1
            do i3 = 1, i2
              nk3 = nk3 + 1
            enddo
          enddo
        enddo
        nd3 = n*(n+1)*(n+2)/(2*3)
        write(*,'(x,i3,x,3i7,3x,3i7)')n,nk1,nk2,nk3,nd1,nd2,nd3
      enddo
      end
C>
C> @}
c $Id$
