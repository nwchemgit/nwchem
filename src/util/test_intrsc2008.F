      program test_intrsc2008
      implicit none
      integer i,j,k
      integer util_popcnt
      do i=-1000000,1000000
        j = util_popcnt(i)
        k = popcnt(i)
        if (j.ne.k) then
            print*,'1 error:',i,j,k
        endif
      enddo
      do i=-63,63
        j = util_popcnt(2**i)
        k = popcnt(2**i)
        if (j.ne.k) then
            print*,'2 error:',i,j,k
        endif
      enddo
      end program test_intrsc2008
