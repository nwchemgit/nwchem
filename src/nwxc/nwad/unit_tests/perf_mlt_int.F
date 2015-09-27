      program perf_mlt_int
      implicit none
      integer ii, jj
      integer nn
      parameter(nn = 1000000000)
      double precision aa
      double precision rr
      aa = 1.0d0
      ii = 2
      do jj = 1, nn
        rr = aa*ii
        aa = aa+1.0d0
      enddo
      write(*,*)rr
      end
