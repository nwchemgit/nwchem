      program perf_mlt_dbl
      implicit none
      integer jj
      integer nn
      parameter(nn = 1000000000)
      double precision aa
      double precision bb
      double precision rr
      aa = 1.0d0
      bb = 2.0d0
      do jj = 1, nn
        rr = aa*bb
        aa = aa+1.0d0
      enddo
      write(*,*)rr
      end
