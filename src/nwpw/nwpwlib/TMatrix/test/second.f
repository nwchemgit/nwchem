      SUBROUTINE current_second(T)
      real*8 T
      real*4 dummy,etime
      real*4 s(2)
      real   t2(2)

      dummy = etime(s)
      T = dble(s(1))

      RETURN
      END

