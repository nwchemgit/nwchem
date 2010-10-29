C $Id$
      Subroutine drdy_dsdt (n, g0, xlam, ts, f)
      Implicit none
      Integer n
      Double Precision g0(n),xlam(n),ts,f
C
      Integer i
      Double Precision t
C
C Evaluate ds/dt for Page-McIver gradient following algorithm
      f = 0.0d00
      do i = 1,n
         t = g0(i)*exp(-xlam(i)*ts)
         f = f + t*t
      enddo
      f = sqrt(f)
      return
      end
