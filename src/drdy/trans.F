C $Id$
      Subroutine drdy_trans (iop,n3,amass,x,dx)
C
      Implicit None
C
      Integer iop,n3
      Double Precision x(n3),dx(n3),amass(n3)
C
      Integer i
C
C Converts cartesian to mass-weighted coordinates (iop=1)
c   or mass-weighted to cartesian coordinates (iop=2)
c   also simultaneously converts derivatives
c   Note that amass is the atomic mass divided by the reduced mass.
c
      if (iop.eq.2) then
         do i = 1, n3
            dx(i) = dx(i)*amass(i)
            x(i) = x(i)/amass(i)
         enddo
      else
         do i = 1, n3
            dx(i) = dx(i)/amass(i)
            x(i) = x(i)*amass(i)
         enddo
      endif
      return
      end
