C $Id: trans.F,v 1.2 2001-01-08 15:52:01 windus Exp $
      Subroutine drdy_trans (iop,n3,amass,x,dx)
C
      Implicit None ! Double Precision (A-H,O-Z)
C
      Integer iop
      Integer n3
      Double Precision x(n3),dx(n3),amass(n3)
C
C Converts cartesian to mass-weighted coordinates (iop=1)
c   or mass-weighted to cartesian coordinates (iop=2)
c   also simultaneously converts derivatives
c   Note that amass is the atomic mass divided by the reduced mass.
c
*
      integer i
*
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
