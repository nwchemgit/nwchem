      Subroutine drdy_tranf (iop,ndim,n3,amass,f)
C
      Implicit Double Precision (A-H,O-Z)
C
      Dimension f(ndim,ndim),amass(n3)
C
C Converts force constant matrix from cartesians to mass-weighted (iop=1)
c   or mass-weighted to cartesian coordinates (iop=2)
c   Note that amass is the sqrt root of the ratio of the atomic and reduced masses.
c
      if (iop.eq.1) then
         do i = 1,n3
            do j = 1,n3
               f(j,i) = f(j,i)/(amass(i)*amass(j))
            enddo
         enddo
      else
         do i = 1,n3
            do j = 1,n3
               f(j,i) = f(j,i)*amass(i)*amass(j)
            enddo
         enddo
      endif
      return
      end
