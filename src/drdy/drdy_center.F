C $Id$
      Subroutine drdy_center (natom,n3,amass,x,xcm)
C
C Move center of mass of x to origin
C
C   x     - vector of mass-weighted coordinates
C   amass - vector of the sqrt root of the ratio of the atomic and reduced masses
C   xcm   - scratch array
C
      Implicit none 
      Integer natom,n3
      Double Precision x(n3),xcm(3),amass(n3)
*
      Integer i,ii,j
      Double Precision t,xmtot
C
C Note: x contains mass-weighted coordinates.  Need to first calculate
C  COM in plain cartesian and convert it to mass-scaled.  amass contains
C  sqrt(mi/redm).  On exit x array is still mass-scaled.
C
      do i = 1,3
         xcm(i) = 0.0d00
      enddo
      xmtot = 0.0d00
      do i = 1,natom
         ii = 3*(i-1)
         t = amass(ii+1)
         xmtot = xmtot + t*t
         do j = 1,3
            xcm(j) = xcm(j) + t*x(ii+j)
         enddo
      enddo
C xcm/xmtot is COM in plain cartesians
C Translate mass-scaled coordinates so new center of mass is at origin
      i = 0
      do j = 1,natom
         do ii = 1,3
            i = i + 1
            x(i) = x(i) - amass(i)*xcm(ii)/xmtot
         enddo
      enddo
      return
      end
