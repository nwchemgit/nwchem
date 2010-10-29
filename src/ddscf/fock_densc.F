      subroutine fock_density_screen(nfock,
     $     i1, j1, k1, l1, i2, j2, k2, l2,
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk, denmax)
      implicit none
c     
c     $Id$
c     
      integer nfock
      integer i1, j1, k1, l1, i2, j2, k2, l2
      integer ilo, jlo, klo, llo, ihi, jhi, khi, lhi
      double precision denmax
c     
      double precision dij(nfock, i1:i2, j1:j2)
      double precision dik(nfock, i1:i2, k1:k2)
      double precision dli(nfock, l1:l2, i1:i2)
      double precision djk(nfock, j1:j2, k1:k2)
      double precision dlj(nfock, l1:l2, j1:j2)
      double precision dlk(nfock, l1:l2, k1:k2)
c     
      integer ifock, i, j, k, l
      double precision sum
c     
c     OLD code was
c     Return the maximum value over the density matrices
c     of the sum of absolute values within the specified
c     range of each block
c
c     CURRENT code is
c     Return the max of sqrt(sum of squares) over the density in the specified 
c     blocks (which is rotationally invariant).
c     
c     Should put in better code for nfock != 1 (but it usually is)
c
      denmax = 0.0d0
      do ifock = 1, nfock          ! Assumption is nfock is usually one
         sum = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum = sum + dij(ifock,i,j)**2
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do i = ilo, ihi
               sum = sum + dik(ifock,i,k)**2
*               sum = sum + abs(dik(ifock,i,k))
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do i = ilo, ihi
            do l = llo, lhi
               sum = sum + dli(ifock,l,i)**2
*               sum = sum + abs(dli(ifock,l,i))
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do j = jlo, jhi
               sum = sum + djk(ifock,j,k)**2
*               sum = sum + abs(djk(ifock,j,k))
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do j = jlo, jhi
            do l = llo, lhi
               sum = sum + dlj(ifock,l,j)**2
*               sum = sum + abs(dlj(ifock,l,j))
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do l = llo, lhi
               sum = sum + dlk(ifock,l,k)**2
*               sum = sum + abs(dlk(ifock,l,k))
            enddo
         enddo
         sum = sqrt(sum)
         denmax = max(denmax, sum)
c     
      enddo
c     
      end
