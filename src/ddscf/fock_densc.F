      subroutine fock_density_screen(nfock,
     $     i1, j1, k1, l1, i2, j2, k2, l2,
     $     ilo, jlo, klo, llo, ihi, jhi, khi, lhi,
     $     dij, dik, dli, djk, dlj, dlk, denmax)
      implicit none
c     
c     $Id: fock_densc.F,v 1.2 1996-06-27 22:42:34 d3g681 Exp $
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
c     Return the maximum value over the density matrices
c     of the sum of absolute values within the specified
c     range of each block
c     
c     Should put in better code for nfock != 1 (but it usually is)
c
      denmax = 0.0d0
      do ifock = 1, nfock          ! Assumption is nfock is usually one
         sum = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum = sum + abs(dij(ifock,i,j))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do i = ilo, ihi
               sum = sum + abs(dik(ifock,i,k))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do i = ilo, ihi
            do l = llo, lhi
               sum = sum + abs(dli(ifock,l,i))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do j = jlo, jhi
               sum = sum + abs(djk(ifock,j,k))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do j = jlo, jhi
            do l = llo, lhi
               sum = sum + abs(dlj(ifock,l,j))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
         sum = 0.0d0
         do k = klo, khi
            do l = llo, lhi
               sum = sum + abs(dlk(ifock,l,k))
            enddo
         enddo
         denmax = max(denmax, sum)
c     
      enddo
c     
      end
