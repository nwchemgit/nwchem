      Subroutine grid_beckew(xyzw, xyz, Rij, rq, p, ictr, nq)
c
C$Id: grid_beckew.F,v 1.2 1999-04-29 16:36:54 d3h449 Exp $
c
      implicit none
c
#include "cdft.fh"
c
      integer nq      
      double precision xyzw(4,nq), xyz(3,*)
c
      double precision Rij(*)
      integer ictr
      double precision p(nq,*)
      double precision rq(nq,*)  ! sum of p(n) [ output]
c
      integer i, j, n, ij
      integer iind, jind, istep, iti
      double precision mu, nu
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision aaij, rrij
c
      double precision TOLL,toll_scr
      parameter(TOLL=0.95d0)
c     
c     Evaluate Becke space partitioning weight. Then, incorporate it
c     with weights from the single-center quadratures to form the total 
c     multi-center quadrature weight.
c     
      toll_scr=max((tol_rho*1d-5),1d-33)
      do 20 i = 1, ncenters
c
         x0i = xyz(1,i)
         y0i = xyz(2,i)
         z0i = xyz(3,i)
c
         do 10 n = 1,nq
c
            x = xyzw(1,n) - x0i
            y = xyzw(2,n) - y0i
            z = xyzw(3,n) - z0i
c
            rq(n,i) = sqrt(x*x + y*y + z*z) 
c
            p(n,i) = 1.D0
c
   10    continue
   20 continue
c
      do 50 i = 1, ncenters
c         iind = iandex(i)
         iind = i
         iti = iatype(iind)
         istep = (iind*(iind-1))/2
         do 40 j = 1, i
c
            if (i.eq.j) goto 40
c            jind = iandex(j)
            jind = j
c
            ij = istep + jind
            rrij = Rij(ij)
c          
            aaij = aratio(iti,iatype(jind))
c
            if (aaij.eq.-1234.d0)then
               write(*,*)' i, j, iind, jind, ij, aaij ',
     &                     i, j, iind, jind, ij, aaij
               call errquit('trouble in beckew; aratio .eq. -1234',911)
            endif
c
            if (rrij.eq.-1.0d0)then
               write(*,*)' i, j, iind, jind, ij, rrij ',
     &                     i, j, iind, jind, ij, rrij
               call errquit('trouble in beckew; rij .eq. = -1',911)
            endif
c          
            do 30 n = 1, nq
c
c               if (p(n,j).eq.0.d0.or.p(n,i).eq.0.d0) goto 30
               mu = (rq(n,i) - rq(n,j))*rrij
c
               if (mu.gt.TOLL) then
                  p(n,i)=0.d0
               elseif (mu.lt.-TOLL) then
                  p(n,j)=0.d0
               else
c
c                 i>j
c
c                 atomic size adjustment 
c
                  nu = mu + aaij*(1.D0 - mu*mu)
c
                  nu = 0.5D0*nu*(3.D0 - nu*nu) ! 1
                  nu = 0.5D0*nu*(3.D0 - nu*nu) ! 2
                  nu = 0.5D0*nu*(3.D0 - nu*nu) ! 3
c              
                  sk = 0.5D0 - 0.5D0*nu
c
                  p(n,i) = p(n,i)*sk
c
c                 flip i j
c
                  p(n,j) = p(n,j)*(1.d0-sk)
c
               endif
   30       continue
   40    continue
   50 continue
c
c     compute sum of Becke weights for normalization
c
      call dfill(nq*ncenters,0.d0,rq,1)
c      
      do i = 1, ncenters
         call daxpy(nq, 1.d0, p(1,i), 1, rq(1,1), 1)
      enddo
c
      do 90 n = 1, nq
         xyzw(4,n) = (p(n,ictr)*xyzw(4,n))/rq(n,1)
   90 continue
      return
      end
