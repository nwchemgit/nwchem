*      call anl_fit_init
*      end
      double precision function anl_fit(n, l, x)
*
* $Id$
*
      implicit none
#include "canlfit.fh"
      integer n, l
      double precision x
c
c     Return the value of Anl(x) computed from the table
c     or asymptotic formula.  See anl_fit_init.
c
      integer nlcnt, i, k
      double precision d, value, rx, rx2
c
      if (x .lt. anl_xcut(n)) then
         nlcnt = anl_mapnl(l,n)
         i = int(anl_scale*x)
         d = x - anl_table(0,i,nlcnt)
         value = anl_table(anl_order+1,i,nlcnt)
         do k = anl_order,1,-1
            value = anl_table(k,i,nlcnt) + d*value
         enddo
         anl_fit = value
      else
         rx = 1.0d0/x
         value = anl_fac(l,n)*rx
         rx2 = rx*rx
         do i = 1, l
            value = value * rx2
         enddo
         anl_fit = value        ! anl_fac(l,n)/x**(2*l+1)
      endif
c
      end
      subroutine anl_fit_group(nmax, x, values)
      implicit none
#include "canlfit.fh"
      integer nmax
      double precision x
      double precision values(*)
c
c     Return the value of Anl(x) computed from the table
c     or asymptotic formula for all n and l up to nmax.
c
c     No error checking is being performed for reasons of speed.
c
c     Several values are returned, as follows.
c
c     ind = 1
c     do n = 0, nmax
c     .  do l = n, 0, -2
c     .     values(ind) = anl_fit(n, l, x)
c     .     ind = ind + 1
c     .  enddo
c     end
c
      integer i, k, n, l, ind
      double precision d, d2, d3, value, rx, rx2
      double precision rx2l1(0:anl_maxn)
c     
      if (x .lt. anl_xcut(nmax)) then
         i = int(anl_scale*x)
         d = x - anl_table(0,i,1)
         if (anl_order .eq. 3) then
            d2 = d*d
            d3 = d*d2
            do ind = 1, anl_nlsum(nmax)
               values(ind) = anl_table(1,i,ind) + d*anl_table(2,i,ind) +
     $              d2*anl_table(3,i,ind) + d3*anl_table(4,i,ind)
            enddo
         else
            do ind = 1, anl_nlsum(nmax)
               value = anl_table(anl_order+1,i,ind)
               do k = anl_order,1,-1
                  value = anl_table(k,i,ind) + d*value
               enddo
               values(ind) = value
            enddo
         endif
      else
         rx = 1.0d0/x
         rx2 = rx*rx
         rx2l1(0) = rx           ! x^-(2l+1)
         do l = 1, nmax
            rx2l1(l) = rx2l1(l-1)*rx2
         enddo
c
         ind = 1
         do n = 0, nmax
            do l = n, 0, -2
               values(ind) = anl_fac(l,n)*rx2l1(l)
               ind = ind + 1
            enddo
         enddo
      endif
c
      end
      subroutine anl_fit_group_0(nmax, x, values)
      implicit none
#include "canlfit.fh"
      integer nmax
      double precision x
      double precision values(*)
c
c     Return the value of Anl(x) computed from the table
c     or asymptotic formula for all n and l up to nmax.
c
c     HARDWIRED FOR NMAX = 0
c
c     No error checking is being performed for reaons of speed.
c
c     Several values are returned, as follows.
c
c     ind = 1
c     do n = 0, nmax
c     .  do l = n, 0, -2
c     .     values(ind) = anl_fit(n, l, x)
c     .     ind = ind + 1
c     .  enddo
c     end
c
      integer i, k
      double precision d, d2, d3
c
      if (x .lt. anl_xcut(0)) then
         i = int(anl_scale*x)
         d = x - anl_table(0,i,1)
         if (anl_order .eq. 3) then
            d2 = d*d
            d3 = d2*d
            values(1) = anl_table(1,i,1) + d*anl_table(2,i,1) +
     $           d2*anl_table(3,i,1) + d3*anl_table(4,i,1)
         else
            values(1) = anl_table(anl_order+1,i,1)
            do k = anl_order,1,-1
               values(1) = anl_table(k,i,1) + d*values(1)
            enddo
         endif
      else
         values(1) = anl_fac(0,0)/x
      endif
c
      end
      subroutine anl_fit_group_1(nmax, x, values)
      implicit none
#include "canlfit.fh"
      integer nmax
      double precision x
      double precision values(*)
c
c     Return the value of Anl(x) computed from the table
c     or asymptotic formula for all n and l up to nmax.
c
c     HARDWIRED FOR NMAX = 1
c
c     No error checking is being performed for reaons of speed.
c
c     Several values are returned, as follows.
c
c     ind = 1
c     do n = 0, nmax
c     .  do l = n, 0, -2
c     .     values(ind) = anl_fit(n, l, x)
c     .     ind = ind + 1
c     .  enddo
c     end
c
      integer i, k
      double precision d, d2, d3, rx, rx2, rx3
c
      if (x .lt. anl_xcut(1)) then
         i = int(anl_scale*x)
         d = x - anl_table(0,i,1)
         if (anl_order .eq. 3) then
            d2 = d*d
            d3 = d2*d
            values(1) = anl_table(1,i,1) + d*anl_table(2,i,1) +
     $           d2*anl_table(3,i,1) + d3*anl_table(4,i,1)
            values(2) = anl_table(1,i,2) + d*anl_table(2,i,2) +
     $           d2*anl_table(3,i,2) + d3*anl_table(4,i,2)
         else
            values(1) = anl_table(anl_order+1,i,1)
            values(2) = anl_table(anl_order+1,i,2)
            do k = anl_order,1,-1
               values(1) = anl_table(k,i,1) + d*values(1)
               values(2) = anl_table(k,i,2) + d*values(2)
            enddo
         endif
      else
         rx = 1.0d0/x
         rx2 = rx*rx
         rx3 = rx2*rx
         values(1) = anl_fac(0,0)*rx
         values(2) = anl_fac(1,1)*rx3
      endif
c
      end
      subroutine anl_fit_group_2(nmax, x, values)
      implicit none
#include "canlfit.fh"
      integer nmax
      double precision x
      double precision values(*)
c
c     Return the value of Anl(x) computed from the table
c     or asymptotic formula for all n and l up to nmax.
c
c     HARDWIRED FOR NMAX = 2
c
c     No error checking is being performed for reaons of speed.
c
c     Several values are returned, as follows.
c
c     ind = 1
c     do n = 0, nmax
c     .  do l = n, 0, -2
c     .     values(ind) = anl_fit(n, l, x)
c     .     ind = ind + 1
c     .  enddo
c     end
c
      integer i, k, ind
      double precision d, d2, d3, rx, rx2, rx3
c
      if (x .lt. anl_xcut(2)) then
         ind = 1
         i = int(anl_scale*x)
         d = x - anl_table(0,i,1)
         if (anl_order .eq. 3) then
            d2 = d*d
            d3 = d2*d
            values(1) = anl_table(1,i,1) + d*anl_table(2,i,1) +
     $           d2*anl_table(3,i,1) + d3*anl_table(4,i,1)
            values(2) = anl_table(1,i,2) + d*anl_table(2,i,2) +
     $           d2*anl_table(3,i,2) + d3*anl_table(4,i,2)
            values(3) = anl_table(1,i,3) + d*anl_table(2,i,3) +
     $           d2*anl_table(3,i,3) + d3*anl_table(4,i,3)
            values(4) = anl_table(1,i,4) + d*anl_table(2,i,4) +
     $           d2*anl_table(3,i,4) + d3*anl_table(4,i,4)
         else
            values(1) = anl_table(anl_order+1,i,1)
            values(2) = anl_table(anl_order+1,i,2)
            values(3) = anl_table(anl_order+1,i,3)
            values(4) = anl_table(anl_order+1,i,4)
            do k = anl_order,1,-1
               values(1) = anl_table(k,i,1) + d*values(1)
               values(2) = anl_table(k,i,2) + d*values(2)
               values(3) = anl_table(k,i,3) + d*values(3)
               values(4) = anl_table(k,i,4) + d*values(4)
            enddo
         endif
      else
         rx = 1.0d0/x
         rx2 = rx*rx
         rx3 = rx2*rx
         values(1) = anl_fac(0,0)*rx
         values(2) = anl_fac(1,1)*rx3
         values(3) = anl_fac(2,2)*rx3*rx2
         values(4) = anl_fac(0,2)*rx
      endif
c
      end
      subroutine anl_fit_init
      implicit none
#include "errquit.fh"
#include "canlfit.fh"
c
      integer n, l, i, nlcnt, k
      double precision f(0:anl_npt), h, x, value, anl, pi, dfac, 
     $     anl_fit, maxerr, maxrelerr, test
      logical initialized
      external anl, anl_fit
      data initialized /.false./

      if (initialized) return
      initialized = .true.
c
c     Initialize interpolation tables used for 0..xhi
c
      h = anl_xhi/dble(anl_npt)
      anl_scale = 1.0d0/h
c
      do n = 0, anl_maxn
         anl_xcut(n) = anl_xhi          ! Initial value - reset below
         do l = 0, anl_maxn
            anl_mapnl(l,n) = 999999999 ! Make invalid access disastrous
            anl_fac(l,n) = 1d+300
         enddo
      enddo
c
      maxerr = 0.0d0
      maxrelerr = 0.0d0
      nlcnt = 0
      do n = 0,anl_maxn
         do l = n,0,-2
            nlcnt = nlcnt + 1
            if (nlcnt .gt. anl_maxnl) call errquit
     $           (' anlfit: too many functions?', nlcnt, INPUT_ERR)
            anl_mapnl(l,n) = nlcnt
            do i = 0,anl_npt
               f(i) = anl(n,l,dble(i)*h)
            enddo
            call interp_table_make(anl_npt, anl_order, 0.0d0, 
     $           anl_xhi, f, 
     $           anl_table(0,0,anl_mapnl(l,n)))
c
c     Verify accuracy of interpolation at midpoints of regions
c     ... make sure we have either ANL_ACC sig. fig. or an absolute
c     accuracy of ANL_ACC
c
            do i = 0,anl_npt-1
               x = h*(0.5d0+dble(i)) 
c$$$               d = x - anl_table(0,i,nlcnt)
c$$$               value = anl_table(anl_order+1,i,nlcnt)
c$$$               do k = anl_order,1,-1
c$$$                  value = anl_table(k,i,nlcnt) + d*value
c$$$               enddo
               value = anl_fit(n,l,x)
               test = anl(n,l,x)
               maxerr = max(maxerr,abs(value-test))
               maxrelerr = max(maxrelerr,abs(value-test)/test)
               if (abs(value-test) .gt. 
     $              max(anl_acc,anl_acc*value)) then
                  write(6,1) n, l, x, value, anl(n,l,x),
     $                 value-anl(n,l,x)
 1                format(' fit error ', 2i5, 3f20.14, 1p, d9.1)
                  call errquit('anlfit: bad iterpolation',i, INPUT_ERR)
               endif
            enddo
         enddo
         anl_nlsum(n) = nlcnt
      enddo
c      write(6,*) ' MAXERR ', maxerr
c      write(6,*) ' MAXRELERR ', maxrelerr
**      write(6,*) ' used ', nlcnt, ' tables out of ', anl_maxnl
c      
c     Initialize prefactors for asymptotic form used beyond xhi
c
c     (n+l+1)!! * sqrt(Pi) / [2^((n+l)/2+2) * (2*l+1) * x^(2*l+1)]
c
      pi = 4.0d0*atan(1.0d0)
      do n = 0, anl_maxn
         do l = n, 0, -2
            dfac = 1.0d0
            do k = (n+l+1),1,-2
               dfac = dfac * dble(k)
            enddo
            anl_fac(l,n) = dfac * sqrt(pi) / 
     $           ((2.0d0*dble(l)+1.0d0)*2.0d0**((n+l)/2+2))
c
c     Verify asymptotic form at xhi
c
            value = anl(n,l,anl_xhi)
            if (abs(value-anl_fit(n,l,anl_xhi)).gt.anl_acc) then
               write(6,2) n, l, dfac, value, anl_fac(l,n), 
     $              anl_fit(n,l,anl_xhi), 
     $              anl_fit(n,l,anl_xhi)-value
 2             format(' prefac err ', 2i5, 4f20.14,1p,d12.4)
               call errquit('anlfit: bad asymptote ', n*1000+l,
     &       INPUT_ERR)
            endif
         enddo
c
c     Find where asymptotic form becomes accurate.  Worst case
c     is always l=0/2.
c      
         l = mod(n,2)
         do x = 3.0d0, anl_xhi, 0.1d0
            if (abs(anl_fac(l,n)/x**(2*l+1)-anl(n,l,x)).le.anl_acc)then
               anl_xcut(n) = x
               goto 333
            endif
         enddo
         call errquit('aln_fit_init: failed to find asymptote?',n,
     &       UNKNOWN_ERR)
 333     continue
      enddo
c
      end
      double precision function anl(n,l,x)
      implicit none
      integer n, l
      double precision x
c
c     Compute the radial part of the potential due to a
c     general gaussian function times a real, solid, spherical harmonic
c
c     A(x) = 0.5*x**(n-l) (n+l+1)F((n+l)/2,x**2) + (n-l)*x**l*I(n-l-1)/2
c
      integer k
      double precision ik, xkm1, xsq, fnl
      double precision fm_slow_but_accurate
      external fm_slow_but_accurate

      if (mod(n-l+1,2) .ne. 1) stop 543
      if (n.lt.l) then 
         stop 544
      else if (n .eq. l) then
         ik = 0.0d0
      else
         ik = 1.0d0
         xsq = x*x
         xkm1 = xsq
         do k = 3, n-l-1, 2
            ik = 0.5d0*dble(k-1)*ik + xkm1
            xkm1 = xkm1 * xsq
         enddo
         ik = ik * exp(-x*x) * 0.5d0
*     write(6,*) ' n-l-1, x, ik ', n-l-1, x, ik
      endif
      
      fnl = fm_slow_but_accurate((n+l)/2,x*x)
      if (n.eq.l) then
        anl = (0.5d0*dble(n+l+1)*fnl) / dble(l+l+1)
      else
        anl = (0.5d0*(x**(n-l))*dble(n+l+1)*fnl + 
     $         0.5d0*dble(n-l)*ik) / dble(l+l+1)
      endif
c
**      call fmvector(1,1,x*x,fm)
**      anl = 0.5d0*3.0d0*fm(2)

      end
      double precision function fm_slow_but_accurate(m,t)
      implicit none
      integer m
      double precision t
c
c     Return Fm(t) computed by downward recursion so
c     that it has full precision.
c
      double precision et, twot, fp
      integer j
c
      if (t .gt. 690.0d0) then
         et = 0.0d0
      else
         et = exp(-t)
      endif
      twot = 2.0d0 * t
c
      fp = 0.0d0
      do j = 200, m, -1
         fp = (twot*fp + et) / dble(j+j+1)
      enddo
c
      fm_slow_but_accurate = fp
c
      end
      subroutine interp_table_make(n, order, xlo, xhi, f, table)
      implicit none
      integer n, order
      double precision xlo, xhi, table(0:order+1,0:n-1), f(0:n)
c
c     Given a tabulation of f at n+1 evenly spaced points in [xlo:xhi]
c     return in table+base+scale the info necessary to compute a piecewise 
c     interpolating polynomial approximation.
c
c     The approximation may be computed from the table at a point
c     xlo <= x <= xhi by the following
c
c     scale = dble(n)/(xhi-xlo)
c     k = int((x-base)*scale)
c     d = x-table(0,k)
c     f = table(order+1,k)
c     do i = order,1,-1
c     .   f = table(i,k) + d*f
c     enddo
c
      integer i, j, jlo, jhi
      integer ordermax
      parameter (ordermax = 10)
      double precision x(ordermax+1), c(ordermax+1), y(ordermax+1),
     $     work(2*ordermax+2), h, xx
c
      h = (xhi - xlo) / dble(n)
c
*      do i = 0, n
*         write(6,*) ' f ', i, f(i)
*      enddo
c
      do i = 0,n-1
         xx = xlo + dble(i)*h ! Reference point for interpolation
c
c     Figure out which points we're going to interpolate
c
         jlo = i - order/2
         jhi = jlo + order
         if (jlo .lt. 0) then
            jlo = 0
            jhi = jlo + order
         else if (jhi .gt. n) then
            jhi = n
            jlo = jhi - order
         endif
         do j = jlo, jhi
            x(j-jlo+1) = (xlo + dble(j)*h) - xx
            y(j-jlo+1) = f(j)
         enddo
c
         call dplint(order+1,x,y,c)
         call dpolcf(0.0d0,order+1,x,c,table(1,i),work)
         table(0,i) = xx
c
*         write(6,*) ' Table  ', i, jlo, jhi, xx
*         write(6,1) (table(j,i),j=0,order+1)
* 1       format(f10.4,2x,10(f12.6))
c
      enddo
c
      end
c$$$      implicit double precision (a-h, o-z)
c$$$      double precision x(4), c(4), y(4), d(4), work(8)
c$$$      xx = 0.0d0
c$$$      n = 4
c$$$      x(1) = 0
c$$$      x(2) = 1
c$$$      x(3) = 2
c$$$      x(4) = 3
c$$$      
c$$$      y(1) = 0
c$$$      y(2) = 3    
c$$$      y(3) = 1    
c$$$      y(4) = 3    
c$$$
c$$$      call dplint(n,x,y,c)
c$$$      call dpolcf(xx,n,x,c,d,work)
c$$$
c$$$      write(6,*) d
c$$$      end
c$$$
      SUBROUTINE DPOLCF (XX, N, X, C, D, WORK)
C***BEGIN PROLOGUE  DPOLCF
C***PURPOSE  Compute the coefficients of the polynomial fit (including
C            Hermite polynomial fits) produced by a previous call to
C            POLINT.
C***LIBRARY   SLATEC
C***CATEGORY  E1B
C***TYPE      DOUBLE PRECISION (POLCOF-S, DPOLCF-D)
C***KEYWORDS  COEFFICIENTS, POLYNOMIAL
C***AUTHOR  Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C        Subroutine DPOLCF computes the coefficients of the polynomial
C     fit (including Hermite polynomial fits ) produced by a previous
C     call to DPLINT.  The coefficients of the polynomial, expanded
C     about XX, are stored in the array D. The expansion is of the form
C     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... +
C                                                  D(N)*((Z-XX)**(N-1)).
C     Between the call to DPLINT and the call to DPOLCF the variable N
C     and the arrays X and C must not be altered.
C
C     *****  INPUT PARAMETERS
C      *** All TYPE REAL variables are DOUBLE PRECISION ***
C
C     XX   - The point about which the Taylor expansion is to be made.
C
C     N    - ****
C            *     N, X, and C must remain unchanged between the
C     X    - *     call to DPLINT and the call to DPOLCF.
C     C    - ****
C
C     *****  OUTPUT PARAMETER
C      *** All TYPE REAL variables are DOUBLE PRECISION ***
C
C     D    - The array of coefficients for the Taylor expansion as
C            explained in the abstract
C
C     *****  STORAGE PARAMETER
C
C     WORK - This is an array to provide internal working storage. It
C            must be dimensioned by at least 2*N in the calling program.
C
C
C     **** Note - There are two methods for evaluating the fit produced
C     by DPLINT. You may call DPOLVL to perform the task, or you may
C     call DPOLCF to obtain the coefficients of the Taylor expansion and
C     then write your own evaluation scheme. Due to the inherent errors
C     in the computations of the Taylor expansion from the Newton
C     coefficients produced by DPLINT, much more accuracy may be
C     expected by calling DPOLVL as opposed to writing your own scheme.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890213  DATE WRITTEN
C   891006  Cosmetic changes to prologue.  (WRB)
C   891024  Corrected KEYWORD section.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DPOLCF
C
      INTEGER I,IM1,K,KM1,KM1PI,KM2N,KM2NPI,N,NM1,NMKP1,NPKM1
      DOUBLE PRECISION C(*),D(*),PONE,PTWO,X(*),XX,WORK(*)
C***FIRST EXECUTABLE STATEMENT  DPOLCF
      DO 10010 K=1,N
      D(K)=C(K)
10010 CONTINUE
      IF (N.EQ.1) RETURN
      WORK(1)=1.0D0
      PONE=C(1)
      NM1=N-1
      DO 10020 K=2,N
      KM1=K-1
      NPKM1=N+K-1
      WORK(NPKM1)=XX-X(KM1)
      WORK(K)=WORK(NPKM1)*WORK(KM1)
      PTWO=PONE+WORK(K)*C(K)
      PONE=PTWO
10020 CONTINUE
      D(1)=PTWO
      IF (N.EQ.2) RETURN
      DO 10030 K=2,NM1
      KM1=K-1
      KM2N=K-2+N
      NMKP1=N-K+1
      DO 10030 I=2,NMKP1
      KM2NPI=KM2N+I
      IM1=I-1
      KM1PI=KM1+I
      WORK(I)=WORK(KM2NPI)*WORK(IM1)+WORK(I)
      D(K)=D(K)+WORK(I)*D(KM1PI)
10030 CONTINUE
      RETURN
      END
      SUBROUTINE DPLINT (N, X, Y, C)
C***BEGIN PROLOGUE  DPLINT
C***PURPOSE  Produce the polynomial which interpolates a set of discrete
C            data points.
C***LIBRARY   SLATEC
C***CATEGORY  E1B
C***TYPE      DOUBLE PRECISION (POLINT-S, DPLINT-D)
C***KEYWORDS  POLYNOMIAL INTERPOLATION
C***AUTHOR  Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C        Subroutine DPLINT is designed to produce the polynomial which
C     interpolates the data  (X(I),Y(I)), I=1,...,N.  DPLINT sets up
C     information in the array C which can be used by subroutine DPOLVL
C     to evaluate the polynomial and its derivatives and by subroutine
C     DPOLCF to produce the coefficients.
C
C     Formal Parameters
C     *** All TYPE REAL variables are DOUBLE PRECISION ***
C     N  - the number of data points  (N .GE. 1)
C     X  - the array of abscissas (all of which must be distinct)
C     Y  - the array of ordinates
C     C  - an array of information used by subroutines
C     *******  Dimensioning Information  *******
C     Arrays X,Y, and C must be dimensioned at least N in the calling
C     program.
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPLINT
      INTEGER I,K,KM1,N
      DOUBLE PRECISION DIF,C(*),X(*),Y(*)
C***FIRST EXECUTABLE STATEMENT  DPLINT
      IF (N .LE. 0) GO TO 91
      C(1)=Y(1)
      IF(N .EQ. 1) RETURN
      DO 10010 K=2,N
      C(K)=Y(K)
      KM1=K-1
      DO 10010 I=1,KM1
C     CHECK FOR DISTINCT X VALUES
      DIF = X(I)-X(K)
      IF (DIF .EQ. 0.0) GO TO 92
      C(K) = (C(I)-C(K))/DIF
10010 CONTINUE
      RETURN
   91 write(6,*) 'SLATEC', 'DPLINT', 'N IS ZERO OR NEGATIVE.', 2, 1
      RETURN
   92 write(6,*) 'SLATEC', 'DPLINT',
     +   'THE ABSCISSAS ARE NOT DISTINCT.', 2, 1
      RETURN
      END

