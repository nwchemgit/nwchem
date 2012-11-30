c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Canonicalizes the integral labels
C>
C> This routine canonicalizes the integral labels of a 2-electron
C> integral \f$ (ij|kl) \f$ such that on output:
C> \f{eqnarray*}{
C>   i &\geq& j \\\\
C>   k &\geq& l \\\\
C>   ij = i(i-1)/2+j &\geq& kl = k(k-1)/2+l 
C> \f}
C>
c:tex-% part of the internal API routines
c:tex-\subsection{int\_canon}
c:tex-This routine canonicalizes integral lables such that:
c:tex-$i \geq j$, $k \geq l$, and $ij \geq kl$
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_canon(i, j, k, l, ii, jj, kk, ll)
c:tex-\end{verbatim}
      implicit none
c:tex-\begin{verbatim}
      integer i   !< [Input] i-label
      integer j   !< [Input] j-label
      integer k   !< [Input] k-label
      integer l   !< [Input] l-label
      integer ii  !< [Output] canonical i-label
      integer jj  !< [Output] canonical j-label
      integer kk  !< [Output] canonical k-label
      integer ll  !< [Output] canonical l-label
c:tex-\end{verbatim}
      integer iii,jjj
#if defined(DEBUG)
#include "stdio.fh"
      logical warn
#endif
c
c     return in ii,jj,kk,ll the cannonically ordered i,j,k,l
c
c     (i>=j) >= (k>=l)
c
#if defined(DEBUG)
      warn = .false.
#endif
      if (i .ge. j) then
        ii = i
        jj = j
      else
#if defined(DEBUG)
        warn = .true.
#endif
        ii = j
        jj = i
      endif
      if (k .ge. l) then
        kk = k
        ll = l
      else
#if defined(DEBUG)
        warn = .true.
#endif
        kk = l
        ll = k
      endif
      if (ii.lt.kk .or. (ii.eq.kk .and. jj.lt.ll)) then
#if defined(DEBUG)
        warn = .true.
#endif
        iii = ii
        jjj = jj
        ii = kk
        jj = ll
        kk = iii
        ll = jjj
      endif
#if defined(DEBUG)
      if (warn) then
        write(luout,10000)i,j,k,l
        write(luout,10001)ii,jj,kk,ll
      endif
10000 format(1x,'old:',4(i5))
10001 format(1x,'new:',4(i5))
#endif
c
      end
C> @}
