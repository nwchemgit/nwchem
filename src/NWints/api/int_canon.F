c $Id: int_canon.F,v 1.1 1994-05-05 19:07:41 d3e129 Exp $
      subroutine int_canon(i, j, k, l, ii, jj, kk, ll)
      implicit none
      integer i,j,k,l,ii,jj,kk,ll,iii,jjj
#if defined(DEBUG)
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
        write(6,10000)i,j,k,l
        write(6,10001)ii,jj,kk,ll
      endif
10000 format(1x,'old:',4(i5))
10001 format(1x,'new:',4(i5))
#endif
c
      end
