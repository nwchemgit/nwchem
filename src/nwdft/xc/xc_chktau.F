      logical function xc_chktau()
C
C $Id: xc_chktau.F,v 1.5 2006-02-24 21:50:53 edo Exp $ 
C
      implicit none
 
c
c     see if any metaGGA Functional using TAU is not null
c
 
 
#include "cdft.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
 
 
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25)+xfac(21)+cfac(27) +
     +     xfac(28) + xfac(29) +
     +     cfac(31) + cfac(32) + cfac(33) + cfac(34) + cfac(35)
      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      return
      end
      subroutine xc_getxcfac(rtdb)
      implicit none
      integer rtdb
c     load cfac and xfac into cdft.fh
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,
     &   xfac))call errquit(' xcgetxcfac ',0,0)
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc,
     &   cfac))call errquit(' xcgetxcfac ',0,0)
      return
      end



