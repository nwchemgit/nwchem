      logical function xc_chkgrad()
*
* $Id: xc_chkgrad.F,v 1.18 2002-04-05 21:37:45 sohirata Exp $
*
      implicit none
#include "cdft.fh"
      
c
c     see if any Gradient Corrected Functional is not null
c
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(3) + xfac(4) + xfac(5) + xfac(6) + xfac(7) +
     + xfac(8) + xfac(9) + xfac(10) + xfac(11) +xfac(12)+ xfac(13)+
     + xfac(14)+xfac(15)+xfac(16)+xfac(17)+
     +     cfac(2) + cfac(4) + cfac(5) +cfac(13)+
     +     cfac(13)+cfac(14)+cfac(15) +cfac(16) + cfac(17)+cfac(18)+
     + cfac(19) + cfac(20)+cfac(21)+cfac(22)+cfac(23)
      if (abs(tot).gt.eps)then
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      endif
      return
      end
      logical function xc_gotxc()
      implicit none
#include "cdft.fh"
      double precision tot
      integer i
c
c     return true if there is any XC functional to compute 
c     (other than HF Exch)
c
      tot = -abs(xfac(1))
      do i = 1, numfunc
         tot = tot + abs(xfac(i)) + abs(cfac(i))
      enddo
      xc_gotxc = tot.gt.1d-8
      return
      end
      integer function xc_getnqmx()
      implicit none
#include "cdft.fh"
      xc_getnqmx=nqmax
      return
      end
      integer function xc_getipol()
      implicit none
#include "cdft.fh"
      xc_getipol=ipol
      return
      end
      logical function xc_got2nd()
      implicit none
#include "cdft.fh"
c     return true for functional without 2nds
      xc_got2nd=.not.(
c     hcth functionals
     .     xfac(4).ne.0d0.or.
     .     xfac(5).ne.0d0.or.
     .     xfac(6).ne.0d0.or.
     .     xfac(10).ne.0d0.or.
     .     xfac(11).ne.0d0.or.
     .     xfac(12).ne.0d0.or.
     .     xfac(13).ne.0d0.or.
     .     xfac(14).ne.0d0.or.
     .     xfac(15).ne.0d0.or.
     .     cfac(13).ne.0d0.or.
     .     cfac(14).ne.0d0.or.
     .     cfac(15).ne.0d0.or.
     .     cfac(16).ne.0d0.or.
     .     cfac(17).ne.0d0.or.
     .     cfac(18).ne.0d0.or.
     .     cfac(19).ne.0d0.or.
     .     cfac(20).ne.0d0.or.
     .     cfac(21).ne.0d0.or.
c     cperdew91
     .     cfac(5).ne.0d0)

      return
      end

