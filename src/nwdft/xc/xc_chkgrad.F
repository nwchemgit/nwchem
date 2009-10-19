      logical function xc_chkgrad()
*
* $Id: xc_chkgrad.F,v 1.38 2008-09-29 20:25:06 edo Exp $
*
      implicit none
#include "cdft.fh"
#include "zora.fh"
#include "steric.fh"
c
c     see if any Gradient Corrected Functional is not null
c
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(3) + xfac(4) + xfac(5) + xfac(6) + xfac(7) +
     + xfac(8) + xfac(9) + xfac(10) + xfac(11) +xfac(12)+ xfac(13)+
     + xfac(14)+xfac(15)+xfac(16)+xfac(17)+xfac(18) +xfac(19) +
     +     xfac(22)+
     +     xfac(30)+xfac(31)+ xfac(33) + xfac(39) + cfac(29) +
     +     cfac(2) + cfac(4) + cfac(5) + cfac(12) +
     +     cfac(13)+cfac(14)+cfac(15) +cfac(16) + cfac(17)+cfac(18)+
     + cfac(19) + cfac(20)+cfac(21)+cfac(22)+cfac(23) + cfac(24) +
     + xfac(25) + cfac(25) +cfac(26)+xfac(20)+xfac(21)+cfac(27) +
     + xfac(26) + xfac(27) + xfac(28) + xfac(29) + cfac(31) + 
     + cfac(32) + cfac(33) + cfac(34) + cfac(35) + cfac(36) +
     + xfac(34) + xfac(35) + xfac(36) + xfac(37) +
     + cfac(30) + cfac(38) + cfac(39) + cfac(40) + cfac(41) +
     + xfac(46) + cfac(46) + xfac(47) + cfac(47) + xfac(40) + xfac(41)

c    +     cfac(28) + cfac(36)
      if (abs(tot).gt.eps.or.lb94.or.cs00)then
         xc_chkgrad = .true.
      else if (do_zora) then    ! turn on gradients for zora
         xc_chkgrad = .true.
      else if (lsteric) then ! turn on gradients for steric calculations
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      end if

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
      xc_gotxc = tot.gt.1d-8.or.lb94.or.cs00
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
     .     xfac(17).ne.0d0.or.
     .     xfac(22).ne.0d0.or.
     .     cfac(13).ne.0d0.or.
     .     cfac(14).ne.0d0.or.
     .     cfac(15).ne.0d0.or.
     .     cfac(16).ne.0d0.or.
     .     cfac(17).ne.0d0.or.
     .     cfac(18).ne.0d0.or.
     .     cfac(19).ne.0d0.or.
     .     cfac(20).ne.0d0.or.
     .     cfac(21).ne.0d0.or.
     .     cfac(22).ne.0d0.or.
     .     cfac(23).ne.0d0.or.
     .     cfac(28).ne.0d0.or.
     .     cfac(36).ne.0d0)
      return
      end

