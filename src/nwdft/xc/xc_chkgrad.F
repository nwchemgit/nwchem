      logical function xc_chkgrad(xfac,cfac)
*
* $Id: xc_chkgrad.F,v 1.9 2000-07-21 20:40:42 edo Exp $
*
      implicit none
      double precision xfac(*) ! Exchange functional coeffs [input]
      double precision cfac(*) ! Correlat functional coeffs [input]
      
c
c     see if any Gradient Corrected Functional is not null
c
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(3) + xfac(4) + xfac(5) + xfac(6) + xfac(7) +
     + xfac(8) + xfac(9) + xfac(10) + xfac(11) +xfac(12)+ xfac(13)+
     + xfac(14)+
     +     cfac(2) + cfac(4) + cfac(5) +cfac(12)+cfac(13)+
     +     cfac(13)+cfac(14)+cfac(15) +cfac(16) + cfac(17)+cfac(18)+
     + cfac(19) + cfac(20)
      if (abs(tot).gt.eps)then
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      endif
      return
      end
