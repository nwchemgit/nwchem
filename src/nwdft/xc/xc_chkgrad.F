      logical function xc_chkgrad(xfac,cfac)
      implicit none
      double precision xfac(*) ! Exchange functional coeffs [input]
      double precision cfac(*) ! Correlat functional coeffs [input]
      
c
c     see if any Gradient Corrected Functional is not null
c
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(3) + xfac(4) + xfac(5) + xfac(6) + xfac(7) +
     + xfac(8) + xfac(9) +
     +     cfac(2) + cfac(4) + cfac(5) +cfac(12)
      if (abs(tot).gt.eps)then
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      endif
      return
      end
