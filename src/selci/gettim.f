      subroutine selci_gettim(cpud,elapsd)
      implicit double precision (a-h,o-z)
      real*4 etime,tt(2),elast,timediff
      external etime
      integer last,ienter
      save ienter,last,elast
      data ienter/0/
      if(ienter.eq.0) then
        ienter=1
        call cputm(last)
        elast = etime(tt)
      endif
      call cputm(now)
      elapsd = dfloat(now-last)*0.01d0
      timediff = etime(tt) - elast
      cpud = dble(timediff)
      end
