      subroutine selci_gettim(cpud,elapsd)
      implicit double precision (a-h,o-z)
      real*4 elast,timediff
      external util_cpusec
      integer last,ienter
      save ienter,last,elast
      data ienter/0/
      if(ienter.eq.0) then
        ienter=1
        call cputm(last)
        elast = util_cpusec()
      endif
      call cputm(now)
      elapsd = dfloat(now-last)*0.01d0
      timediff = util_cpusec() - elast
      cpud = dble(timediff)
      end
