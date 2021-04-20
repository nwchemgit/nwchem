      subroutine  gw_donewton(xin,xout,fx,idfx,xlow,xupp,fxlow,fxupp,
     $                        bracket,thresh)
      implicit none
      logical          :: bracket
      double precision :: xin, xout, xlow, xupp
      double precision :: fx, idfx, fxlow, fxupp
      double precision :: thresh

      logical          :: pole
      double precision :: z

      z = -idfx
      pole = .false.

      ! when dfx is positive, we are sitting on a pole
      ! perhaps we should do somehting to get out of it
      if (idfx.gt.-1d0) then
        pole = .true.
      endif  


      !take full newton step for large enough z
      if (z.gt.0.2d0 .and. z.le.1d0 ) then 
        xout = xin - fx*idfx

      ! when we have bracketed the solution, the
      ! golden section or regula falsi method might work best
      elseif (bracket) then

        !use golden section for the others
        if (fx.eq.fxupp) then
          xout = xupp - 0.61803d0*(xupp-xlow)
        else
          xout = xlow + 0.61803d0*(xupp-xlow)
        endif

      ! use a small newton step for the rest
      elseif (z.gt.0.05d0) then
        xout = xin - 1.0d0*fx*idfx

      else
        xout = xin - 1.0d0*fx*idfx

      endif


      !check that energy is inside bracket
      if (bracket) then
        if ( xout.gt.xupp ) then
          xout = xupp - 1.1d0*thresh
        elseif (xout.lt.xlow ) then
          xout = xlow + 1.1d0*thresh
        endif
      endif

      end subroutine
