      subroutine  gw_donewton(xin,xout,fx,dfx,xlow,xupp,fxlow,fxupp,
     $                        bracket,thresh)
      implicit none
      logical          :: bracket
      double precision :: xin, xout, xlow, xupp
      double precision :: fx, dfx, fxlow, fxupp
      double precision :: thresh

      logical          :: pole
      double precision :: z

      z = -1d0/dfx
      pole = .false.

      ! when dfx is positive, we are sitting on a pole
      ! perhaps we should do somehting to get out of it
      if (dfx.gt.-1d0) then
        pole = .true.
      endif  


      !take full newton step for large enough z
      if (z.gt.0.3d0 .and. z.le.1d0 ) then 
        xout = xin - fx/dfx

      ! when we have bracketed the solution, the
      ! golden section or regula falsi method might work best
      elseif (bracket) then

        !use golden section for the others
        if (fx.eq.fxupp) then
          xout = xupp - 0.61803d0*(xupp-xlow)
        else
          xout = xlow + 0.61803d0*(xupp-xlow)
        endif

      ! take smaller step for intermediate z since
      ! the linear approximation might not work in
      ! this case
      elseif (z.gt.0.1d0 .and. z.le.0.3d0 ) then
        xout = xin - 0.61803d0*fx/dfx

      ! use a small newton step for the rest
      else
        xout = xin - 0.10d0*fx/dfx

      endif


      !check that energy is inside bracket
      if (bracket) then
        if ( xout.gt.xupp ) then
          xout = xupp - 1.1*thresh
        elseif (xout.lt.xlow ) then
          xout = xlow + 1.1*thresh
        endif
      endif

      end subroutine
