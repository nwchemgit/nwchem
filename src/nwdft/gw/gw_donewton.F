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
      if (z.gt.0.4d0 .and. z.le.1d0 ) then 
        xout = xin - fx*idfx

      ! when we have bracketed the solution, the
      ! golden section or regula falsi method might work best
      elseif (bracket) then

        xout = (xlow*fxupp - xupp*fxlow)/(fxupp-fxlow)

      ! use a small newton step for the rest
      else
        xout = xin - 0.7d0*fx*idfx

      endif


      !check that energy is inside bracket
      if (bracket) then
        if ( xout.gt.xupp .or. xout.lt.xlow ) then
          xout = (xlow*fxupp - xupp*fxlow)/(fxupp-fxlow)
        endif
      endif

      end subroutine
