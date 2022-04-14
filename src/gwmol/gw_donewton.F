      subroutine  gw_donewton(xin,xout,fx,idfx,xlow,xupp,fxlow,fxupp,
     $                        bracket,thresh,radius,iter)
C
C     Decide the step of the augement Newton algorithm (Figure 1) of
C     JCTC 17, 7504
C

      implicit none
      integer          :: iter
      logical          :: bracket
      double precision :: xin, xout, xlow, xupp
      double precision :: fx, idfx, fxlow, fxupp
      double precision :: thresh, radius

      logical          :: pole
      double precision :: z,step

      z = -1.0d0/idfx
      step = z*fx
      pole = .false.

      !
      ! when dfx is positive, we are sitting on a pole
      ! perhaps we should do somehting to get out of it
      !
      if (idfx.gt.-1d0) then
        pole = .true.
      endif  

      !
      !take full newton step for large enough z
      !
      if (z.gt.0.4d0 .and. z.le.1d0 ) then 
        xout = xin + step

      !
      ! when we have bracketed the solution, the
      ! golden section or regula falsi method might work best
      !
      elseif (bracket) then

        if (mod(iter,3).eq.0 .and. fxupp*fxlow.lt.0d0 ) then
          xout = (xlow*fxupp - xupp*fxlow)/(fxupp-fxlow)
        elseif (mod(iter,3).eq.1) then
          xout = xlow + 0.6180d0*(xupp-xlow)
        else
          xout = xupp - 0.6180d0*(xupp-xlow)
        endif

      !
      ! Take an scaled Newton step for small Z
      !
      elseif (z.gt.0.1d0) then
        xout = xin + 0.60d0*step

      !
      ! If Z is greater than 1.0 or very small we have to suspect
      ! numerical issues with the current derivative. Take a scaled
      ! steepest decent step 
      !
      else
        xout = xin + radius*fx

      endif


      !
      ! check that energy is inside bracket. If it is, then better take
      ! Golden section or Regular falsi method.
      !
      if (bracket) then
        if ( xout.gt.xupp .or. xout.lt.xlow ) then
          if (fxupp*fxlow.lt.0d0) then
            xout = (xlow*fxupp - xupp*fxlow)/(fxupp-fxlow)
          else
            xout = xlow + 0.6180d0*(xupp-xlow)
          endif
        endif
      endif

      end subroutine
