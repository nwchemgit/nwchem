      subroutine dimqm_fldpot(d_xyz, q_xyz, q_chrg, fld)
c      use err_func
      implicit none
#include "stdio.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c      Input Variables
      double precision d_xyz(3)
      double precision q_xyz(3)
      double precision q_chrg
      double precision fld(3)
c      Local Variables
      double precision r(3)
      double precision dist, dist3
      double precision screen
      double precision screen_fld
c
c      Distance between points
c     
c      if(ldebug) write(luout,*) "Start FldPot"
      r(:) = d_xyz(:) - q_xyz(:)
      dist = SQRT(DOT_PRODUCT(r, r))
c      write(luout,*) "Distance:",dist
c      write(luout,*) "Charge:",q_chrg
c      write(luout,*) "R:",r
c
c      Determine screening factor based off screening type
c      
      screen_fld = ONE
      if(scrnType == ERFSCRN) then
        screen     = erf(scrnFactor*dist)
        screen_fld = screen * screen
      else if(scrnType == EXPSCRN) then
        screen = ONE - EXP(-scrnFactor*dist)
        screen_fld = screen * screen * screen
      else
        if(dist < 1.0d-12) dist = 1.0d-12
      end if
c
      dist3 = dist * dist * dist
c      write(luout,*) "Screen:", screen_fld
c
c      Add E-field to previous value
c
      fld(1) = fld(1) + screen_fld * q_chrg * r(1) / dist3
      fld(2) = fld(2) + screen_fld * q_chrg * r(2) / dist3
      fld(3) = fld(3) + screen_fld * q_chrg * r(3) / dist3
c      write(luout,*) "Field:",fld
c
c      if(ldebug) write(luout,*) "End FldPot"
      end subroutine dimqm_fldpot
