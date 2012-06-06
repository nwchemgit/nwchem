c
c $Id: util_io.F 19707 2010-10-29 17:59:36Z d3y133 $
c
      subroutine util_length_scale(aunit,scale)
      implicit none
#include "util_params.fh"      
#include "errquit.fh"
      double precision scale
      character*(*) aunit
c
      character*5 test

      test = aunit
      select case (test)
      case('au','angst','bohr')
        scale = cang2au
      case('nm','nanom')
        scale = cnm2au
      case('pm','picom')
        scale = cpm2au
      case default
        call errquit("unknown units: "//test,0,-1)
      end select
   

      end subroutine
