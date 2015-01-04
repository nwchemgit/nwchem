c
c $Id$
c
      subroutine util_length_scale(aunit,scale)
      implicit none
#include "util_params.fh"      
#include "errquit.fh"
      double precision scale
      character*(*) aunit
c
      character*3 test

      test = aunit
      select case (test)
      case('au','boh')
        scale = 1.0d0
      case('ang')
        scale = cang2au
      case('nm','nan')
        scale = cnm2au
      case('pm','pic')
        scale = cpm2au
      case default
        call errquit("unknown units: "//test,0,-1)
      end select
   

      end subroutine

