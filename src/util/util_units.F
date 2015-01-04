c
c $Id: util_length_scale.F 24349 2013-06-23 04:30:58Z d3y133 $
c
      subroutine util_convert_units(unitin,unitout,scale)
      implicit none
#include "util_params.fh"      
#include "errquit.fh"
      double precision scale
      character*(*) unitin,unitout
c
      character*3 test1,test2
      character*7 test

      test1 = unitin
      test2 = unitout
      test = trim(test1)//"_"//trim(test2)
      write(*,*) "test=",test
      select case (test)
      case('au_boh')
        scale = 1.0d0
      case('ang_boh')
        scale = cang2au
      case('nm_au')
        scale = cnm2au
      case('nan_au')
        scale = cnm2au
      case('pm_au')
        scale = cpm2au
      case('pic_au')
        scale = cpm2au
      case('au_ang')
        scale = cau2ang
      case('kcal_au')
        scale = ckcal2au
      case('au_kcal')
        scale = cau2kcal
      case default
        call errquit("unknown units: "//test,0,-1)
      end select
   

      end subroutine
