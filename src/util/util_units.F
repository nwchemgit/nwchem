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
      character*20 msg

      test1 = unitin
      test2 = unitout
      if (test1 .eq. test2) then
        scale = 1.0
        return
      end if
      test = trim(test1)//"_"//trim(test2)
c     write(*,*) "test=",test
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
      case('ang_au')
        scale = cang2au
      case('au_ang')
        scale = cau2ang
      case('kca_au')
        scale = ckcal2au
      case('au_kca')
        scale = cau2kcal
      case default
        msg = "unknown units: "//trim(test)
        call errquit(msg,0,-1)
      end select
   

      end subroutine
