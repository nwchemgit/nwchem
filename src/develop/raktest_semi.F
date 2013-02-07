#if defined(SEMI)
      subroutine raktest_semi(rtdb)
* $Id$
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
c#include "se_defs.fh"
c
c::passed
      integer rtdb ![input] rtdb handle
c::local
      character*30 mode_name
      integer se_mode
c
c in input deck put: 
*    set "semi:mode" "AM1"
      mode_name = 'MNDO'
      if (.not.rtdb_cget(rtdb,'semi:mode',1,mode_name))
     &    mode_name='MNDO'

      if (inp_compare(.false.,'am1',mode_name)) then
        write(luout,*)' am1 mode '
      elseif (inp_compare(.false.,'mndo',mode_name)) then
        write(luout,*)' am1 mode '
      else
        write(luout,*)' unknown mode '
      endif
      write(luout,*)' mode is ',mode_name
c--------------------------------------------
      if(mode_name.eq.'AM1')se_mode=SE_AM1
      if(mode_name.eq.'MNDO') se_mode=SE_MNDO
      call se_init(rtdb,se_mode)
c--------------------------------------------
c
      end
#else
      subroutine raktest_semi(rtdb)
      implicit none
      integer rtdb
      end
#endif
