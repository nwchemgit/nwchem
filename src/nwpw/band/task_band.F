*
* $Id: task_band.F,v 1.1 2003-12-02 19:41:44 bylaska Exp $
*

      logical function task_band(rtdb)
      implicit none
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
c
      integer rtdb
      logical value,dummy
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 1)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 'structure' /

c     *** define external functions ****
      logical  band_structure
      external band_structure
c
c
      value = .false.
      task_band = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task band: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100 ) ind
      call errquit('task band: unimplemented directive', ind, INPUT_ERR)


c
c  band_structure
c
 100  value = .true.
      
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> BAND parallel Module - band_structure <<<<"  
      end if
      value = rtdb_parallel(.true.)
      value = band_structure(rtdb)
      call ga_sync()
      goto 9999



 9999 continue
      dummy = rtdb_parallel(.true.)
      task_band = value

      return
      end
