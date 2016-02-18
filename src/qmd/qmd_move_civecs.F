c
c     qmd_move_civecs
c
      subroutine qmd_move_civecs(rtdb,ipol,mult)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
c
#include "qmd_common.fh"
c
      integer rtdb             ! RTDB handle
      integer ipol
      integer mult
c
c     local variables
      logical vector_exists
      logical status
      integer len_fn_civecs           ! length of civecs file name
      character*255 fn_civecs         ! file name for TDDFT vectors
      character*255 fn_civecs_current ! file name for current CI vectors
      character*255 fn_civecs_old     ! file name for previous CI vectors
c
      character*30 pname
c
c     preliminaries
      pname='qmd_move_civecs: '
c
c     CI Vectors files for restricted singlets
      if (rtdb_cget(rtdb,'tddft:civecs',1,fn_civecs)) then
        len_fn_civecs = inp_strlen(fn_civecs)
        if (ipol.eq.1) then
          fn_civecs_current = fn_civecs(1:len_fn_civecs)//"_singlet"
        else
          fn_civecs_current = fn_civecs(1:len_fn_civecs)
        end if
        len_fn_civecs = inp_strlen(fn_civecs_current)
        fn_civecs_old = fn_civecs_current(1:len_fn_civecs)//"_previous"
        if (ga_nodeid().eq.0) 
     $    call util_file_copy(fn_civecs_current,fn_civecs_old)
      else
        call errquit(pname//'no saved tddft results?',0,CALC_ERR)
      endif
      call ga_sync()
c
      return
      end
