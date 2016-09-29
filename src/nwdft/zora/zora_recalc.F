c
c     determine if zora correction needs to be recalculated
      logical function zora_recalc(rtdb)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer rtdb
c
      character*32 operation
      character*30 pname
c
c     preliminaries
      pname = "zora_dorecal: "
c
      zora_recalc = .false.
      if (.not.rtdb_get(rtdb, 'dft:zora_recalc', mt_log, 1,
     &       zora_recalc)) zora_recalc = .false.
      if (zora_recalc) return
c
c     get operation
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $ operation = 'energy'
c
      zora_recalc = (operation.eq.'opt').or.(operation.eq.'optimize')
     & .or. (operation.eq.'saddle').or. (operation.eq.'neb').or. 
     &  (operation.eq.'string')
c
      return
      end
