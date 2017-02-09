      subroutine dimqm_printAtoms(rtdb, xyz)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "global.fh"
c
c   Input variables
      integer rtdb
      double precision xyz(3, nDIM)
c
c   Local variables
      integer tx(nDIM)
      character elements(2,nDIMTypes)
      character*2 ele
      character aele(2)
      character*256 f
      integer i, itype
      integer l_name, k_name
      logical stat
      character*50 dd
      character*50 d

      if(ga_nodeid().eq.0) then
        stat = rtdb_parallel(.false.)
c
c   ===============
c   Allocate Arrays
c   ===============
c
c   Read in atom names
      if(.not.rtdb_get(rtdb,'dimpar:name',mt_byte,nDIMTypes*2,
     $                 elements))
     $  call errquit('dimqm_print get names failed',1,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'dimpar:typeindex',mt_int,nDIM,tx))
     $  call errquit('dimqm_print get tx failed',1,RTDB_ERR)
      write(luout,*) "Elements", elements
      dd = 
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
c
c
c   Print atoms and their coordinates
      write(luout,'(/1x,a)') dd
      write(luout,'(1x,a)')  "DIM COORDINATES"
      write(luout,'(1x,a)')  d
      f = '(17X,A1,11X,A1,11X,A1)'
      write (luout,f) 'X', 'Y', 'Z'
      f = '(1X,A)'
      write (luout,f) d
      f = '(1X,I5,1X,2A1,3F12.5)'
      do i = 1, nDIM
        itype = tx(i)
        write(luout,f) i, elements(1:2,itype), xyz(:,i)
      end do 
      call util_flush(LuOut)
        stat = rtdb_parallel(.true.)
      end if

      end subroutine dimqm_printAtoms
