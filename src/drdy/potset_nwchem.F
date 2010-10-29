C $Id$
      subroutine potset_nwchem(rtdb)
*
* check for basis sets and theoretical method
*
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "rtdb.fh"
*::passed:
      integer rtdb
*::local:
      character*32 theory
      logical okay
*
      okay = .true.
*
      okay = okay.and.bas_rtdb_in(rtdb)
      if (.not.okay) then
        write(luout,*)' no basis set info on rtdb '
        call util_flush(luout)
        call errquit('potset_nwchem: fatal error ',911, RTDB_ERR)
      endif
      if (.not.rtdb_cget(rtdb,'drdy:theory',1,theory)) then
        write(luout,*)' default theory for drdy set to scf '
        call util_flush(luout)
        theory = 'scf'
        if (.not.rtdb_cput(rtdb,'drdy:theory',1,theory))
     &        call errquit('potset_nwchem: could not set theory',911,
     &       RTDB_ERR)
      endif
      end



