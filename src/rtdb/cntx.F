      program main
C$Id: cntx.F,v 1.3 1995-02-02 23:21:56 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "context.fh"
#include "rtdb.fh"
c
      integer rtdb
      character*32 name
c
      if (.not. rtdb_open('crap.db', 'empty', rtdb)) stop 1
      if (.not. rtdb_cput(rtdb, 'fred', 1, 'george')) stop 2
      if (.not. rtdb_cget(rtdb, 'fred', 1, name)) stop 3
      write(6,*) name
      name = ' '
      if (.not. context_rtdb_match(rtdb, 'fred', name)) stop 4
      write(6,*) name
c
      end
