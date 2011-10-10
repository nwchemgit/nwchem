      subroutine argos_test
c
      implicit none
c
#include "argos.fh"
c
      itest=itest-1
c
      if(me.ne.0) return
c
      if(ntype.eq.1) then
      write(lfntst,1000) etot
 1000 format('Energy           = ',1pe12.3)
      endif
c
      if(ntype.eq.2) then
      write(lfntst,1100) stime,temp,volume,pres,etot
 1100 format('Time             = ',0pf9.3,/,
     + 'Temperature      = ',0pf8.2,/,
     + 'Volume           = ',0pf8.2,/,
     + 'Pressure         = ',1pe12.2,/,
     + 'Energy           = ',1pe12.3,/)
      endif
c
      if(ntype.eq.3) write(lfntst,1200) isdit+icgit,etot
 1200 format('Iteration        = ',i10,/,
     + 'Energy           = ',1pe12.3,/)
c
      if(itest.eq.0) close(lfntst)
c
      return
      end
c $Id$
