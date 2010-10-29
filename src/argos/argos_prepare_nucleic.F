      logical function argos_prepare_nucleic(sgmnam)
c
c $Id$
c
      implicit none
c
#include "argos_prepare_common.fh"
c
      character*5 sgmnam
c
      argos_prepare_nucleic=.true.
c
      if(sgmnam(1:3).eq.'A  '.or.sgmnam(1:3).eq.'+A '.or.
     + sgmnam(1:3).eq.'C  '.or.sgmnam(1:3).eq.'+C '.or.
     + sgmnam(1:3).eq.'G  '.or.sgmnam(1:3).eq.'+G '.or.
     + sgmnam(1:3).eq.'I  '.or.sgmnam(1:3).eq.'+I '.or.
     + sgmnam(1:3).eq.'T  '.or.sgmnam(1:3).eq.'+T '.or.
     + sgmnam(1:3).eq.'U  '.or.sgmnam(1:3).eq.'+U ') return
      if(ffield(1:5).eq.'amber') then
      if(sgmnam(1:3).eq.'DA '.or.sgmnam(1:3).eq.'RA '.or.
     + sgmnam(1:3).eq.'DC '.or.sgmnam(1:3).eq.'RC '.or.
     + sgmnam(1:3).eq.'DG '.or.sgmnam(1:3).eq.'RG '.or.
     + sgmnam(1:3).eq.'DT '.or.sgmnam(1:3).eq.'RU ') return
      endif
      if(ffield(1:6).eq.'charmm') then
      endif
c
      argos_prepare_nucleic=.false.
      return
      end
