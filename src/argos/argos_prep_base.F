      integer function argos_prep_base(isgm,isegm,csa,msa,nsa)
c
c $Id$
c
      implicit none
c
#include "inp.fh"
c
      logical argos_prep_str_replace
      external argos_prep_str_replace
c
      integer msa,nsa,isgm
      integer isegm(msa)
      character*16 csa(msa)
c
      integer i,ibase
c
c     determine base
c
c     0: none
c     1: adenine
c     2: guanine
c     3: cytosine
c     4: uracil
c     5: thymine
c
      ibase=0
c
      do 1 i=1,nsa
      if(isgm.eq.isegm(i)) then
      if(csa(i)(1:3).eq.'DA ') ibase=1
      if(csa(i)(1:3).eq.'DA_') ibase=1
      if(csa(i)(1:3).eq.'DAT') ibase=1
      if(csa(i)(1:3).eq.'DG ') ibase=2
      if(csa(i)(1:3).eq.'DG_') ibase=2
      if(csa(i)(1:3).eq.'DGT') ibase=2
      if(csa(i)(1:3).eq.'DC ') ibase=3
      if(csa(i)(1:3).eq.'DC_') ibase=3
      if(csa(i)(1:3).eq.'DCT') ibase=3
      if(csa(i)(1:3).eq.'RU ') ibase=4
      if(csa(i)(1:3).eq.'RU_') ibase=4
      if(csa(i)(1:3).eq.'RUT') ibase=4
      if(csa(i)(1:3).eq.'DT ') ibase=5
      if(csa(i)(1:3).eq.'DT_') ibase=5
      if(csa(i)(1:3).eq.'DTT') ibase=5
      goto 2
      endif
    1 continue
    2 continue
c
      argos_prep_base=ibase
      return
      end
