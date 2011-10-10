      subroutine argos_cafe_copya(ia,nami,wgti,numa,nam,wgt,num,ityp)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia
      character*6 nami(24),nam(mat,nparms)
      real*8 wgti(24),wgt(mat,mset)
      integer numa(24),num(mat,nparms),ityp(mat,nparms)
c
      integer i
c
      if(.not.lfree) then
      do 1 i=1,nparms
      nam(ia,i)=nami(i)
      wgt(ia,i)=wgti(i)
      num(ia,i)=numa(i)
      ityp(ia,i)=0
      if(nam(ia,i)(6:6).eq.'Q') ityp(ia,i)=1
      if(nam(ia,i)(6:6).eq.'H') ityp(ia,i)=2
    1 continue
      else
      do 2 i=1,3
      nam(ia,i)=nami(i)
      wgt(ia,i)=wgti(i)
      num(ia,i)=numa(i)
      ityp(ia,i)=0
      if(nam(ia,i)(6:6).eq.'Q') ityp(ia,i)=1
      if(nam(ia,i)(6:6).eq.'H') ityp(ia,i)=2
    2 continue
      wgt(ia,4)=wgt(ia,3)-wgt(ia,2)
      wgt(ia,5)=wgt(ia,2)
      wgt(ia,6)=wgt(ia,3)
      endif
c
      return
      end
c $Id$
