      subroutine argos_cafe_rddrs(lfntop,idrs,rdrs)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfntop
      integer idrs(ndrs,2)
      real*8 rdrs(ndrs,6)
c
      integer i,j
c
      if(me.eq.0) then
      do 1 i=1,ndrs
      read(lfntop,1000) (idrs(i,j),j=1,2),(rdrs(i,j),j=1,5)
 1000 format(2i10,3f12.6,/,20x,2e12.5)
      rdrs(i,6)=zero
    1 continue
      endif
c
      call ga_sync()
      call ga_brdcst(mcf_57,idrs,ma_sizeof(mt_int,ndrs*2,mt_byte),0)
      call ga_brdcst(mcf_58,rdrs,ma_sizeof(mt_dbl,ndrs*6,mt_byte),0)
c
      return
      end
c $Id$
