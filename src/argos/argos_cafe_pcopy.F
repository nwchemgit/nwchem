      subroutine argos_cafe_pcopy(id,ipari,pari,idimi,idimp,
     + ipar,par,idimt,ip2log,ip3log,ithlog)
c
#include "argos_cafe_common.fh"
c
      integer id,idimi,idimp,idimt
      integer ipari(idimi),ipar(idimt,idimi)
      real*8 pari(4,24),par(idimt,idimp,mset)
      logical ip2log,ip3log,ithlog
c
      do 1 i=1,idimi
      ipar(id,i)=ipari(i)
    1 continue
c
      if(.not.lfree) then
      do 2 i=1,idimp
      do 3 j=1,nparms
      par(id,i,j)=pari(i,j)
    3 continue
    2 continue
      else
      do 4 i=1,idimp
      do 5 j=1,3
      par(id,i,j)=pari(i,j)
    5 continue
      par(id,i,4)=par(id,i,3)-par(id,i,2)
      par(id,i,5)=par(id,i,2)
      par(id,i,6)=par(id,i,3)
      if(abs(pari(i,1)-pari(i,2)).gt.tiny) ip2log=.true.
      if(abs(pari(i,1)-pari(i,3)).gt.tiny) ip3log=.true.
      if(abs(pari(i,2)-pari(i,3)).gt.tiny) ithlog=.true.
    4 continue
      endif
c
      return
      end
c $Id$
