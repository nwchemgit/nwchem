      integer function ana_index(card,sgmnam)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*80 card,atom
      character*16 sgmnam(nsa)
c
      integer i,j,k,isgm
c
      i=index(card,':')-1
      if(i.le.0) call errquit('Error in atom name',0)
      read(card(1:i),'(i)') isgm
      atom=card(i+2:80)
      j=index(atom,' ')-1
      if(atom(1:1).eq.'_') atom(1:1)=' '
c
      do 1 i=1,nsa
      read(sgmnam(i)(12:16),'(i5)') k
      if(isgm.eq.k) then
      if(sgmnam(i)(6:5+j).eq.atom(1:j)) then
      ana_index=i
      return
      endif
      endif
    1 continue
c
      call errquit('Error in atom name',1)
      ana_index=0
c
      return
      end
