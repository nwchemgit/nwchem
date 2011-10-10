      subroutine argos_space_prtcnt(lcnt)
c
      implicit none
c
#include "argos_space_common.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lcnt(np,3)
c
      integer i,j
c
      do 1 i=1,np
      lcnt(i,1)=0
      lcnt(i,2)=0
      lcnt(i,3)=0
    1 continue
c
      lcnt(me+1,1)=mbxloc
      lcnt(me+1,2)=nwmloc*nwa
      lcnt(me+1,3)=nsaloc
c
      if(np.gt.1) call ga_igop(msp_08,lcnt,3*np,'+')
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' Initial distribution p:b(w+s)',/)
      write(lfnout,1001) (i-1,(lcnt(i,j),j=1,3),i=1,np)
 1001 format(4(3x,i4,':',i5,'(',i7,'+',i7,')'))
      write(lfnout,1002) mwm,msa
 1002 format(/,' Dimension workarrays solvent ',i6,/,
     + 22x,'solute  ',i6)
      endif
c
      return
      end
c $Id$
