      subroutine qhop_pindex(msa,nsa,isgm,ihop,ips)
c
c     in:  msa  : array dimension
c          nsa  : number of atoms
c          isgm : segment index
c          ihop : flag for hoppability
c
c     out: ips  : protonation index
c
      implicit none
c
#include "qhop_common.fh"
c
      integer msa,nsa
      integer isgm(msa),ips(msa),ihop(msa)
c
      integer i
c
      do 1 i=1,nsa
      if(ihop(i).lt.0)then
      ips(i)=arat(i-ptseq(isgm(i),2)+ptarat(ptseq(isgm(i),1)),3)
      else
      ips(i)=0
      endif
    1 continue
c  
      return
      end
