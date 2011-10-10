      subroutine qhop_pindex(msa,nsa,isgm,ihop,ips)
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
      integer msa,nsa
      integer isgm(msa),ips(msa),ihop(msa)
c
      call qhop_pindex2(msa,nsa,isgm,ihop,ips,
     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
     + dbl_mb(i_par+28*maxpar))
c
      return
      end
      subroutine qhop_pindex2(msa,nsa,isgm,ihop,ips,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy)
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
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
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
c $Id$
