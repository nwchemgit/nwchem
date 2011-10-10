      SUBROUTINE ao_tiles(atpart2,nalength2,a2length2,imaxsize2)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "bas.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer atpart2
      integer nalength2(200)
      integer a2length2(200)
      integer imaxsize2
c
      integer ix,i,length
      integer mu,mu_lo,mu_hi,mu_range
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
      do i=1,200
         nalength2(i) = 0
         a2length2(i) = 0
      enddo
c
      length = 0 !current length of a given block (lenght<= imaxsize)
      ix     = 1 ! number of block
      a2length2(ix)=0
      do mu = 1,nshells
         if (.not.bas_cn2bfr(ao_bas_han,mu,mu_lo,mu_hi))
     1   call errquit('tce_ao2e: basis fn range problem 1',0,
     2   BASIS_ERR)
         mu_range = mu_hi - mu_lo + 1
         length = length + mu_range
         if (length.gt.imaxsize2) then
             nalength2(ix)=length
             atpart2 =  ix
             ix=ix+1
             a2length2(ix)=mu
             length = 0
          endif
      enddo  !loop over mu-shells
c
      if (length.gt.0) then
         nalength2(ix)=length
         a2length2(ix+1)=nshells
         atpart2=ix
      endif
      return
      end
c $Id$
