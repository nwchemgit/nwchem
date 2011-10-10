      subroutine argos_diana_rama1(iram,xtrj,trt)
c
      implicit none
c
#include "argos_diana_params.fh"
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer iram(msgm,7)
      real*8 xtrj(ntrj,3,nfrdim),trt(6,nfrdim)
c
      real*8 argos_diana_torsion
      external argos_diana_torsion
c
      real*8 phi,psi,rhb1,rhb2,dist
      integer i,j,k,ihb1,ihb2,ndxram
      logical lphi,lpsi
c
      write(filram,1000) strram(1:index(strram,' ')-1),
     + '.ram '
 1000 format(a,a)
      open(unit=lfnram,file=filram(1:index(filram,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnram)
      ndxram=me*nfrdim*(73*nsgm+20)
      call fseek(lfnram,ndxram,0)
      do 1 k=1,nfrme
c      write(filram,1000) strram(1:index(strram,' ')-1),
c     + me*nfrdim+k,'.ram '
c 1000 format(a,i5.5,a)
      call argos_diana_rama2(iram,xtrj(1,1,k),trt(1,k))
    1 continue
c
      close(unit=lfnram)
c
      return
      end
c $Id$
