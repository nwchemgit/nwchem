      subroutine argos_diana_equivalent(sgmnam,imol,iequi,itmp)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      character*16 sgmnam(msa)
      integer imol(msa),iequi(msa),itmp(nsa,2)
c
      integer i,ism,jsm,nsm,nsmi,nsmj,isa,jsa
c
      if(me.eq.0) then
c
      do 1 i=1,nsa
      iequi(i)=0
      itmp(imol(i),2)=i
    1 continue
      nsm=imol(nsa)
      itmp(1,1)=1
      do 2 i=2,nsm
      itmp(i,1)=itmp(i-1,2)+1
    2 continue
c
      do 3 ism=1,nsm-1
      nsmi=itmp(ism,2)-itmp(ism,1)+1
      do 4 jsm=ism+1,nsm
      if(nsmi.eq.itmp(jsm,2)-itmp(jsm,1)+1) then
      do 5 i=1,nsmi
      isa=itmp(ism,1)+i-1
      jsa=itmp(jsm,1)+i-1
      if(sgmnam(isa)(1:10).ne.sgmnam(jsa)(1:10)) goto 4
    5 continue
      jsa=itmp(jsm,1)-1
      do 6 i=itmp(ism,1),itmp(ism,2)
      iequi(i)=jsa+i-itmp(ism,1)+1
    6 continue
      goto 3
      endif
    4 continue
    3 continue
c
c      write(*,'(3i7,1x,a)') (i,imol(i),iequi(i),sgmnam(i),i=1,nsa)
c
      endif
      if(np.gt.1) then
      call ga_brdcst(mag_d06,iequi,nsa*ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      return
      end
c $Id$
