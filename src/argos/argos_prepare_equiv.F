      logical function argos_prepare_equiv(csa,isml,isgm,msa,nsa,iequi,
     + msm,nsm)
c
c $Id$
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
c
      integer msa,msm,nsa,nsm
      character*16 csa(msa)
      integer isml(msa),isgm(msa),iequi(msm,4)
c
      integer i,j,ia,ja,number
c
      do 1 i=1,msm
      iequi(i,1)=0
      iequi(i,2)=0
      iequi(i,3)=0
      iequi(i,4)=0
    1 continue
c
      do 2 i=1,nsa
      if(iequi(isml(i),1).eq.0) iequi(isml(i),1)=i
      iequi(isml(i),2)=i
    2 continue
c
      do 3 i=1,nsm-1
      do 4 j=i+1,nsm
      if(iequi(j,3).eq.0) then
      if(iequi(i,2)-iequi(i,1).eq.iequi(j,2)-iequi(j,1)) then
      do 5 ia=iequi(i,1),iequi(i,2)
      ja=ia+iequi(j,1)-iequi(i,1)
      if(csa(ia).ne.csa(ja)) goto 4
    5 continue
      iequi(j,3)=i
      endif
      endif
    4 continue
    3 continue
c
      do 6 i=1,msm
      iequi(i,4)=iequi(i,1)
      iequi(i,1)=0
      iequi(i,2)=0
    6 continue
c
      do 7 i=1,nsa
      if(iequi(isml(i),1).eq.0) iequi(isml(i),1)=isgm(i)
      iequi(isml(i),2)=isgm(i)
    7 continue
c
      do 8 i=1,nsm
      number=1
      do 9 j=1,nsm
      if(iequi(j,3).eq.i) number=number+1
    9 continue
      if(number.gt.1) write(*,1000) number,i
 1000 format(/,i5,' equivalent solute molecules ',i5)
    8 continue
c
      argos_prepare_equiv=.true.
      return
c
c 9999 continue
c      argos_prepare_equiv=.false.
      return
      end
