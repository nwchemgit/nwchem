      logical function pre_disres(card,card2,isegm,csa,msa,nsa,
     + inoe,dnoe,mnoe,nnoe)
c
c $Id$
c
      implicit none
c
#include "inp.fh"
c
      logical str_replace
      external str_replace
c
      integer mnoe,nnoe,msa,nsa
      integer isegm(msa),inoe(2,mnoe)
      real*8 dnoe(5,mnoe)
      character*16 csa(msa)
      character*80 card,card2
      integer i,ifr,ito
      character*255 target,atomi,atomj
      integer ndx,iatom,jatom,isgm,jsgm
c
      card=card(8:80)
c
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find atom',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      read(target(1:ndx-1),*) jsgm
      atomj=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomj,'_',' '))
     + call md_abort('str_replace',9999)
      iatom=0
      jatom=0
      do 1 i=1,nsa
      if(isgm.eq.isegm(i).and.atomi(1:4).eq.csa(i)(11:14)) iatom=i
      if(jsgm.eq.isegm(i).and.atomj(1:4).eq.csa(i)(11:14)) jatom=i
    1 continue
      if(iatom.eq.0.or.jatom.eq.0)
     + call md_abort('pre_disres error',9999)
      nnoe=nnoe+1
      if(nnoe.gt.mnoe) call md_abort('Increase mnoe',9999)
      inoe(1,nnoe)=iatom
      inoe(2,nnoe)=jatom
      read(card2,1000) (dnoe(i,nnoe),i=1,5)
 1000 format(3f12.6,2e12.5)
c
      pre_disres=.true.
      return
c
 9999 continue
      pre_disres=.false.
      return
      end
