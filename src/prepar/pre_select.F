      logical function pre_select(card,isegm,csa,msa,nsa,
     + maxgrp,maxatm,igroup,lgroup)
c
c $Id: pre_select.F,v 1.8 2001-05-11 20:10:59 d3j191 Exp $
c
      implicit none
c
#include "inp.fh"
c
      logical str_replace
      external str_replace
c
      integer maxgrp,maxatm
      integer igroup(maxgrp,maxatm),lgroup(maxgrp)
c
      integer msa,nsa
      integer isegm(msa)
      character*16 csa(msa)
      character*80 card
c
c
      integer i,j,ifr,ito,igrp
      character*255 target,atomi
      integer ndx,iatom,isgm
      character*4 trgt
c
      card=card(7:80)
c
c     read group number
c
      read(card,'(i5)',err=9999) igrp
      if(igrp.le.0) call errquit('Group number error',igrp)
      if(igrp.gt.maxgrp) call errquit('Group number error',igrp)
c
      card=card(7:80)
c
c     read atom name
c
      ifr=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call errquit('Unable to find group number',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      if(ndx.eq.0) then
      read(target,*,err=5) isgm
      atomi='????'
      goto 6
    5 continue
      isgm=0
      ndx=1
      atomi=target(1:4)
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
    6 continue
      else
      read(target(1:ndx-1),*,err=9999) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call errquit('str_replace',9999)
      endif
      iatom=0
      do 1 i=1,nsa
      trgt=csa(i)(11:14)
      do 4 j=1,4
      if(atomi(j:j).eq.'?') trgt(j:j)='?'
    4 continue
      if((isgm.eq.0.or.isgm.eq.isegm(i)).and.
     + (ndx.eq.0.or.atomi(1:4).eq.trgt)) then
      iatom=iatom+1
      lgroup(igrp)=lgroup(igrp)+1
      if(lgroup(igrp).gt.maxatm)
     + call errquit('Increase group maxatm',maxatm)
      igroup(igrp,lgroup(igrp))=i
      endif
    1 continue
      if(iatom.eq.0) then
      write(*,3) igrp,isgm,atomi(1:4)
    3 format(' Unable to find in select',i5,' atom ',i4,':',a4)
      call errquit('Atom in select not found:',isgm)
      endif
c
      pre_select=.true.
      return
c
 9999 continue
      pre_select=.false.
      return
      end
