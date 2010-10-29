      logical function argos_prepare_select(card,isegm,ismol,csa,msa,
     + nsa,
     + maxgrp,maxatm,igroup,lgroup,mgroup)
c
c $Id$
c
      implicit none
c
#include "inp.fh"
#include "mafdecls.fh"
c
      logical str_replace
      external str_replace
c
      integer maxgrp,maxatm
      integer igroup(maxgrp,maxatm),lgroup(maxgrp),mgroup(maxgrp)
c
      integer msa,nsa
      integer isegm(msa),ismol(msa)
      character*16 csa(msa)
      character*80 card
c
c
      integer i,j,ifr,ito,igrp
      character*255 target,atomi
      integer ndx,iatom,isgm
      character*4 trgt
      integer i_list,l_list
c
      card=card(7:80)
c
c     read group number
c
      read(card,'(i5)',err=9999) igrp
      if(igrp.le.0) call md_abort('Group number error',igrp)
      if(igrp.gt.maxgrp) call md_abort('Group number error',igrp)
c
      card=card(7:80)
c
      ndx=index(card,' ')
      target=card(1:ndx-1)
c
      if(.not.ma_push_get(mt_int,nsa,'list',l_list,i_list))
     + call md_abort('Failed to allocate list',0)
c
      call argos_prepare_atom_selection(csa,isegm,nsa,card,
     + int_mb(i_list))
c
c
      iatom=0
      do 1 i=1,nsa
      if(int_mb(i_list-1+i).eq.1) then
      do 2 j=1,lgroup(igrp)
      if(igroup(igrp,j).eq.i) goto 1
 2    continue
      iatom=iatom+1
      lgroup(igrp)=lgroup(igrp)+1
      if(iatom.eq.1) mgroup(igrp)=ismol(i)
      if(mgroup(igrp).ne.ismol(i)) mgroup(igrp)=0
      if(lgroup(igrp).gt.maxatm)
     + call md_abort('Increase group maxatm',maxatm)
      igroup(igrp,lgroup(igrp))=i
      endif
 1    continue
c
      write(*,'(a,i5,a,i5,a)') ' Selection group',igrp,' contains',
     + lgroup(igrp),' atoms'
c
      if(.not.ma_pop_stack(l_list))
     + call md_abort('Failed to deallocate list',0)
c
      argos_prepare_select=.true.
c
      return
c
 9999 continue
      argos_prepare_select=.false.
c
      stop
c      return
      end
