      integer function argos_prepare_atom_index(card,isegm,csa,msa,nsa)
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
      integer msa,nsa
      integer isegm(msa)
      character*16 csa(msa)
      character*80 card
c
      integer i,ifr,ito
      character*255 target,atomi
      integer ndx,iatom,isgm
c
c     read atom name
c
      ifr=0
      iatom=0
      if(.not.inp_strtok(card,' ',ifr,ito))
     + call md_abort('Unable to find group number',9999)
      target=card(ifr:ito)//' '
      ndx=index(target,':')
      if(ndx.eq.0) then
      read(target,*,err=9999) isgm
      atomi='    '
      else
      read(target(1:ndx-1),*,err=9999) isgm
      atomi=target(ndx+1:index(target,' ')-1)//'      '
      if(.not.str_replace(atomi,'_',' '))
     + call md_abort('str_replace',9999)
      endif
      iatom=0
      do 1 i=1,nsa
      if(isgm.eq.isegm(i).and.
     + (ndx.eq.0.or.atomi(1:4).eq.csa(i)(11:14))) then
      iatom=i
      goto 2
      endif
    1 continue
    3 format(' Unable to find  atom ',i4,':',a4)
      call md_abort('Atom in atom_index not found:',isgm)
    2 continue
c
      argos_prepare_atom_index=iatom
      return
c
 9999 continue
      argos_prepare_atom_index=iatom
      return
      end
