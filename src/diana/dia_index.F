      subroutine dia_index(sgmnam,ia,ja)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*80 atom
      character*16 sgmnam(nsa)
      character*10 word
      integer ia,ja
c
      integer i,j,k,isgm
c
      ia=0
      ja=0
c
      if(card(1:1).ne.'W'.and.card(1:1).ne.'w') then
c
      i=index(card,':')-1
      if(i.le.0.or.i.gt.10) call md_abort('Error in atom name',0)
      word='          '
      word(11-i:10)=card(1:i)
      read(word,'(i10)') isgm
      atom=card(i+2:80)
      j=index(atom,' ')-1
      if(atom(1:1).eq.'_') atom(1:1)=' '
c
      do 1 i=1,nsa
      read(sgmnam(i)(12:16),'(i5)') k
      if(isgm.eq.k) then
      if(sgmnam(i)(6:10).eq.atom(1:5)) then
      ia=i
      return
      endif
      endif
    1 continue
      print*,'index for atom is ',ia
c
      else
c
      i=index(card,':')-1
      if(i.le.1) call md_abort('Error in atom name',0)
      if(i.eq.2) then
      read(card(2:i),'(i1)') ia
      elseif(i.eq.3) then
      read(card(2:i),'(i2)') ia
      elseif(i.eq.4) then
      read(card(2:i),'(i3)') ia
      elseif(i.eq.5) then
      read(card(2:i),'(i4)') ia
      elseif(i.eq.6) then
      read(card(2:i),'(i5)') ia
      else
      call md_abort('Error in dia_index',me)
      endif
      ia=ia-nsgm
      if(ia.le.0) call md_abort('Error in dia_index',me)  
      atom=card(i+2:80)
      j=index(atom,' ')-1
      if(atom(1:1).eq.'_') atom(1:1)=' '
c
      do 2 i=1,nwa
      if(wnam(i)(6:10).eq.atom(1:5)) then
      ja=i
      return
      endif
    2 continue
c
      endif
c
      write(atom,'(a,a,a)') 'Atom ',card(1:index(card,' ')-1),
     + ' could not be found' 
c
      call md_abort(atom,0)
c
      return
      end
c $Id$
