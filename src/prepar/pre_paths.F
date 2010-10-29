      logical function pre_paths()
c
c $Id$
c
      implicit none
c
      logical util_nwchemrc_get
      external util_nwchemrc_get
c
#include "pre_common.fh"
c
      integer i,j,len,lend,ifc
      character*255 key,value
c
      ndirpar=9
      do 1 i=1,mdirpar
      nfilpar(i)=0
      write(key,'(a,a,i1,a)') ffield(1:index(ffield,' ')-1),'_',i,' '
      if(util_nwchemrc_get(key,value)) then
      len=index(value,' ')-1
      if(len.gt.0) then
      lend=len
      if(len.gt.4) then 
      if(value(len-3:len).eq.'.par') then
      do 2 j=1,len
      if(value(j:j).eq.'/') lend=j
    2 continue
      endif
      endif
      dirpar(i)=value(1:lend)
      if(value(lend:lend).ne.'/') dirpar(i)(lend+1:lend+1)='/'
      if(len.eq.lend) then
      nfilpar(i)=1
      filpar(i,1)=ffield(1:index(ffield,' ')-1)//'.par '
      else
      nfilpar(i)=0
      ifc=lend+1
      do 3 j=lend+1,len
      if(value(j:j).eq.',') then
      nfilpar(i)=nfilpar(i)+1
      filpar(i,nfilpar(i))=value(ifc:j-1)//' '
      ifc=j+1
      elseif(j.eq.len) then
      nfilpar(i)=nfilpar(i)+1
      filpar(i,nfilpar(i))=value(ifc:len)//' '
      endif
    3 continue
      endif
      endif
      endif
    1 continue
c
      if(nfilpar(mdirpar).eq.0) then
      nfilpar(mdirpar)=1
      dirpar(mdirpar)='./'
      filpar(mdirpar,1)=ffield(1:index(ffield,' ')-1)//'.par'
      endif
c
c      do 4 i=1,mdirpar
c      if(nfilpar(i).gt.0) then
c      do 5 j=1,nfilpar(i)
c      write(*,'(a,a)') dirpar(i)(1:index(dirpar(i),' ')-1),
c     + filpar(i,j)(1:index(filpar(i,j),' ')-1)
c    5 continue
c      endif
c    4 continue
c
      pre_paths=.true.
      return
      end
