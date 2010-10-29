      logical function argos_prep_paths2(irtdb)
c
c $Id$
c
      implicit none
c
      logical util_nwchemrc_get
      external util_nwchemrc_get
c
#include "rtdb.fh"
#include "argos_prep_common.fh"
c
      integer irtdb
      integer i,j,len,lend,ifc
      character*255 key,value,string
      logical add
c
c      print*,'ARGOS_PREP_PATHS2'
      ndirpar=9
      do 1 i=1,mdirpar
      add=.false.
      write(key,'(a,i1)') 'prep:dir_',i
      if(rtdb_cget(irtdb,key,1,value)) then
      len=index(value,' ')-1
      if(len.gt.0) then
      lend=len
      if(len.gt.4) then
      if(value(len-3:len).eq.'.par') then
      do 2 j=1,len
      if(value(j:j).eq.'/') then
      lend=j
      if(j.lt.len.and.value(j+1:j+1).eq.'+') add=.true.
      endif
    2 continue
      endif
      endif
      string=value(1:lend)
      if(value(lend:lend).ne.'/') then
      string(lend+1:lend+1)='/'
      endif
c
      if(add) then
      if(dirpar(i)(1:index(string,' ')).ne.string(1:index(string,' ')))
     + call md_abort('Parameter file added from different directory',0)
      lend=lend+1
      else
      dirpar(i)=string(1:index(string,' '))
      endif
c
      if(len.eq.lend) then
      nfilpar(i)=1
      filpar(i,1)=ffield(1:index(ffield,' ')-1)//'.par '
      else
      if(.not.add) nfilpar(i)=0
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
c      do 4 i=1,mdirpar
c      if(nfilpar(i).gt.0) then
c      do 5 j=1,nfilpar(i)
c      write(*,'(a,a)') dirpar(i)(1:index(dirpar(i),' ')-1),
c     + filpar(i,j)(1:index(filpar(i,j),' ')-1)
c    5 continue
c      endif
c    4 continue
c
      argos_prep_paths2=.true.
      return
      end
