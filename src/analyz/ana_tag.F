      subroutine ana_tag(sgmnam)
c
c $Id: ana_tag.F,v 1.1 1999-03-15 21:27:48 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
c
      call ana_settag(sgmnam,byte_mb(i_tag),dbl_mb(i_xref))
c
      return
      end
      subroutine ana_settag(sgmnam,tag,x)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      character*24 ttag,tag(nsa,2)
      real*8 x(nsa,3)
c
      integer i,j,k,itag
      character*255 filtag,fpref
      integer is,js,it
      character*5 cs,ct
      real*8 rtag,rtag2,r2
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      do 1 i=1,nsa
      tag(i,1)='                        '
      tag(i,2)='                        '
    1 continue
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,fpref))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      filtag=fpref(1:index(fpref,' ')-1)//'.tag '
      open(unit=69,file=filtag(1:index(filtag,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(unit=69)
c
  100 continue
      read(69,6901,end=9998) itag,ttag,is,js,rtag,cs
 6901 format(i3,a24,2i7,f12.6,a5)
      write(*,6901) itag,ttag,is,js,rtag,cs
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
      read(sgmnam(i),'(5x,a5,i6)') ct,it
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    4 continue
      if((is.eq.0.or.(is.le.it.and.js.ge.it)).and.
     + (cs.eq.'     '.or.cs.eq.ct)) then
      tag(i,itag)=ttag
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 5 k=1,nsa
      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
      if(r2.le.rtag2) tag(k,itag)=ttag
    5 continue
      endif
      endif
    3 continue
      goto 100
c
 9998 continue
      close(unit=69)
c
 9999 continue
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d06,tag,48*nsa,0)
      endif
c
      return
      end
