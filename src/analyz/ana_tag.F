      subroutine ana_tag(card,sgmnam,x)
c
c $Id: ana_tag.F,v 1.7 2000-08-17 22:14:00 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*80 card
      character*16 sgmnam(nsa)
      real*8 x(nsa,3)
c
      call ana_settag(card,sgmnam,byte_mb(i_tag),x)
c
      return
      end
      subroutine ana_settag(card,sgmnam,tag,x)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*80 card
      character*16 sgmnam(nsa)
      character*24 ttag,tag(nsa,2)
      real*8 x(nsa,3)
c
      integer i,j,k,itag
      integer is,js,it
      character*5 cs,ct
      real*8 rtag,rtag2,r2
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card,6901,end=9999) itag,ttag,is,js,rtag,cs
 6901 format(7x,i3,a24,2i7,f12.6,a5)
      write(*,6902) itag,ttag(1:index(ttag,' ')-1),
     + cs(1:index(cs,' ')),is,js,rtag
 6902 format(' Label ',i2,' set to ',a,' for atoms ',a,
     + ' in segments',
     + i7,' through ',i7,' and atoms within ',f12.6, ' nm')
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
cx      print*,sgmnam(i)
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
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d06,byte_mb(i_tag),48*nsa,0)
      endif
c
      return
c
 9999 continue
c
      call errquit('Error in settag',0)
      return
      end
      subroutine ana_initag(tag)
c
c     implicit none
c
#include "ana_common.fh"
c
      character*24 tag(nsa,2)
c
      integer i
c
      do 1 i=1,nsa
      tag(i,1)='                        '
      tag(i,2)='                        '
    1 continue
c
      return
      end
