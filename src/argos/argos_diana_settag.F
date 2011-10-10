      subroutine argos_diana_settag(sgmnam,tag,val,x)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      character*24 ttag,tag(nsa,2)
      real*8 x(nsa,3),val(nsa,2)
c
      integer i,j,k,itag
      integer is,js,it
      character*11 cato
      character*5 cs,ct,csgm,cst
      real*8 rtag,rtag2,r2,rval
c
c     do selection on node 0
c
      if(me.eq.0) then
c
      read(card,6901,end=9999) itag,ttag,rval,is,js,rtag,cato
 6901 format(7x,i3,a24,f12.6,2i7,f12.6,a11)
      csgm='     '
      k=index(cato,':')-1
      if(k.le.0) k=0
      if(k.gt.0) then
      if(k.le.5) csgm(1:k)=cato(1:k)
      k=k+1
      endif
      cs(1:5)=cato(k+1:k+5)
      if(is.ne.0.and.js.ne.0) then
      if(is.gt.0) then
      write(*,6902) itag,ttag(1:index(ttag,' ')-1),
     + cs(1:index(cs,' ')),is,js,rtag
 6902 format(' Label ',i2,' set to ',a,' for atoms ',a,
     + ' in segments',
     + i7,' through ',i7,' and atoms within ',f12.6, ' nm')
      else
      write(*,6903) itag,ttag(1:index(ttag,' ')-1),
     + cs(1:index(cs,' ')),rtag
 6903 format(' Label ',i2,' set to ',a,' for atoms ',a,
     + ' in all segments and atoms within ',f12.6, ' nm')
      endif
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 3 i=1,nsa
      read(sgmnam(i),'(a5,a5,i6)') cst,ct,it
      do 4 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
      if(csgm(j:j).eq.'?') cst(j:j)='?'
    4 continue
      if((is.eq.-1.or.(is.le.it.and.js.ge.it)).and.
     + (cs.eq.'     '.or.cs.eq.ct).and.
     + (csgm.eq.'     '.or.csgm.eq.cst)) then
      tag(i,itag)=ttag
      val(i,itag)=rval
      if(rtag.gt.0.0d0) then
      rtag2=rtag*rtag
      do 5 k=1,nsa
      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
      if(r2.le.rtag2) then
      tag(k,itag)=ttag
      val(k,itag)=rval
      endif
    5 continue
      endif
      endif
    3 continue
c
      else
      write(*,6904) itag,ttag(1:index(ttag,' ')-1),
     + cs(1:index(cs,' '))
 6904 format(' Label ',i2,' set to ',a,' for solvent atoms ',a)
c     + ' in solvent and solvent atoms within ',f12.6, ' nm')
      if(cs(1:1).eq.'_') cs(1:1)=' '
      do 6 i=1,nwa
      read(wnam(i),'(5x,a5,i6)') ct,it
      do 7 j=1,5
      if(cs(j:j).eq.'?') ct(j:j)='?'
    7 continue
      if(cs.eq.'     '.or.cs.eq.ct) then
      wtag(i,itag)=ttag
      wval(i,itag)=rval
c      if(rtag.gt.0.0d0) then
c      rtag2=rtag*rtag
c      do 8 k=1,nsa
c      r2=(x(k,1)-x(i,1))**2+(x(k,2)-x(i,2))**2+(x(k,3)-x(i,3))**2
c      if(r2.le.rtag2) tag(k,itag)=ttag
c    8 continue
c      endif
      endif
    6 continue
 
      endif
c
      endif
c
c     broadcast to all nodes
c
#if !defined(WIN32)
      if(np.gt.1) then
      call ga_brdcst(mag_d06,byte_mb(i_tag),48*nsa,0)
      call ga_brdcst(mag_d06,dbl_mb(i_tag),16*nsa,0)
      endif
#endif
c
      return
c
 9999 continue
c
      call md_abort('Error in settag',0)
      return
      end
c $Id$
