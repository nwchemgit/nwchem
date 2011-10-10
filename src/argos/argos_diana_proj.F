      subroutine argos_diana_proj(ndex,vct,isel,xtrj,xprj,xa,val)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer ndex(nsel)
      real*8 vct(*),val(nfrdim),pnorm,vval
      real*8 xtrj(ntrj,3,nfrdim),xprj(ntrj,3,nfrdim),xa(nsa,3)
      integer isel(nsa)
c
      integer i,j,k,l,m
      real*8 value
c
c      write(*,'(3f12.6)') ((xtrj(5,i,k),i=1,3),k=1,5)
c      write(*,'(5i5)') (ndex(i),i=1,5)
c      write(*,'(6f12.6)') (vct(k),k=1,3*nsel)
c      write(*,'(3f12.6)') (xa(5,i),i=1,3)
      do 1 k=1,nfrme
      val(k)=0.0d0
      pnorm=0.0d0
      m=0
      do 2 i=1,nsel
      vval=0.0d0
      do 3 j=1,3
      m=m+1
      val(k)=val(k)+(xtrj(ndex(i),j,k)-xa(ndex(i),j))*vct(m)
      pnorm=pnorm+vct(m)*vct(m)
    3 continue
    2 continue
      val(k)=val(k)/sqrt(pnorm)
    1 continue
c      write(*,'(5f12.6)') (val(k),k=1,nfr)
c
      do 4 k=1,nfrme
      m=0
      do 5 i=1,nsel
      do 6 j=1,3
      m=m+1
      xprj(ndex(i),j,k)=val(k)*vct(m)+xa(ndex(i),j)
c      write(*,'(4i5,3f12.6)') k,i,j,ndex(i),val(k),vct(m),xa(ndex(i),j)
    6 continue
    5 continue
    4 continue
c
c     write the vector
c
      if(me.eq.0) then
      open(unit=lfntcl,file=filtcl(1:index(filtcl,' ')-1),
     + form='formatted',status='unknown',err=9999)
      k=0
      l=0
      do 15 i=1,nsel
      value=0.0d0
      do 16 j=1,3
      value=value+(dbl_mb(i_cov+k))**2
      k=k+1
   16 continue
      l=l+1
      write(lfntcl,5555)
     + (10.0*xa(ndex(i),j),j=1,3),
     + (dbl_mb(i_cov+k-4+j),j=1,3)
 5555 format('draw vector2 {',3f12.6,'} {',3f12.6,'} 1.0 10')
      write(lfnamp,'(i5,f12.6)') i,sqrt(value)
   15 continue
      close(unit=lfntcl)
 9999 continue
      endif
c      write(*,'(3f12.6)') ((xprj(5,i,k),i=1,3),k=1,5)
c
c      if(valmin.gt.valmax) then
c      call argos_diana_wrtpdb(lfnmax,filmax,byte_mb(i_snam),isel,xp,
c     + isel,xp,.true.,.false.,int_mb(i_ndxw))
c      valmax=val
c      call argos_diana_wrtpdb(lfnmin,filmin,byte_mb(i_snam),isel,xp,
c     + isel,xp,.true.,.false.,int_mb(i_ndxw))
c      valmin=val
c      elseif(val.lt.valmin) then
c      call argos_diana_wrtpdb(lfnmin,filmin,byte_mb(i_snam),isel,xp,
c     + isel,xp,.true.,.false.,int_mb(i_ndxw))
c      valmin=val
c      elseif(val.gt.valmax) then
c      call argos_diana_wrtpdb(lfnmax,filmax,byte_mb(i_snam),isel,xp,
c     + isel,xp,.true.,.false.,int_mb(i_ndxw))
c      valmax=val
c      endif
c
      return
      end
c $Id$
