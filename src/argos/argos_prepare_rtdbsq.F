      logical function argos_prepare_rtdbsq(irtdb,igeom,lfnout,sysnam,
     + lseq,cseq,mseq,nseq,latm,catm,xatm,matm,natm)
c
c $Id$
c
c     function to read sequence from the rtdb
c
c     in  : irtdb
c           igeom
c           lfnout      = output file logical file number
c           mseq        = dimension of the sequence list
c           nseq        = expected length of the sequence list
c
c     out : lseq(1,*)   = segment numbers
c           cseq(1,*)   = segment names
c           nseq        = actual length of the sequence list
c           latm(1,*)   = index to sequence list
c           catm(1,*)   = atom name
c           xatm(1:3,*) = atomic coordinates
c
      implicit none
c
#include "geom.fh"
c
      real*8 cau2nm
      parameter(cau2nm=5.29177249d-02)
c
      integer irtdb,igeom
      integer lfnout,mseq,nseq,matm,natm
      integer lseq(6,mseq),latm(5,matm)
      character*6 catm(3,matm)
      character*5 cseq(2,mseq)
      character*80 sysnam
      real*8 xatm(3,matm)
      real*8 c(3),q
      character*16 tag
      integer i,num,numh
      character*10 symbol,element
      integer iatn
c
      do 1 i=1,natm
      if(.not.geom_cent_get(igeom,i,tag,c,q))
     + call md_abort('argos_prepare_rtdbsq: geom_cent_get failed',i)
      latm(1,i)=1
      latm(2,i)=0
      latm(3,i)=0
      latm(4,i)=0
      latm(5,i)=0
      xatm(1,i)=cau2nm*c(1)
      xatm(2,i)=cau2nm*c(2)
      xatm(3,i)=cau2nm*c(3)
      if(.not.geom_tag_to_element(tag,symbol,element,iatn))
     + call md_abort('geom_tag_to_element failed',9999)
      latm(2,1)=iatn
      catm(1,i)='      '
      if(symbol(2:2).eq.' ') then
      catm(1,i)(2:2)=symbol(1:1)
      else
      catm(1,i)(1:2)=symbol(1:2)
      endif
      catm(2,i)='      '
      catm(3,i)='      '
    1 continue
      cseq(1,1)=sysnam(1:5)
      cseq(2,1)='     '
      lseq(1,1)=1
      lseq(2,1)=0
      lseq(3,1)=1
      lseq(3,2)=natm+1
c
      num=0
      numh=0
      do 2 i=1,natm
      if(catm(1,i)(1:2).ne.' H') then
      num=num+1
      if(num.lt.10) then
      write(catm(1,i)(3:3),'(i1.1)') num
      else
      if(num.lt.100) then
      write(catm(1,i)(3:4),'(i2.2)') num
      else
      write(catm(1,i)(3:5),'(i3.3)') num
      endif
      endif
      else
      numh=numh+1
      if(numh.lt.10) then
      write(catm(1,i)(3:3),'(i1.1)') numh
      else
      if(numh.lt.100) then
      write(catm(1,i)(3:4),'(i2.2)') numh
      else
      write(catm(1,i)(3:5),'(i3.3)') numh
      endif
      endif
 
      endif
    2 continue
c
      argos_prepare_rtdbsq=.true.
      return
c
 9999 continue
c
      argos_prepare_rtdbsq=.false.
      return
      end
