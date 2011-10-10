      logical function argos_prep_master(lfnfrg,lfnout,nparms,nparmsp,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd,izmat,zmat,mzmat,nzmat,
     + dir,fils,ndx,ndx2)
c
      implicit none
c
#include "util.fh"
#include "argos_prep_common.fh"
c
      integer lfnfrg,lfnout,matm,natm,mbnd,nbnd,mzmat,nzmat
      integer nparms,nparmsp
      integer latm(6,matm),lbnd(2,mbnd),izmat(mzmat,4)
      integer ndx(matm),ndx2(matm)
      character*6 catm(nparmsp,matm)
      real*8 qatm(nparms,2,matm),zmat(mzmat,3)
      character*255 dir,fils(nparmsp)
c
      character*255 filfrg,card
      character*10 date,time
c
      integer i,j,k,l,length,number
c
      call swatch(date,time)
c
      do 1 i=1,nparms
      write(filfrg,'(a,a)') dir(1:index(dir,' ')-1),
     + fils(i)(1:index(fils(i),' ')-1)//'.frg '
      length=index(filfrg,' ')-1
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='unknown',err=9999)
      write(lfnfrg,1000) fils(i)(1:index(fils(i),' ')),date,time
 1000 format('$',a,/,'# ',2a10)
      number=0
      do 2 j=1,natm
      ndx2(j)=0
      if(catm(i+1,j)(6:6).ne.'D') then
      number=number+1
      ndx(number)=j
      ndx2(j)=number
      endif
    2 continue
      write(lfnfrg,1002) number,1,1,nzmat,fils(i)(1:index(fils(i),' '))
 1002 format(4i5,/,a)
      do 3 k=1,number
      j=ndx(k)
      write(lfnfrg,1003) k,catm(1,j),catm(i+1,j),
     + (latm(l,j),l=1,5),qatm(i,1,j),qatm(i,2,j)
 1003 format(i5,2a6,5i5,2f12.6)
    3 continue
      do 4 k=1,nbnd
      if(ndx2(lbnd(1,k)).gt.0.and.ndx2(lbnd(2,k)).gt.0) then
      write(lfnfrg,1004) ndx2(lbnd(1,k)),ndx2(lbnd(2,k))
 1004 format(2i5)
      endif
    4 continue
      close(unit=lfnfrg,status='keep')
    1 continue
 
c
  999 continue
c
      argos_prep_master=.true.
      return
c
 9999 continue
      argos_prep_master=.false.
      return
      end
c $Id$
