      logical function pre_master(lfnfrg,lfnout,nparms,nparmsp,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd,izmat,zmat,mzmat,nzmat,
     + dir,fils,ndx)
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfnfrg,lfnout,matm,natm,mbnd,nbnd,mzmat,nzmat
      integer nparms,nparmsp
      integer latm(6,matm),lbnd(2,mbnd),izmat(mzmat,4),ndx(matm)
      character*6 catm(nparmsp,matm)
      real*8 qatm(nparms,2,matm),zmat(mzmat,3)
      character*255 dir,fils(nparmsp)
c
      character*255 filfrg,card
c
      integer i,j,k,length,number
c
      do 1 i=1,nparms
      write(filfrg,'(a,a)') dir(1:index(dir,' ')-1),
     + fils(i)(1:index(fils(i),' ')-1)//'.frg '
      length=index(filfrg,' ')-1
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='unknown',err=9999)
      write(lfnfrg,1000) fils(i)(1:index(fils(i),' '))
 1000 format('$',a)
      number=0
      do 2 j=1,natm
      if(catm(i+1,j)(6:6).ne.'D') then
      number=number+1
      ndx(number)=j
      endif
    2 continue
      write(lfnfrg,1002) number,1,1,nzmat
 1002 format(4i5)
      do 3 k=1,number
      j=ndx(k)
      write(lfnfrg,1003) k,catm(1,j),catm(i+1,j),
     + latm(3,j),latm(4,j),latm(1,j),latm(2,j),latm(6,j),
     + qatm(i,1,j),qatm(i,2,j)
 1003 format(i5,2a6,5i5,2f12.6)
    3 continue
      close(unit=lfnfrg,status='keep')
    1 continue
c
  999 continue
c
      pre_master=.true.
      return
c
 9999 continue
      pre_master=.false.
      return
      end
