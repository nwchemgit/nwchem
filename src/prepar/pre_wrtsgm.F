      logical function pre_wrtsgm(lfnout,lfnsgm,filsgm,nparms,nparmsp,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd,lang,mang,nang,
     + ldih,mdih,ndih,limp,mimp,nimp,izmat,zmat,mzmat,nzmat,ipardef)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfnout,lfnsgm,nparms,nparmsp,ipardef
      integer matm,natm,mbnd,nbnd,mang,nang,mdih,ndih,mimp,nimp
      integer latm(6,matm),lbnd(2,mbnd),lang(3,mang),ldih(4,mdih)
      integer limp(4,mimp)
      character*6 catm(nparmsp,matm)
      real*8 qatm(nparms,2,matm)
      integer mzmat,nzmat
      integer izmat(mzmat,4)
      real*8 zmat(mzmat,3)
      character*255 filsgm
      character*10 date,time
c
      integer i,j,length
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='new',err=9999)
c
      call swatch(date,time)
c
      write(lfnsgm,2000) date,time
 2000 format('# This is an automatically generated segment file',/,
     + '# ',2a10)
c
      write(lfnsgm,2001) 4.6
 2001 format(f12.6)
c
      write(lfnsgm,2002) natm,nbnd,nang,ndih,nimp,nzmat,nparms,ipardef
 2002 format(8i5)
c
      write(lfnsgm,2003) (0.0,j=1,nparms)
 2003 format(f12.6)
c
      do 1 i=1,natm
      write(lfnsgm,2004) i,catm(1,i),(latm(j,i),j=1,5)
 2004 format(i5,a6,5i5)
      do 201 j=1,nparms
      write(lfnsgm,2005) catm(1+j,i),qatm(j,1,i),qatm(j,2,i)
 2005 format(5x,a6,2f12.6)
  201 continue
    1 continue
c
      do 2 i=1,nbnd
      write(lfnsgm,2006) i,(lbnd(j,i),j=1,2),0,0
 2006 format(5i5)
      do 202 j=1,nparms
      write(lfnsgm,2007) 0.0d0,0.0d0
 2007 format(f12.6,e12.5)
  202 continue
    2 continue
c
      do 3 i=1,nang
      write(lfnsgm,2008) i,(lang(j,i),j=1,3),0,0
 2008 format(6i5)
      do 203 j=1,nparms
      write(lfnsgm,2009) 0.0d0,0.0d0
 2009 format(f10.6,e12.5)
  203 continue
    3 continue
c
      do 4 i=1,ndih
      write(lfnsgm,2010) i,(ldih(j,i),j=1,4),0,0
 2010 format(7i5)
      do 204 j=1,nparms
      write(lfnsgm,2011) 0,0.0d0,0.0d0
 2011 format(i3,f10.6,e12.5)
  204 continue
    4 continue
c
      do 5 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      write(lfnsgm,2012) i,limp(2,i),limp(3,i),limp(1,i),limp(4,i),
     + 0,0
 2012 format(7i5)
      do 205 j=1,nparms
      write(lfnsgm,2013) 0,0.0d0,0.0d0
 2013 format(i3,f10.6,e12.5)
  205 continue
      else
      write(lfnsgm,2012) i,(limp(j,i),j=1,4),0,0
      do 305 j=1,nparms
      write(lfnsgm,2013) 0,0.0d0,0.0d0
  305 continue
      endif
    5 continue
c
      do 6 i=1,nzmat
      write(lfnsgm,2014) i,(izmat(i,j),j=1,4),(zmat(i,j),j=1,3)
 2014 format(5i5,3f12.6)
    6 continue
c
      close(unit=lfnsgm)
c
      if(util_print('sequence',print_medium)) then
      write(lfnout,1000) filsgm(1:length)
 1000 format(' Created segment',t40,a)
      endif
c
      pre_wrtsgm=.true.
      return
c
 9999 continue
      pre_wrtsgm=.false.
      return
      end
