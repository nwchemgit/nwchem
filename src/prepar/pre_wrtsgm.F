      logical function pre_wrtsgm(lfnout,ffield,lfnsgm,filsgm,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd,lang,mang,nang,
     + ldih,mdih,ndih,limp,mimp,nimp,izmat,zmat,mzmat,nzmat)
c
c $Id: pre_wrtsgm.F,v 1.11 2000-09-15 22:21:54 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnsgm
      integer matm,natm,mbnd,nbnd,mang,nang,mdih,ndih,mimp,nimp
      integer latm(5,matm),lbnd(2,mbnd),lang(3,mang),ldih(4,mdih)
      integer limp(4,mimp)
      character*6 catm(2,matm)
      real*8 qatm(2,matm)
      integer mzmat,nzmat
      integer izmat(mzmat,4)
      real*8 zmat(mzmat,3)
      character*255 filsgm
      character*80 ffield
c
      integer i,j,length
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='new',err=9999)
c
      write(lfnsgm,2000)
 2000 format('# This is an automatically generated segment file',/,
     + '#')
c
      write(lfnsgm,2001) natm,nbnd,nang,ndih,nimp,nzmat
 2001 format(5i5,36x,i5)
c
      do 1 i=1,natm
      write(lfnsgm,2002) i,catm(1,i),(catm(2,i),j=1,3),latm(3,i),
     + latm(4,i),latm(1,i),latm(2,i),(qatm(1,i),qatm(2,i),j=1,3)
 2002 format(i5,a6,3a6,4i4,/,3(2f12.6))
    1 continue
c
      do 2 i=1,nbnd
      write(lfnsgm,2003) i,(lbnd(j,i),j=1,2),0,0,(0.0d0,j=1,6)
 2003 format(i5,4i5,/,3(f12.6,e12.5))
    2 continue
c
      do 3 i=1,nang
      write(lfnsgm,2004) i,(lang(j,i),j=1,3),0,0,(0.0d0,j=1,6)
 2004 format(i5,5i5,/,3(f10.6,e12.5))
    3 continue
c
      do 4 i=1,ndih
      write(lfnsgm,2005) i,(ldih(j,i),j=1,4),0,0,(0,0.0d0,0.0d0,j=1,3)
 2005 format(i5,6i5,/,3(i3,f10.6,e12.5))
    4 continue
c
      do 5 i=1,nimp
      if(ffield(1:5).eq.'amber') then
      write(lfnsgm,2006) i,limp(2,i),limp(3,i),limp(1,i),limp(4,i),
     + 0,0,(0,0.0d0,0.0d0,j=1,3)
      else
      write(lfnsgm,2006) i,(limp(j,i),j=1,4),0,0,(0,0.0d0,0.0d0,j=1,3)
      endif
 2006 format(i5,6i5,/,3(i3,f10.6,e12.5))
    5 continue
c
      do 6 i=1,nzmat
      write(lfnsgm,2007) i,(izmat(i,j),j=1,4),(zmat(i,j),j=1,3)
 2007 format(5i5,3f12.6)
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
