      subroutine pre_wrtpmf(lfnpmf,filpmf,numpmf,maxpmf,ipmf,rpmf,
     + maxgrp,maxatm,igroup,lgroup)
c
      implicit none
c
      integer lfnpmf
      character*(*) filpmf
c
      integer numpmf,maxpmf,maxgrp,maxatm
      integer ipmf(maxpmf,5),rpmf(maxpmf,4)
      integer igroup(maxgrp,maxatm),lgroup(maxgrp)
c
      integer i,j,k,m
c
      m=0
      k=0
      do 5 i=1,numpmf
      do 6 j=2,ipmf(i,1)+2
      m=max(m,lgroup(ipmf(i,j)))
    6 continue
      if(ipmf(i,1).eq.5) k=k+lgroup(ipmf(i,2))-1
    5 continue
c
c     write pmf
c
      if(numpmf.gt.0) then
      open(unit=lfnpmf,file=filpmf(1:index(filpmf,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(lfnpmf,1000) numpmf+k,m
 1000 format(2i5)
      do 1 i=1,numpmf
      if(ipmf(i,1).eq.0) call errquit('Unknown pmf type',ipmf(i,1))
c
c     distance pmf
c
      if(ipmf(i,1).eq.1) then
      write(lfnpmf,1001) 2,(rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,3)
 1001 format(i5,2f12.6,2e12.5,4i5)
      do 2 k=2,3
      if(lgroup(ipmf(i,k)).le.0)
     + call errquit('No atoms selected in group',ipmf(i,k))
      write(lfnpmf,1002) (igroup(ipmf(i,k),j),j=1,lgroup(ipmf(i,k)))
 1002 format(10i7)
    2 continue
      endif
c
c     angle pmf
c
      if(ipmf(i,1).eq.2) then
      write(lfnpmf,1001) 3,(rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,4)
      do 3 k=2,4
      if(lgroup(ipmf(i,k)).le.0)
     + call errquit('No atoms selected in group',ipmf(i,k))
      write(lfnpmf,1002) (igroup(ipmf(i,k),j),j=1,lgroup(ipmf(i,k)))
    3 continue
      endif
c
c     torsion pmf
c
      if(ipmf(i,1).eq.3) then
      write(lfnpmf,1001) 4,(rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,5)
      do 4 k=2,5
      if(lgroup(ipmf(i,k)).le.0)
     + call errquit('No atoms selected in group',ipmf(i,k))
      write(lfnpmf,1002) (igroup(ipmf(i,k),j),j=1,lgroup(ipmf(i,k)))
    4 continue
      endif
c
c     alignment pmf
c
      if(ipmf(i,1).eq.5) then
      do 7 k=1,lgroup(ipmf(i,2))
      write(lfnpmf,1001) 5,(rpmf(i,j),j=1,4),1,lgroup(ipmf(i,2))
      if(lgroup(ipmf(i,2)).le.0)
     + call errquit('No atoms selected in group',ipmf(i,2))
      write(lfnpmf,1002) igroup(ipmf(i,2),k)
      write(lfnpmf,1002) (igroup(ipmf(i,2),j),j=1,lgroup(ipmf(i,2)))
    7 continue
      endif
c
    1 continue
c
      close(unit=lfnpmf)
      endif
c
      return
c
 9999 continue
      call errquit('Could not open filpmf',0)
      return
      end
