      subroutine argos_cafe_plist(nums,xs,fs,isgan,isatt,
     + ipmf,jpmf,rpmf,xpmf,wpmf,upmf,lfnpmf,lrec,iprpmf)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 argos_cafe_atmass
      external argos_cafe_atmass
c
      integer nums,lfnpmf,iprpmf
      real*8 xs(msa,3),fs(msa,3)
      integer isgan(msa),isatt(msa)
      integer ipmf(numpmf,8),jpmf(numpmf,4,npmfa)
      real*8 rpmf(3,6,numpmf),xpmf(4,3,numpmf),wpmf(4,numpmf)
      real*8 et(2),ft(4,3),upmf(numpmf)
      integer ip(5)
      real*8 rv(2),chg(2)
      integer it(2),ign(4),iml(4),idt(4),ich(4)
      logical lrec
c
      integer i,j,k,l,m
      real*8 wm,wmtot
      logical skip
c
      integer lseq(2),isgm(2)
c
      call timer_start(20)
      epmf=zero
      skip=.true.
c
      chg(1)=zero
c
      do 1 i=1,numpmf
      do 2 j=1,4
      xpmf(j,1,i)=zero
      xpmf(j,2,i)=zero
      xpmf(j,3,i)=zero
      wpmf(j,i)=zero
    2 continue
      upmf(i)=zero
      if(ipmf(i,6).eq.2.and.lpmfc) goto 1
c
      m=0
      if(ipmf(i,1).eq.1) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).eq.5) m=2
      if(ipmf(i,1).eq.6) m=1
      if(ipmf(i,1).eq.8) m=1
      if(ipmf(i,1).eq.9) m=2
      if(m.eq.0) call md_abort('Error in pmf xx',ipmf(i,1))
      skip=.false.
      do 3 j=1,m
      wmtot=0.0d0
      do 4 k=1,ipmf(i,j+1)
      do 5 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=argos_cafe_atmass(l,isatt)
      wmtot=wmtot+wm
      xpmf(j,1,i)=xpmf(j,1,i)+wm*xs(l,1)
      xpmf(j,2,i)=xpmf(j,2,i)+wm*xs(l,2)
      xpmf(j,3,i)=xpmf(j,3,i)+wm*xs(l,3)
      wpmf(j,i)=wpmf(j,i)+wm
      endif
    5 continue
    4 continue
    3 continue
    1 continue
      call timer_stop(20)
c
      if(skip) return
c
      call timer_start(21)
      call ga_dgop(mcf_59,wpmf,4*numpmf,'+')
      call ga_dgop(mcf_73,xpmf,12*numpmf,'+')
      call timer_stop(21)
c
      call timer_start(22)
      do 6 i=1,numpmf
      if(ipmf(i,6).eq.2.and.lpmfc) goto 6
c
      m=0
      if(ipmf(i,1).eq.1) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).eq.5) m=2
      if(ipmf(i,1).eq.6) m=1
      if(ipmf(i,1).eq.8) m=1
      if(ipmf(i,1).eq.9) m=2
      do 222 j=1,m
      wmtot=one/wpmf(j,i)
      wpmf(j,i)=wmtot
      xpmf(j,1,i)=xpmf(j,1,i)*wmtot
      xpmf(j,2,i)=xpmf(j,2,i)*wmtot
      xpmf(j,3,i)=xpmf(j,3,i)*wmtot
  222 continue
c
c     5: alignment in z-direction does not change z coordinates
c
      if(ipmf(i,1).eq.5) xpmf(2,3,i)=xpmf(1,3,i)
c
c     9: z-distance pmf does not change x or y coordinates
c
      if(ipmf(i,1).eq.9) then
      xpmf(2,1,i)=xpmf(1,1,i)
      xpmf(2,2,i)=xpmf(1,2,i)
      endif
c
c     10: z-axis COG pmf does not change z coordinates
c
      if(ipmf(i,1).eq.10) then
      xpmf(2,1,i)=zero
      xpmf(2,2,i)=zero
      xpmf(2,3,i)=xpmf(1,3,i)
      endif
c
      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.5.or.ipmf(i,1).eq.9.or.
     + ipmf(i,1).eq.10) then
      do 7 j=1,2
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 8 k=1,3
      ft(j,k)=zero
    8 continue
    7 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=0
      call argos_cafe_fsb(1,it,1,3,ip,rpmf(1,1,i),rv,2,4,
     + ign,isgm,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.false.,.false.,.false.,me.eq.0,
     + dbl_mb(i_sti),lseq)
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0.and.(iprpmf.lt.0.or.iprpmf.eq.ipmf(i,1)))
     + write(-lfnpmf,1000) i,rv
 1000 format(i5,f12.6,e12.5)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      epmf=epmf+rv(2)
      upmf(ipmf(i,8))=upmf(ipmf(i,8))+rv(2)
      endif
      do 9 j=1,2
      do 10 k=1,ipmf(i,j+1)
      do 11 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=argos_cafe_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   11 continue
   10 continue
    9 continue
      endif
c
      if(ipmf(i,1).eq.2) then
      do 12 j=1,3
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 13 k=1,3
      ft(j,k)=zero
   13 continue
   12 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=3
      ip(4)=0
      call argos_cafe_fsh(1,it,1,3,ip,rpmf(1,1,i),rv,3,4,
     + ign,isgm,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,me.eq.0,
     + dbl_mb(i_sti),lseq)
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Angle:',rv,(ft(3,j),j=1,3)
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0) write(-lfnpmf,1000) i,rv
      epmf=epmf+rv(2)
      upmf(ipmf(i,8))=upmf(ipmf(i,8))+rv(2)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      endif
      do 14 j=1,3
      do 15 k=1,ipmf(i,j+1)
      do 16 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=argos_cafe_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   16 continue
   15 continue
   14 continue
      endif
c
      if(ipmf(i,1).eq.3.or.ipmf(i,1).eq.4) then
      do 17 j=1,4
      ign(j)=j
      iml(j)=1
      idt(j)=0
      ich(j)=1
      do 18 k=1,3
      ft(j,k)=zero
   18 continue
   17 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=3
      ip(4)=4
      ip(5)=0
      if(ipmf(i,1).eq.3) then
      call argos_cafe_fsd(1,it,1,3,ip,rpmf(1,1,i),rv,4,4,ign,
     + isgm,iml,idt,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,me.eq.0,
     + dbl_mb(i_sti),lseq)
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Dih:',rv,(ft(3,j),j=1,3)
      endif
      if(ipmf(i,1).eq.4) then
      call argos_cafe_fso(1,it,1,3,ip,rpmf(1,1,i),rv,4,4,ign,
     + isgm,iml,idt,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,me.eq.0,
     + dbl_mb(i_sti),lseq)
c      write(*,'(a,f12.6,e12.5,3e12.5)') 'Imp:',rv,(ft(3,j),j=1,3)
      endif
      if(me.eq.0) then
      if(lrec.and.lfnpmf.lt.0) write(-lfnpmf,1000) i,rv
      epmf=epmf+rv(2)
      upmf(ipmf(i,8))=upmf(ipmf(i,8))+rv(2)
      if(ipmf(i,6).eq.1) ubias=ubias+rv(2)
      endif
      do 19 j=1,4
      do 20 k=1,ipmf(i,j+1)
      do 21 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      wm=argos_cafe_atmass(l,isatt)
      fs(l,1)=fs(l,1)+ft(j,1)*wm*wpmf(j,i)
      fs(l,2)=fs(l,2)+ft(j,2)*wm*wpmf(j,i)
      fs(l,3)=fs(l,3)+ft(j,3)*wm*wpmf(j,i)
      endif
   21 continue
   20 continue
   19 continue
      endif
c
c     planar pmf
c
      if(ipmf(i,1).eq.6.or.ipmf(i,1).eq.8) then
c
c     for planar pmf copy x and y coordinates of center of mass
c
      if(ipmf(i,1).eq.6) then
      xpmf(2,1,i)=xpmf(1,1,i)
      xpmf(2,2,i)=xpmf(1,2,i)
      endif
c
c     for z-axis pmf set x and y coordinates to zero, and
c                        z coordinate to center of mass
c
      if(ipmf(i,1).eq.8) then
      xpmf(2,1,i)=zero
      xpmf(2,2,i)=zero
      xpmf(2,3,i)=xpmf(1,3,i)
      endif
c
c     loop over all atoms on this processor
c
      do 22 j=1,nums
c
c     loop over all atoms in the selected group
c
      do 23 k=1,ipmf(i,2)
      if(isgan(j).eq.jpmf(i,1,k)) then
c
c     copy z coordinate of the atom in xpmf
c
      if(ipmf(i,1).eq.6) xpmf(2,3,i)=xs(j,3)
      if(ipmf(i,1).eq.8) then
      xpmf(2,1,i)=xs(j,1)
      xpmf(2,2,i)=xs(j,2)
      endif
c
c     initialize arrays for fsb
c
      do 24 l=1,2
      ign(l)=l
      iml(l)=1
      idt(l)=0
      ich(l)=1
      do 25 m=1,3
      ft(l,m)=zero
   25 continue
   24 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=0
      call argos_cafe_fsb(1,it,1,3,ip,rpmf(1,1,i),rv,2,4,
     + ign,isgm,iml,idt,ich,chg,
     + xpmf(1,1,i),ft,et,.true.,.true.,.false.,me.eq.0,
     + dbl_mb(i_sti),lseq)
c
c     sum ubias on every processor since this is not done in replicated mode
c
      ubias=ubias+rv(2)
      epmf=epmf+rv(2)
      upmf(ipmf(i,8))=upmf(ipmf(i,8))+rv(2)
c
c     sum the forces
c
      wm=argos_cafe_atmass(j,isatt)
      fs(j,3)=fs(j,3)+ft(2,3)*wm*wpmf(1,i)
c
c     note that the reaction force is not used since they sum to zero.
c
      endif
   23 continue
   22 continue
      endif
c
    6 continue
c     if(me.eq.0) write(*,'(2i5,e12.5)')
c    +  (i,ipmf(i,8),upmf(ipmf(i,8)),i=1,numpmf)
      call timer_stop(22)
c
      return
      end
c $Id$
