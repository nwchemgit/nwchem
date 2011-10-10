      subroutine argos_space_listbb(ibownr,boxsiz,lbbl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer ibownr(maxbox,3),lbbl(mbbl,mbb2)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer i,j,jbox,jnode,ltemp
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
c
c     Construction of the box-box pairlist
c
      nbbl=0
      do 1 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 2 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 3 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 4 jbx=0,nbx-1
      kbx=jbx-ibx
      jpx=ibownr(jbx+1,1)
c
      if(ibx.le.jbx) then
      ilx=ibx
      ihx=jbx
      else
      ilx=jbx
      ihx=ibx
      endif
c
      dx=zero
      if(ibx.ne.jbx) then
      dx=boxsiz(ihx,1)-boxsiz(ilx+1,1)
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=boxsiz(ilx,1)
      if(ihx.lt.nbx-1) dxtmp=dxtmp-boxsiz(ihx+1,1)+box(1)
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 5 jby=0,nby-1
      kby=jby-iby
      jpy=ibownr(jby+1,2)
c
      if(iby.le.jby) then
      ily=iby
      ihy=jby
      else
      ily=jby
      ihy=iby
      endif
c
      dy=zero
      if(iby.ne.jby) then
      dy=boxsiz(ihy,2)-boxsiz(ily+1,2)
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=boxsiz(ily,2)
      if(ihy.lt.nby-1) dytmp=dytmp-boxsiz(ihy+1,2)+box(2)
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 6 jbz=0,nbz-1
      kbz=jbz-ibz
      jpz=ibownr(jbz+1,3)
c
      if(ibz.le.jbz) then
      ilz=ibz
      ihz=jbz
      else
      ilz=jbz
      ihz=ibz
      endif
c
      dz=zero
      if(ibz.ne.jbz) then
      dz=boxsiz(ihz,3)-boxsiz(ilz+1,3)
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=boxsiz(ilz,3)
      if(ihz.lt.nbz-1) dztmp=dztmp-boxsiz(ihz+1,3)+box(3)
      if(dztmp.lt.dz) dz=dztmp
      if(kbz.gt.0.and.kbz.gt.iabs(kbz-nbz)) kbz=kbz-nbz
      if(kbz.lt.0.and.-kbz.gt.iabs(kbz+nbz)) kbz=kbz+nbz
      endif
      endif
c
      jbox=(jbz*nby+jby)*nbx+jbx
      jnode=(jpz*npy+jpy)*npx+jpx
c
c     determine orientation jbox in relation to ibox
c
c     lside is true if
c
c     i: 0  j: 0  k:  +
c     i: 0  j:  + k:-0+
c     i:  + j:-0+ k:-0+
c
      lside=(kbx.eq.0.and.kby.eq.0.and.kbz.ge.0)
     + .or.(kbx.eq.0.and.kby.gt.0) .or. kbx.gt.0
c
c     determine if ibox is identical to jbox
c
c     lsame=kbx.eq.0.and.kby.eq.0.and.kbz.eq.0
c
c     determine if difference in box numbers is even or odd
c
      leven=2*(iabs(ibox-jbox)/2).eq.iabs(ibox-jbox)
c
c     calculate the distance between the two boxes
c
      if(nbxtyp.eq.1) then
      dist2=
     + (vlat(1,1)*dx/box(1)+vlat(1,2)*dy/box(2)+vlat(1,3)*dz/box(3))**2+
     + (vlat(2,1)*dx/box(1)+vlat(2,2)*dy/box(2)+vlat(2,3)*dz/box(3))**2+
     + (vlat(3,1)*dx/box(1)+vlat(3,2)*dy/box(2)+vlat(3,3)*dz/box(3))**2
      else
      dist2=dx*dx+dy*dy+dz*dz
      endif
c
c     keep half of the box pairs
c
c     this test also appears in argos_space_numbb
c     any changes need to be made in both routines
c
      if((inode.eq.jnode.and.ibox.ge.jbox).or. (inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rbbl*rbbl.gt.dist2) then
      nbbl=nbbl+1
      if(nbbl.gt.mbbl) call md_abort('Box-box list too small',mbbl)
      lbbl(nbbl,1)=jnode
      lbbl(nbbl,2)=jbox
      lbbl(nbbl,3)=ibox
      lbbl(nbbl,4)=0
      endif
      endif
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
      npprev=0
c
      nbbloc=0
      do 7 i=1,nbbl-1
      do 8 j=i+1,nbbl
      if((lbbl(i,1).ne.me.and.lbbl(j,1).eq.me).or.
     + (lbbl(i,1).gt.lbbl(j,1).and.lbbl(i,1).ne.me).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).gt.lbbl(j,2)).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).eq.lbbl(j,2).and.
     + lbbl(i,3).gt.lbbl(j,3))) then
      ltemp=lbbl(i,1)
      lbbl(i,1)=lbbl(j,1)
      lbbl(j,1)=ltemp
      ltemp=lbbl(i,2)
      lbbl(i,2)=lbbl(j,2)
      lbbl(j,2)=ltemp
      ltemp=lbbl(i,3)
      lbbl(i,3)=lbbl(j,3)
      lbbl(j,3)=ltemp
      endif
    8 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    7 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      nrempr=0
      if(nbget.ne.0.and.np.gt.1) then
      nrempr=1
      do 9 i=2,nbbl
      if(lbbl(i,2).ne.lbbl(i-1,2)) nrempr=nrempr+1
    9 continue
      call ga_igop(msp_29,nrempr,1,'max')
      nrempr=min(2*nrempr,nrempr+25)
      endif
c
      if(iand(idebug,2).eq.2) then
      write(lfndbg,8000) (i,(lbbl(i,j),j=1,3),i=1,nbbl)
 8000 format('lbbl',/,(4i5))
      call util_flush(lfndbg)
      endif
c
      return
      end
c $Id$
