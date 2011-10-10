      subroutine argos_space_numbb(ibownr,boxsiz)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
c
      integer ibownr(maxbox,3)
      real*8 boxsiz(maxbox,3)
      logical lside,leven
c
      integer ibx,iby,ibz,ipx,ipy,ipz,ibox,inode
      integer jbx,jby,jbz,kbx,kby,kbz,jpx,jpy,jpz
      integer ilx,ihx,ily,ihy,ilz,ihz
      integer jbox,jnode,mbblb,i
      real*8 dx,dxtmp,dy,dytmp,dz,dztmp,dist2
      character*255 string
c
      if(iand(idebug,2).eq.2) then
      write(lfndbg,'(a,i6)') 'boxsiz in argos_space_numbb'
      write(lfndbg,'(6f12.6)') (boxsiz(i,1),i=1,nbx)
      write(lfndbg,'(6f12.6)') (boxsiz(i,2),i=1,nby)
      write(lfndbg,'(6f12.6)') (boxsiz(i,3),i=1,nbz)
      call util_flush(lfndbg)
      endif
c
c     determine size of box-box pairlist
c     ----------------------------------
c
      mbbl=0
      mbblb=0
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
      if((inode.eq.jnode.and.ibox.ge.jbox).or.(inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rbbl*rbbl.gt.dist2) mbbl=mbbl+1
      endif
      if(rbbl*rbbl.gt.dist2) mbblb=mbblb+1
 6    continue
 5    continue
 4    continue
      endif
 3    continue
 2    continue
 1    continue
c
      mbbl=mbbl+1
      if(mbblb+1.gt.mbbl) mbbl=mbblb+1
c
      if(np.gt.1) call ga_igop(msp_20,mbbl,1,'max')
      if(me.eq.0) then
      if(util_print('distribution',print_debug)) then
      write(lfnout,2001) mbbl
 2001 format(' Dimension of the cell-cell list is ',i7)
      endif
      endif
c
c     determine size of box-box pairlist assuming minimum cell sizes
c     --------------------------------------------------------------
c
      mbbl=0
      mbblb=0
      do 11 ibx=0,nbx-1
      ipx=ibownr(ibx+1,1)
      do 12 iby=0,nby-1
      ipy=ibownr(iby+1,2)
      do 13 ibz=0,nbz-1
      ipz=ibownr(ibz+1,3)
      ibox=(ibz*nby+iby)*nbx+ibx
      inode=(ipz*npy+ipy)*npx+ipx
      if(inode.eq.me) then
      do 14 jbx=0,nbx-1
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
      dx=dble(ihx-ilx-1)*bxmin
      if(npbtyp.gt.0) then
      dxtmp=zero
      if(ilx.gt.0) dxtmp=dble(ilx)*bxmin
      if(ihx.lt.nbx-1) dxtmp=dxtmp+dble(nbx-ihx-1)*bxmin
      if(dxtmp.lt.dx) dx=dxtmp
      if(kbx.gt.0.and.kbx.gt.iabs(kbx-nbx)) kbx=kbx-nbx
      if(kbx.lt.0.and.-kbx.gt.iabs(kbx+nbx)) kbx=kbx+nbx
      endif
      endif
c
      do 15 jby=0,nby-1
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
      dy=dble(ihy-ily-1)*bymin
      if(npbtyp.gt.0) then
      dytmp=zero
      if(ily.gt.0) dytmp=dble(ily)*bymin
      if(ihy.lt.nby-1) dytmp=dytmp+dble(nby-ihy-1)*bymin
      if(dytmp.lt.dy) dy=dytmp
      if(kby.gt.0.and.kby.gt.iabs(kby-nby)) kby=kby-nby
      if(kby.lt.0.and.-kby.gt.iabs(kby+nby)) kby=kby+nby
      endif
      endif
c
      do 16 jbz=0,nbz-1
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
      dz=dble(ihz-ilz-1)*bzmin
      if(npbtyp.eq.1) then
      dztmp=zero
      if(ilz.gt.0) dztmp=dble(ilz)*bzmin
      if(ihz.lt.nbz-1) dztmp=dztmp-dble(nbz-ihz-1)*bzmin
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
      if((inode.eq.jnode.and.ibox.ge.jbox).or.(inode.ne.jnode.and.
     + ((lside.and.leven).or.(.not.lside.and..not.leven)))) then
c
c     keep only those within maximum cutoff distance
c
      if(rbbl*rbbl.gt.dist2) mbbl=mbbl+1
      endif
      if(rbbl*rbbl.gt.dist2) mbblb=mbblb+1
   16 continue
   15 continue
   14 continue
      endif
   13 continue
   12 continue
   11 continue
c
      mbbl=mbbl+1
      if(mbblb+1.gt.mbbl) mbbl=mbblb+1
c
c
c     allocate memory for the box-box list
c
      if(np.gt.1) call ga_igop(msp_20,mbbl,1,'max')
c
      if(me.eq.0) then
      if(util_print('distribution',print_debug)) then
      write(lfnout,2002) mbbl
 2002 format(' Dimension of the cell-cell list is ',i7)
      endif
      endif
c
      if(mbblp.eq.0) then
      mbbl=max(mbbl,mbbreq)
      if(.not.ma_push_get(mt_int,mbb2*mbbl,'bb',l_bb,i_bb))
     + call md_abort('Failed to allocate memory for bb',0)
      mbblp=mbbl
      else
      if(mbbl.gt.mbblp) then
      write(string,1111) mbblp,mbbl,mbbreq
 1111 format('error: lbbl increase from ',i6,' to ',i6,'(',i6,')')
      call md_abort(string,me)
c      call md_abort('lbbl increased beyond allocated memory',me)
      endif
      endif
c
      return
      end
c $Id$
