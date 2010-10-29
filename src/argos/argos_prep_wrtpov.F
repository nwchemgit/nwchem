      logical function argos_prep_wrtpov(lfnout,lfnpov,filpov,iopt,box,
     + num,amass,mat,
     + csa,isat,isgm,xs,vs,msa,nsa,cwa,iwat,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc,slvnam,nxrep,nyrep,nzrep,drep,msb,nsb,idsb,
     + rdist,scale,cpk)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      character*2 argos_prep_atnam
      real*8 atom_radius
      integer argos_prep_atnum
      character*255 atom_color
      external argos_prep_atnam,atom_radius,argos_prep_atnum,atom_color
c
      integer lfnout,lfnpov,iopt,msa,nsa,mwm,nwm,mwa,nwa,mat,mwmc,nwmc
      character*255 filpov
      integer num(mat),isat(msa),isgm(msa),iwat(mwa)
      character*16 cwa(mwa),csa(msa)
      character*3 slvnam
      real*8 amass(mat),xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),box(3)
      integer nxrep,nyrep,nzrep
      real*8 rdist,drep(3),xh(3)
      integer msb,nsb
      integer idsb(2,msb)
      real*8 scale,cpk
c
      integer length,ib,jb
      character*255 color
c
      integer i,j,k
c
      real*8 xbmax,xbmin,ybmax,ybmin,zbmax,zbmin
c
      xbmax=xs(1,1)+0.2d0
      xbmin=xs(1,1)-0.2d0
      ybmax=xs(2,1)+0.2d0
      ybmin=xs(2,1)-0.2d0
      zbmax=xs(3,1)+0.2d0
      zbmin=xs(3,1)-0.2d0
      do 1 i=2,nsa
      xbmax=max(xbmax,xs(1,i))
      xbmin=min(xbmin,xs(1,i))
      ybmax=max(ybmax,xs(2,i))
      ybmin=min(ybmin,xs(2,i))
      zbmax=max(zbmax,xs(3,i))
      zbmin=min(zbmin,xs(3,i))
    1 continue
      if(scale.le.0.0d0) then
      scale=2.0d0/max(xbmax-xbmin,ybmax-ybmin,zbmax-zbmin)
      xbmax=scale*xbmax
      xbmin=scale*xbmin
      ybmax=scale*ybmax
      ybmin=scale*ybmin
      zbmax=scale*zbmax
      zbmin=scale*zbmin
      endif
      call povinc(lfnpov,xbmin,xbmax,ybmin,ybmax,zbmin,zbmax)
c
c     open POV file
c
      length=index(filpov,' ')-1
      open(unit=lfnpov,file=filpov(1:length),form='formatted',
     + status='unknown',err=9999)
      write(lfnpov,8000)
 8000 format('#include "camera.inc"',/,'#include "colors.inc"',/,
     + '#include "plane.inc"')
c
      do 2 i=1,nsa
      k=argos_prep_atnum(csa(i)(11:12))
      color=atom_color(k)
      if(k.gt.0) then
      write(lfnpov,8003) isgm(i),csa(i)
 8003 format('// ',i5,':',a16)
      if(cpk.gt.0.0d0) then
      write(lfnpov,8001) (scale*xs(j,i),j=1,3),scale*cpk*atom_radius(k),
     + color(1:index(color,' ')-1)
      else
      write(lfnpov,8001) (scale*xs(j,i),j=1,3),
     + scale*abs(cpk)*atom_radius(1),
     + color(1:index(color,' ')-1)
      endif
 8001 format('sphere {',/,'  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + 3x,f12.6,/,
     + '  texture { pigment { ',a,' }',/,
     + '            finish  { ambient 0.16 ',/,
     + '                      diffuse 0.48 ',/,
     + '                      phong 1.25 ',/,
     + '                      phong_size 200 } }',/,'}')
      endif
    2 continue
      if(cpk.le.0.0d0) then
      do 3 i=1,nsb
      ib=idsb(1,i)
      jb=idsb(2,i)
      xh(1)=0.5d0*(xs(1,ib)+xs(1,jb))
      xh(2)=0.5d0*(xs(2,ib)+xs(2,jb))
      xh(3)=0.5d0*(xs(3,ib)+xs(3,jb))
      k=argos_prep_atnum(csa(ib)(11:12))
      color=atom_color(k)
      write(lfnpov,8002) (scale*xs(j,ib),j=1,3),
     + (scale*xh(j),j=1,3),scale*abs(cpk)*atom_radius(1),
     + color(1:index(color,' ')-1)
 8002 format('cylinder {',/,'  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + '  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + 3x,f12.6,/,'  open',/,'  texture { pigment { ',a,' } }',/,'}')
      k=argos_prep_atnum(csa(jb)(11:12))
      color=atom_color(k)
      write(lfnpov,8002) (scale*xs(j,jb),j=1,3),
     + (scale*xh(j),j=1,3),scale*abs(cpk)*atom_radius(1),
     + color(1:index(color,' ')-1)
    3 continue
      endif
      close(unit=lfnpov)
c
      if(util_print('files',print_default)) then
      write(lfnout,2000) filpov(1:length)
 2000 format(' Created pov',t40,a,/)
      endif
c
      argos_prep_wrtpov=.true.
      return
c
 9999 continue
      argos_prep_wrtpov=.false.
      return
      end
