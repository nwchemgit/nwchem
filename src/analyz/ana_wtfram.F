      subroutine ana_wtfram(iunit,fmt,x,sgmnam,w,isel,tag,val,
     + iwsel,ndxw)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3),val(nsa,2)
      integer isel(msa),iwsel(mwm),ndxw(nwm)
      character*24 tag(msa,2)
      integer iunit
      character*3 fmt
      character*16 sgmnam(msa)
c
      character*10 today,now
      integer i,j,k,l
      character*2 elemnt
      character*2 atom
c
      character*1 cdummy
      character*7 string
      integer nsb,ib,jb,nss,nwb
c
      character*255 fname
      character*255 color
      real*8 xh(3)
c
      real*8 xbmax,xbmin,ybmax,ybmin,zbmax,zbmin,factor
      integer*4 numato
c
      real*8 ana_atrad
      integer pre_atnum
      character*255 atom_color
      external ana_atrad,pre_atnum,atom_color
c
      real*8 rv(3),rw(3),rotang,rx(3),ry(3),rxi(3),rxj(3),twopi
c
      logical logw,logs
c
      logs=.true.
      logw=.not.lsonly
c
      twopi=8.0d0*tan(1.0d0)
c
      if(me.eq.0) then
c
      if(iunit.eq.lfncop) then
      if(ncop.gt.mcopf.and.mcopf.gt.0) then
      ncop=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      icopf=icopf+1
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),tag,isel,.false.)
      endif
      endif
c
      if(iunit.eq.lfnsup) then
      if(nsup.gt.msupf.and.msupf.gt.0) then
      nsup=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      isupf=isupf+1
      call ana_wthdr(lfnsup,fmtsup,byte_mb(i_snam),tag,isel,.false.)
      endif
      endif
c
      if(fmt.eq.'trj') then
      write(iunit,1000)
 1000 format('frame')
      write(iunit,1001) time,temp,pres,datum,tijd
 1001 format(2f12.6,e12.5,2a10)
c      write(*,'(a,2f12.6,e12.5)') 'frame ',time,temp,pres
      write(iunit,1002) box(1),0.0,0.0,0.0,box(2),0.0,
     + 0.0,0.0,box(3)
 1002 format(3f12.6,/,3f12.6,/,3f12.6)
      write(iunit,1003) lxw.and..not.lsonly,.false.,.false.,.false.,
     + lsx,.false.,.false.,.false.,nwsel,nwa,nsa
 1003 format(8l1,3i10)
      if(lxw.and..not.lsonly) then
      do 1 i=1,nwm
      if(iwsel(i).ne.0) then
      do 2 j=1,nwa
      write(iunit,1004) (w(i,j,k),k=1,3)
 1004 format(3f8.3)
    2 continue
      endif
    1 continue
      endif
      if(lsx) then
      do 3 i=1,nsa
      write(iunit,1005) (x(i,j),j=1,3)
 1005 format(3f8.3)
    3 continue
      endif
      ncop=ncop+1
      endif
c
      if(fmt.eq.'arc') then
      call swatch(today,now)
      write(iunit,2000) time,today,now
 2000 format('title ',t67,f10.3,/,'!DATE ',2a10)
      nss=0
      if(lsx) then
      do 4 i=1,nsa
      if(.not.lesd.or.isel(i).ne.0) then
      string=sgmnam(i)(13:16)//'   '
      if(sgmnam(i)(13:13).eq.' ') string=sgmnam(i)(14:16)//'    '
      if(sgmnam(i)(14:14).eq.' ') string=sgmnam(i)(15:16)//'     '
      if(sgmnam(i)(15:15).eq.' ') string=sgmnam(i)(16:16)//'      '
      if(sgmnam(i)(6:6).eq.' '.or.
     + sgmnam(i)(6:6).eq.'2'.or.sgmnam(i)(6:6).eq.'3'.or.
     + sgmnam(i)(6:6).eq.'4'.or.sgmnam(i)(6:6).eq.'1') then
      write(iunit,2001) sgmnam(i)(6:9),(1.0d1*x(i,j),j=1,3),
     + sgmnam(i)(1:4),string,'?',sgmnam(i)(7:7)//' ',0.0
      else
      write(iunit,2001) sgmnam(i)(6:9),(1.0d1*x(i,j),j=1,3),
     + sgmnam(i)(1:4),string,'?',sgmnam(i)(6:7),0.0
      endif
 2001 format(a4,1x,3f15.9,1x,a4,1x,a7,a1,7x,a2,f7.3)
      endif
    4 continue
      write(iunit,2002)
      read(sgmnam(nsa)(13:16),'(i4)') nss
      endif
      if(lxw.and..not.lsonly) then
      do 44 i=1,nwm
      if(iwsel(i).ne.0) then
      write(string,'(i7)') i+nss
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      if(string(1:1).eq.' ') string=string(2:7)//' '
      do 45 j=1,nwa
      if(wnam(j)(6:6).eq.' '.or.
     + wnam(j)(6:6).eq.'2'.or.wnam(j)(6:6).eq.'3'.or.
     + wnam(j)(6:6).eq.'4'.or.wnam(j)(6:6).eq.'1') then
      write(iunit,2001) wnam(j)(6:9),(1.0d1*w(i,j,k),k=1,3),
     + wnam(j)(1:4),string,'?',wnam(j)(7:7)//' ',0.0
      else
      write(iunit,2001) wnam(j)(6:9),(1.0d1*w(i,j,k),k=1,3),
     + wnam(j)(1:4),string,'?',wnam(j)(6:7),0.0
      endif
   45 continue
      write(iunit,2002)
      endif
   44 continue
      endif
      write(iunit,2002)
 2002 format('end')
      endif
c
      if(fmt.eq.'xyz') then
      write(iunit,*) nsa+nwa*nwm
      write(iunit,*)
      do 22 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,*) elemnt,(1.0d1*x(i,j),j=1,3)
   22 continue
      do 23 i=1,nwm
      if(iwsel(i).ne.0) then
      do 24 j=1,nwa
      atom=wnam(j)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(iunit,2501) atom,(1.0d1*w(i,j,k),k=1,3)
 2501 format(a2,3f12.6)
   24 continue
      endif
   23 continue
      endif
c
      if(fmt.eq.'amb') then
      if(lxw.and..not.lsonly) then
      write(iunit,3000) ((1.0d1*(x(i,j)+0.5d0*box(j)),j=1,3),i=1,nsa),
     + (((1.0d1*(w(ndxw(i),j,k)+0.5d0*box(k)),k=1,3),j=1,nwa),i=1,nwsel)
      else
      write(iunit,3000) ((1.0d1*(x(i,j)+0.5d0*box(j)),j=1,3),i=1,nsa)
      endif
c      write(iunit,3000) (1.0d1*box(j),j=1,3)
 3000 format(10f8.3)
      endif
c
      if(fmt.eq.'crd') then
      if(lxw.and..not.lsonly) then
      write(iunit,3000) ((1.0d1*(x(i,j)),j=1,3),i=1,nsa),
     + (((1.0d1*(w(ndxw(i),j,k)),k=1,3),j=1,nwa),i=1,nwsel)
      else
      write(iunit,3000) ((1.0d1*(x(i,j)),j=1,3),i=1,nsa)
      endif
      write(iunit,3000) (1.0d1*box(j),j=1,3)
      endif
c
      if(fmt.eq.'bam') then
      if(lxw.and..not.lsonly) then
      write(iunit) ((real(1.0d1*(x(i,j)+0.5d0*box(j))),j=1,3),i=1,nsa),
     + (((1.0d1*(w(i,j,k)+0.5d0*box(k)),k=1,3),j=1,nwa),i=1,nwm),
     + (real(1.0d1*box(j)),j=1,3)
      else
      write(iunit) ((real(1.0d1*(x(i,j)+0.5d0*box(j))),j=1,3),i=1,nsa),
     + (real(1.0d1*box(j)),j=1,3)
      endif
      endif
c
      if(fmt.eq.'dcd') then
      if(.not.ldcd) then
      numato=nsa
      if(lxw.and..not.lsonly) numato=nsa+nwa*nwm
      write(iunit) numato
      ldcd=.true.
      endif
      do 72 k=1,3
      if(lxw.and..not.lsonly) then
      write(iunit) (1.0d1*(x(i,k)+0.5d0*box(k)),i=1,nsa),
     + ((1.0d1*(w(i,j,k)+0.5d0*box(k)),j=1,nwa),i=1,nwm)
      else
      write(iunit) (1.0d1*(x(i,k)+0.5d0*box(k)),i=1,nsa)
      endif
   72 continue
      endif
c
      if(fmt.eq.'eci') then
      write(iunit,5000) time
 5000 format('Time ',f12.4)
      write(iunit,5001) nsa
 5001 format('Geometry{',/,i7,/,' ')
      do 5 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,5002) elemnt,(1.0d1*x(i,j),j=1,3)
 5002 format(a2,3f12.4)
    5 continue
      write(iunit,5003)
 5003 format('}')
      endif
c
      if(fmt.eq.'mvm') then
      write(iunit,6001) nsa
 6001 format('Geometry{',/,'num_atoms:',i7,/,'atom_list:')
      do 6 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,6002) elemnt,(1.0d1*x(i,j),j=1,3)
 6002 format(a2,3f12.4)
    6 continue
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
 4401 format(4i10)
      do 7 i=1,nwa+nsa+nwb
      read(44,4402) cdummy
 4402 format(a1)
    7 continue
      write(iunit,6003) nsb
 6003 format('num_bonds:      ',i7,/,'bond_list:')
      do 8 i=1,nsb
      read(44,6004) ib,jb
      write(iunit,6004) ib,jb
 6004 format(2i8)
    8 continue
      close(unit=44)
      write(iunit,6005)
 6005 format('}')
      endif
c
      if(fmt.eq.'frm') then
      if(lsx) then
      write(iunit,7001)
 7001 format(/,'Geometry {')
      do 9 i=1,nsa
      if(isel(i).ne.0) write(iunit,7002) (1.0d1*x(i,j),j=1,3)
 7002 format(3f12.4)
    9 continue
      write(iunit,7003)
 7003 format('}')
      endif
      endif
c
      if(fmt.eq.'pdb') then
      write(fname,'(a,i5.5,a)') root(1:index(root,'.')-1),indx,'.pdb '
      write(*,'(a,a)') 'Writing ',fname
      call ana_wrtpdb(lfnchg,fname,sgmnam,isel,x,iwsel,w,
     + logs,logw,ndxw)
      endif
c
      if(fmt.eq.'pov') then
      xbmax=x(1,1)
      xbmin=x(1,1)
      ybmax=x(1,2)
      ybmin=x(1,2)
      zbmax=x(1,3)
      zbmin=x(1,3)
      do 234 i=2,nsa
      xbmax=max(xbmax,x(i,1))
      xbmin=min(xbmin,x(i,1))
      ybmax=max(ybmax,x(i,2))
      ybmin=min(ybmin,x(i,2))
      zbmax=max(zbmax,x(i,3))
      zbmin=min(zbmin,x(i,3))
  234 continue
      if(lxw.and..not.lsonly) then
      do 235 i=1,nwm
      if(iwsel(i).ne.0) then
      do 236 j=1,nwa
      xbmax=max(xbmax,w(i,j,1))
      xbmin=min(xbmin,w(i,j,1))
      ybmax=max(ybmax,w(i,j,2))
      ybmin=min(ybmin,w(i,j,2))
      zbmax=max(zbmax,w(i,j,3))
      zbmin=min(zbmin,w(i,j,3))
  236 continue
      endif
  235 continue
      endif
      if(scale.le.0.0d0) then
      scale=2.0d0/max(xbmax-xbmin,ybmax-ybmin,zbmax-zbmin)
      xbmax=scale*xbmax
      xbmin=scale*xbmin
      ybmax=scale*ybmax
      ybmin=scale*ybmin
      zbmax=scale*zbmax
      zbmin=scale*zbmin
      if(me.eq.0) write(*,3333) scale
 3333 format(' Scaling set to ',f12.6)
      endif
      call povinc(iunit,xbmin,xbmax,ybmin,ybmax,zbmin,zbmax)
      npov=npov+1
      write(fname,'(a,i5.5,a)') filpov(1:index(filpov,'.pov')-1),
     + npov,'.pov '
      open(unit=iunit,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown')
      write(iunit,8000)
 8000 format('#include "camera.inc"',/,'#include "colors.inc"',/,
     + '#include "plane.inc"')
      rv(1)=0.0d0
      rv(2)=0.0d0
      rv(3)=0.0d0
      rw(1)=0.0d0
      rw(2)=1.0d0
      rw(3)=0.0d0
      if(rangle.gt.0.0d0) rotang=(time/rangle)*twopi
      do 10 i=1,nsa
      rx(1)=x(i,1)
      rx(2)=x(i,2)
      rx(3)=x(i,3)
      if(rangle.gt.0.0d0) then
      ry(1)=x(i,1)
      ry(2)=x(i,2)
      ry(3)=x(i,3)
      call rotate(rv,rw,rotang,ry,rx)
      endif
      if(isel(i).ne.0) then
      write(iunit,8003) sgmnam(i)
 8003 format('// ',a)
      k=pre_atnum(sgmnam(i)(6:7))
      color=atom_color(k)
      factor=scale*val(i,1)
      if(tag(i,2)(1:1).ne.' ') color=tag(i,2)
      if(tag(i,2)(1:5).eq.'atom ') color=atom_color(k) 
      if(k.gt.0) then
      if(tag(i,1)(1:3).eq.'cpk') then
      write(iunit,8001) (scale*rx(j),j=1,3),
     + factor*abs(cpk)*ana_atrad(k),
     + color(1:index(color,' ')-1)
      else
      write(iunit,8001) (scale*rx(j),j=1,3),
     + factor*abs(stick)*ana_atrad(1),
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
      endif
   10 continue
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
      do 11 i=1,nwa+nsa+nwb
      read(44,4402) cdummy
   11 continue
      do 12 i=1,nsb
      read(44,6004) ib,jb
      if(isel(ib).gt.0.and.isel(jb).gt.0) then
      if(rangle.gt.0.0d0) then
      ry(1)=x(ib,1)
      ry(2)=x(ib,2)
      ry(3)=x(ib,3)
      call rotate(rv,rw,rotang,ry,rxi)
      ry(1)=x(jb,1)
      ry(2)=x(jb,2)
      ry(3)=x(jb,3)
      call rotate(rv,rw,rotang,ry,rxj)
      else
      rxi(1)=x(ib,1)
      rxi(2)=x(ib,2)
      rxi(3)=x(ib,3)
      rxj(1)=x(jb,1)
      rxj(2)=x(jb,2)
      rxj(3)=x(jb,3)
      endif
      xh(1)=0.5d0*(rxi(1)+rxj(1))
      xh(2)=0.5d0*(rxi(2)+rxj(2))
      xh(3)=0.5d0*(rxi(3)+rxj(3))
      if(tag(ib,1)(1:3).ne.'cpk') then
      k=pre_atnum(sgmnam(ib)(6:7))
      color=atom_color(k)
      factor=scale*val(ib,1)
      if(tag(ib,2)(1:1).ne.' ') color=tag(ib,2)
      if(tag(ib,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*rxi(j),j=1,3),
     + (scale*xh(j),j=1,3),factor*abs(stick)*ana_atrad(1),
     + color(1:index(color,' ')-1)
 8002 format('cylinder {',/,'  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + '  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + 3x,f12.6,/,'  open',/,'  texture { pigment { ',a,' } }',/,'}')
      endif
      if(tag(jb,1)(1:3).ne.'cpk') then
      k=pre_atnum(sgmnam(jb)(6:7))
      color=atom_color(k)
      factor=scale*val(jb,1)
      if(tag(jb,2)(1:1).ne.' ') color=tag(jb,2)
      if(tag(jb,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*rxj(j),j=1,3),
     + (scale*xh(j),j=1,3),factor*abs(stick)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
      endif
   12 continue
      if(lxw.and..not.lsonly) then
      do 13 i=1,nwm
      if(iwsel(i).ne.0) then
      do 14 j=1,nwa
      if(rangle.gt.0.0d0) then
      ry(1)=w(i,j,1)
      ry(2)=w(i,j,2)
      ry(3)=w(i,j,3)
      call rotate(rv,rw,rotang,ry,rx)
      else
      rx(1)=w(i,j,1)
      rx(2)=w(i,j,2)
      rx(3)=w(i,j,3)
      endif
      k=pre_atnum(wnam(j)(6:7))
      color=atom_color(k)
      factor=scale*wval(j,1)
      if(wtag(j,2)(1:1).ne.' ') color=wtag(j,2)
      if(wtag(j,2)(1:5).eq.'atom ') color=atom_color(k) 
      if(k.gt.0) then
      if(wtag(j,1)(1:3).eq.'cpk') then
      write(iunit,8001) (scale*rx(k),k=1,3),
     + factor*abs(cpk)*ana_atrad(k),
     + color(1:index(color,' ')-1)
      else
      write(iunit,8001) (scale*rx(k),k=1,3),
     + factor*abs(stick)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
      endif
c
   14 continue
      endif
   13 continue
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
      do 17 i=1,nwa+nsa
      read(44,4402) cdummy
   17 continue
      do 15 l=1,nwb
      read(44,6004) ib,jb
      if(pre_atnum(wnam(ib)(6:7)).gt.1.or.
     + pre_atnum(wnam(jb)(6:7)).gt.1.) then
      do 16 i=1,nwm
      if(iwsel(i).ne.0) then
      if(rangle.gt.0.0d0) then
      ry(1)=w(i,ib,1)
      ry(2)=w(i,ib,2)
      ry(3)=w(i,ib,3)
      call rotate(rv,rw,rotang,ry,rxi)
      ry(1)=w(i,jb,1)
      ry(2)=w(i,jb,2)
      ry(3)=w(i,jb,3)
      call rotate(rv,rw,rotang,ry,rxj)
      else
      rxi(1)=w(i,ib,1)
      rxi(2)=w(i,ib,2)
      rxi(3)=w(i,ib,3)
      rxj(1)=w(i,jb,1)
      rxj(2)=w(i,jb,2)
      rxj(3)=w(i,jb,3)
      endif
      xh(1)=0.5d0*(rxi(1)+rxj(1))
      xh(2)=0.5d0*(rxi(2)+rxj(2))
      xh(3)=0.5d0*(rxi(3)+rxj(3))
      if(wtag(ib,1)(1:3).ne.'cpk') then
      k=pre_atnum(wnam(ib)(6:7))
      color=atom_color(k)
      factor=scale*val(ib,1)
      if(wtag(ib,2)(1:1).ne.' ') color=wtag(ib,2)
      if(wtag(ib,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*rxi(j),j=1,3),
     + (scale*xh(j),j=1,3),factor*abs(stick)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
      if(wtag(jb,1)(1:3).ne.'cpk') then
      k=pre_atnum(wnam(jb)(6:7))
      color=atom_color(k)
      factor=scale*val(jb,1)
      if(wtag(jb,2)(1:1).ne.' ') color=wtag(jb,2)
      if(wtag(jb,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*rxj(j),j=1,3),
     + (scale*xh(j),j=1,3),factor*abs(stick)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
      endif
   16 continue
      endif
   15 continue
      endif
      close(unit=iunit)
      endif
c
      endif
c
      if(fmt.eq.'esp') then
      if(me.eq.0) then
      write(filchg,'(a,i5.5,a)') root(1:index(root,'.')-1),indx,'.xyz '
      write(filplt,'(a,i5.5,a)') root(1:index(root,'.')-1),indx,'.plt '
      write(filepz,'(a,i5.5,a)') root(1:index(root,'.')-1),indx,'.epz '
      call ana_wrtxyz(lfnchg,filchg,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.not.lsonly,int_mb(i_ndxw))
      write(filchg,'(a,i5.5,a)')
     + root(1:index(root,'.')-1),indx,'.xyzq '
      endif
      call ana_pltgrd(.true.,int_mb(i_isel))
      endif
c
      return
c
 9999 call md_abort('Oops',0)
      return
      end
