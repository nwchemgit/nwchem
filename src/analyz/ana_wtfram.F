      subroutine ana_wtfram(iunit,fmt,x,sgmnam,w,isel,tag,iwsel,ndxw)
c
c $Id: ana_wtfram.F,v 1.30 2000-08-16 20:22:14 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3)
      integer isel(msa),iwsel(mwm),ndxw(nwm)
      character*24 tag(msa,2)
      integer iunit
      character*3 fmt
      character*16 sgmnam(msa)
c
      character*10 today,now
      integer i,j,k
      character*2 elemnt
c
      character*1 cdummy
      character*7 string
      integer nsb,ib,jb,nss,nwb
c
      character*255 fname
      character*255 color
      real*8 xh(3)
c
      real*8 xbmax,xbmin,ybmax,ybmin,zbmax,zbmin
      integer*4 numato
c
      real*8 ana_atrad
      integer pre_atnum
      character*255 atom_color
      external ana_atrad,pre_atnum,atom_color
c
      if(me.eq.0) then
c
      if(iunit.eq.lfncop) then
      if(ncop.gt.mcopf.and.mcopf.gt.0) then
      ncop=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      icopf=icopf+1
      call ana_wthdr(lfncop,fmtcop,byte_mb(i_snam),tag,isel)
      endif
      endif
c
      if(iunit.eq.lfnsup) then
      if(nsup.gt.msupf.and.msupf.gt.0) then
      nsup=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      isupf=isupf+1
      call ana_wthdr(lfnsup,fmtsup,byte_mb(i_snam),tag,isel)
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
      write(iunit,1003) lxw.and..not.lsonly,.false.,lsx,.false.,
     + nwsel,nwa,nsa
 1003 format(4l1,3i10)
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
      endif
c
      if(fmt.eq.'arc') then
      call swatch(today,now)
      write(iunit,2000) time,today,now
 2000 format('title ',t67,f10.3,/,'!DATE ',2a10)
      nss=0
      if(lsx) then
      do 4 i=1,nsa
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
      write(fname,'(a,i3.3,a)') filpov(1:index(filpov,'.pov')-1),
     + npov,'.pov '
      open(unit=iunit,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown')
      write(iunit,8000)
 8000 format('#include "camera.inc"',/,'#include "colors.inc"')
      do 10 i=1,nsa
      if(isel(i).ne.0) then
      k=pre_atnum(sgmnam(i)(6:7))
      color=atom_color(k)
      if(tag(i,2)(1:1).ne.' ') color=tag(i,2)
      if(tag(i,2)(1:5).eq.'atom ') color=atom_color(k) 
      if(k.gt.0) then
      if(tag(i,1)(1:3).eq.'cpk') then
      write(iunit,8001) (scale*x(i,j),j=1,3),
     + scale*abs(cpk)*ana_atrad(k),
     + color(1:index(color,' ')-1)
      else
      write(iunit,8001) (scale*x(i,j),j=1,3),
     +scale*abs(cpk)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
 8001 format('sphere {',/,'  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + 3x,f12.6,/,'  texture { pigment { ',a,' } }',/,'}')
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
      xh(1)=0.5d0*(x(ib,1)+x(jb,1))
      xh(2)=0.5d0*(x(ib,2)+x(jb,2))
      xh(3)=0.5d0*(x(ib,3)+x(jb,3))
      if(tag(ib,1)(1:3).ne.'cpk') then
      k=pre_atnum(sgmnam(ib)(6:7))
      color=atom_color(k)
      if(tag(ib,2)(1:1).ne.' ') color=tag(ib,2)
      if(tag(ib,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*x(ib,j),j=1,3),
     + (scale*xh(j),j=1,3),scale*abs(cpk)*ana_atrad(1),
     + color(1:index(color,' ')-1)
 8002 format('cylinder {',/,'  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + '  <',f12.6,',',f12.6,',',f12.6,'>,',/,
     + 3x,f12.6,/,'  open',/,'  texture { pigment { ',a,' } }',/,'}')
      endif
      if(tag(jb,1)(1:3).ne.'cpk') then
      k=pre_atnum(sgmnam(jb)(6:7))
      color=atom_color(k)
      if(tag(jb,2)(1:1).ne.' ') color=tag(jb,2)
      if(tag(jb,2)(1:5).eq.'atom ') color=atom_color(k) 
      write(iunit,8002) (scale*x(jb,j),j=1,3),
     + (scale*xh(j),j=1,3),scale*abs(cpk)*ana_atrad(1),
     + color(1:index(color,' ')-1)
      endif
      endif
   12 continue
      close(unit=iunit)
      endif
c
      if(fmt.eq.'esp') then
      write(filchg,'(a,i3.3,a)') root(1:index(root,'.')-1),nwrit,'.xyz '
      write(filplt,'(a,i3.3,a)') root(1:index(root,'.')-1),nwrit,'.plt '
      call ana_wrtxyz(lfnchg,filchg,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.not.lsonly,int_mb(i_ndxw))
      write(filchg,'(a,i5.5,a)')
     + root(1:index(root,'.')-1),nwrit,'.xyzq '
      call ana_pltgrd(.true.)
      endif
c
      endif
c
      return
c
 9999 call errquit('Oops',0)
      return
      end
